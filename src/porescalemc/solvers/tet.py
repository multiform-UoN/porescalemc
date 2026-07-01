"""Tetrahedral finite-element diffusion solver adapter.

Two solve modes, both based on the same homogenisation cell problem
-------------------------------------------------------------------
periodic=True (default)
    Proper periodic cell problem.  Periodic node pairs are identified by
    coordinate matching on the Gmsh-periodic mesh; DOF condensation via a
    sparse restriction matrix R collapses slave nodes onto masters:

        K_r = R^T K R,    f_r^(j) = R^T f^(j).

    The corrector χ^(j) solves  K_r χ_r = −f_r^(j)  and is gauge-fixed
    by subtracting its mean.

periodic=False
    Neumann (no-flux) cell problem on the non-periodic mesh.  The same
    stiffness K and load f^(j) are used without condensation; since K has
    a 1-D null space (constants), the gauge is fixed identically by
    mean-subtraction after the solve.

Both modes use the homogenisation formula

    D_ij = φ_f δ_ij + (1/V) χ^(j) · f^(i),

where f_i^(j) = ∫_fluid ∂φ_i/∂x_j dΩ is assembled element-wise over the
P1 fluid tetrahedra and V is the total box volume.

The returned QoI vector is [D_xx, D_yy, D_zz, mesh_porosity, num_elements].
The full 3×3 tensor is stored in diagnostics under keys D_xy, D_yz, etc.
"""

from __future__ import annotations

import logging

import numpy as np
from scipy.sparse import coo_matrix, eye as speye
from scipy.sparse.linalg import spsolve

from porescalemc.geometry.grains import Packing
from porescalemc.solvers.base import SolverProtocol

_log = logging.getLogger(__name__)

try:
    import gmsh  # type: ignore
    from porescalemc.geometry.mesher import (
        GmshPackingMesher,
        _GMSH_AVAILABLE,
        _physical_group_tags,
    )
    _HAS_GMSH = bool(_GMSH_AVAILABLE)
except ImportError:
    _HAS_GMSH = False


TET_DIFFUSION_QOI_NAMES = (
    "diffusivity_x",
    "diffusivity_y",
    "diffusivity_z",
    "mesh_porosity",
    "num_elements",
)


class TetDiffusionSolver(SolverProtocol):
    """Gmsh tetrahedral P1-FEM effective diffusivity solver.

    Returns ``[D_xx, D_yy, D_zz, mesh_porosity, num_elements]``.

    With ``periodic=True`` (default) this solves the full periodic cell
    problem via DOF condensation; with ``periodic=False`` it uses pure
    Neumann BCs (no-flux on all boundaries) — both modes fix the gauge
    by mean-subtraction so the system is always well-posed.
    """

    qoi_names = TET_DIFFUSION_QOI_NAMES

    def __init__(
        self,
        packing: Packing | None = None,
        mesh_size: float = 0.05,
        periodic: bool = True,
        order: int = 1,
        algo_3d: int = 4,
        verbosity: int = 0,
    ):
        if not _HAS_GMSH:
            raise ImportError("TetDiffusionSolver requires the optional gmsh dependency.")
        if int(order) != 1:
            raise ValueError("TetDiffusionSolver currently supports only linear P1 tetrahedra.")
        self.packing = packing
        self.mesh_size = float(mesh_size)
        self.periodic = periodic
        self.order = int(order)
        self.algo_3d = int(algo_3d)
        self.verbosity = int(verbosity)
        self.mesher: GmshPackingMesher | None = None
        self._result: np.ndarray | None = None
        self._diagnostics: dict[str, float] = {}

    def setup(self, packing: Packing) -> None:
        self.packing = packing
        self.mesher = GmshPackingMesher(
            mesh_size=self.mesh_size,
            periodic=self.periodic,
            order=self.order,
            algo_3d=self.algo_3d,
            verbosity=self.verbosity,
        )
        self._result = None
        self._diagnostics = {}

    def solve(self) -> np.ndarray:
        if self.packing is None or self.mesher is None:
            raise RuntimeError("TetDiffusionSolver.setup must be called first.")
        self.mesher.build(self.packing)
        porosity = float(self.mesher.get_mesh_porosity())
        n_elements = float(self.mesher.get_num_elements())

        D_tensor = _solve_cell_problem(self.packing.box, use_periodic=self.periodic)
        diffusivities = [float(D_tensor[i, i]) for i in range(3)]

        self._result = np.array([*diffusivities, porosity, n_elements], dtype=float)
        self._diagnostics = {
            "diffusivity_x": diffusivities[0],
            "diffusivity_y": diffusivities[1],
            "diffusivity_z": diffusivities[2],
            "mesh_porosity": porosity,
            "num_elements": n_elements,
            "mesh_size": self.mesh_size,
            "mode": "periodic" if self.periodic else "neumann",
        }
        for i in range(3):
            for j in range(3):
                self._diagnostics[f"D_{'xyz'[i]}{'xyz'[j]}"] = float(D_tensor[i, j])
        return self._result

    def diagnostics(self) -> dict[str, float]:
        return dict(self._diagnostics)

    def close(self) -> None:
        if self.mesher is not None:
            self.mesher.finalize()
            self.mesher = None


# ---------------------------------------------------------------------------
# Shared cell-problem solver (periodic or Neumann)
# ---------------------------------------------------------------------------

def _solve_cell_problem(
    box: tuple[float, float, float],
    use_periodic: bool = True,
) -> np.ndarray:
    """Assemble and solve the homogenisation cell problem; return 3×3 D_eff.

    Both modes:
    1. Assemble K (P1 stiffness on fluid tets) and F (load, shape n×3).
    2. Optionally condense periodic DOFs via restriction matrix R.
    3. Regularise the null space with eps·I so the system is non-singular.
    4. Solve K_sys χ = −F_sys for each axis; gauge-fix with mean-subtraction.
    5. Compute D_ij = φ_f·δ_ij + (1/V) χ^(j)·F^(i).
    """
    K, F, node_coords, fluid_vol = _assemble_fluid_p1_system()
    Lx, Ly, Lz = box
    V = Lx * Ly * Lz
    phi_f = fluid_vol / V
    n_full = len(node_coords)

    if use_periodic:
        pairs = _find_periodic_node_pairs(node_coords, box)
        dof_condensed, n_reduced = _condense_dofs(n_full, pairs)
        R = coo_matrix(
            (np.ones(n_full, dtype=float),
             (np.arange(n_full), dof_condensed)),
            shape=(n_full, n_reduced),
        ).tocsr()
        K_sys = (R.T @ K @ R).tocsr()
        F_sys = (R.T @ F)          # (n_reduced, 3)
        mode = "periodic"
    else:
        R = None
        K_sys = K.tocsr()
        F_sys = F                  # (n_full, 3)
        mode = "neumann"

    # Null-space regularisation: K has exactly one zero eigenvalue (constants).
    # eps·I shifts it just enough to make spsolve converge; mean-subtraction
    # afterwards removes the gauge ambiguity from the solution.
    n_sys = K_sys.shape[0]
    eps = float(K_sys.diagonal().max()) * 1e-10
    K_reg = K_sys + speye(n_sys, format="csr") * eps

    chi_full = np.zeros((n_full, 3), dtype=float)
    for j in range(3):
        chi_r = spsolve(K_reg, -F_sys[:, j])
        chi_r -= chi_r.mean()                          # fix gauge
        chi_full[:, j] = (R @ chi_r) if R is not None else chi_r

    D_eff = np.zeros((3, 3), dtype=float)
    for i in range(3):
        for j in range(3):
            D_eff[i, j] = phi_f * float(i == j) + float(chi_full[:, j] @ F[:, i]) / V

    _log.debug(
        "FEM cell problem (%s): phi_f=%.3f, D_eff diag=[%.4f, %.4f, %.4f]",
        mode, phi_f, D_eff[0, 0], D_eff[1, 1], D_eff[2, 2],
    )
    return D_eff


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

def _assemble_fluid_p1_system() -> tuple[object, np.ndarray, np.ndarray, float]:
    """Assemble K, load matrix F, node coords, and fluid volume.

    Returns
    -------
    K : csr_matrix  (n × n)  — P1 Laplace stiffness on fluid tets
    F : ndarray     (n × 3)  — load vectors; F[i,j] = ∫_fluid ∂φ_i/∂x_j dΩ
    node_coords : ndarray (n × 3)
    fluid_volume : float
    """
    node_tags, coord_flat, _ = gmsh.model.mesh.getNodes()
    all_coords = np.asarray(coord_flat, dtype=float).reshape(-1, 3)
    tag_to_global = {int(t): i for i, t in enumerate(node_tags)}

    fluid_entities = _physical_group_tags(3, "fluid")
    if not fluid_entities:
        raise RuntimeError("No 'fluid' physical group found in the Gmsh model.")

    all_tet_nodes_list: list[np.ndarray] = []
    for ent in fluid_entities:
        etypes, _, enode_tags = gmsh.model.mesh.getElements(3, ent)
        for etype, flat in zip(etypes, enode_tags):
            if int(etype) == 4:  # linear tet
                arr = np.asarray(flat, dtype=np.int64).reshape(-1, 4)
                all_tet_nodes_list.append(arr)

    if not all_tet_nodes_list:
        raise RuntimeError("No linear tetrahedra found in the fluid physical group.")

    all_tet_nodes = np.vstack(all_tet_nodes_list)  # (n_tets, 4)
    used_tags = sorted(set(all_tet_nodes.ravel().tolist()))
    compact = {tag: i for i, tag in enumerate(used_tags)}
    node_coords = np.array([all_coords[tag_to_global[t]] for t in used_tags], dtype=float)
    n = len(used_tags)

    rows: list[int] = []
    cols: list[int] = []
    K_data: list[float] = []
    F = np.zeros((n, 3), dtype=float)
    fluid_volume = 0.0

    for tet_tags in all_tet_nodes:
        local = [compact[int(t)] for t in tet_tags]
        x = node_coords[local]                        # (4, 3)
        A_mat = np.column_stack((np.ones(4), x))      # (4, 4)
        try:
            Ainv = np.linalg.inv(A_mat)
        except np.linalg.LinAlgError:
            continue
        detA = float(np.linalg.det(A_mat))
        volume = abs(detA) / 6.0
        if volume <= 0.0:
            continue

        # grads[i] = grad(φ_i) = Ainv[1:, i]  →  shape (4, 3)
        grads = Ainv[1:, :].T
        ke = volume * (grads @ grads.T)   # (4, 4) element stiffness
        fe = volume * grads               # (4, 3): fe[i,j] = ∫ ∂φ_i/∂x_j

        fluid_volume += volume
        for a, ga in enumerate(local):
            F[ga] += fe[a]
            for b, gb in enumerate(local):
                rows.append(ga)
                cols.append(gb)
                K_data.append(float(ke[a, b]))

    K = coo_matrix((K_data, (rows, cols)), shape=(n, n)).tocsr()
    return K, F, node_coords, fluid_volume


# ---------------------------------------------------------------------------
# Periodic DOF condensation helpers
# ---------------------------------------------------------------------------

def _find_periodic_node_pairs(
    node_coords: np.ndarray,
    box: tuple[float, float, float],
    tol_rel: float = 1e-6,
) -> list[tuple[int, int]]:
    """Return (slave_idx, master_idx) pairs for all three face-pair directions.

    Slave nodes sit on the positive faces (x+, y+, z+); masters are the
    matching nodes on the negative faces.  Chains (corner/edge nodes on
    multiple positive faces) are resolved so all slaves ultimately map to
    a single negative-face master.
    """
    n = len(node_coords)
    dof = np.arange(n, dtype=np.intp)

    for axis, L in enumerate(box):
        tol = tol_rel * max(L, 1.0)
        lo_idx = np.where(node_coords[:, axis] <= tol)[0]
        hi_idx = np.where(node_coords[:, axis] >= L - tol)[0]
        if lo_idx.size == 0 or hi_idx.size == 0:
            continue

        other = [a for a in range(3) if a != axis]
        lo_trans = node_coords[lo_idx][:, other]  # (n_lo, 2)

        for i in hi_idx:
            y_i = node_coords[i, other[0]]
            z_i = node_coords[i, other[1]]
            dist = np.abs(lo_trans[:, 0] - y_i) + np.abs(lo_trans[:, 1] - z_i)
            j_local = int(np.argmin(dist))
            if dist[j_local] < tol:
                dof[i] = lo_idx[j_local]

    # Resolve chains: at most 3 hops in 3-D
    for _ in range(3):
        dof = dof[dof]

    return [(int(i), int(dof[i])) for i in range(n) if dof[i] != i]


def _condense_dofs(
    n_full: int,
    pairs: list[tuple[int, int]],
) -> tuple[np.ndarray, int]:
    """Map full node indices to reduced DOF indices, collapsing slave→master.

    Returns ``(dof_condensed, n_reduced)`` where ``dof_condensed[i]`` is the
    condensed index for node ``i``.
    """
    canon = np.arange(n_full, dtype=np.intp)
    for slave, master in pairs:
        canon[slave] = master
    for _ in range(3):
        canon = canon[canon]

    unique_vals = sorted(set(canon.tolist()))
    remap = {old: new for new, old in enumerate(unique_vals)}
    dof_condensed = np.array([remap[int(c)] for c in canon], dtype=np.intp)
    return dof_condensed, len(unique_vals)
