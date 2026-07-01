"""Tetrahedral finite-element diffusion solver adapter.

Two solve modes
---------------
periodic=False
    Unit-gradient Dirichlet conduction test between opposite faces.
    Fast and simple; does not exploit the periodic cell structure.

periodic=True (default)
    Proper periodic homogenisation cell problem.

    After building a periodic Gmsh mesh, periodic node pairs on opposite
    faces are identified by coordinate matching.  DOF condensation via a
    sparse restriction matrix R (n_full × n_reduced) collapses slave nodes
    onto their masters, giving the reduced system

        K_r  =  R^T K R,    f_r^(j) = R^T f^(j)

    where K is the standard P1 Laplace stiffness assembled on fluid tets
    and f_i^(j) = ∫_fluid (∂φ_i/∂x_j) dΩ is the homogenisation load.

    The corrector χ^(j) solves  K_r χ_r = −f_r^(j)  (one solve per axis),
    and the effective diffusivity tensor follows from

        D_ij = φ_f δ_ij + (1/V) χ^(j) · f^(i).

The returned QoI vector is [D_xx, D_yy, D_zz, mesh_porosity, num_elements].
The full 3×3 tensor is stored in diagnostics as "D_tensor".
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
    problem; with ``periodic=False`` it falls back to the simpler
    unit-gradient Dirichlet test between opposite faces.
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

        if self.periodic:
            D_tensor = _solve_periodic_diffusion(self.packing.box)
            diffusivities = [float(D_tensor[i, i]) for i in range(3)]
            self._diagnostics = {
                "diffusivity_x": diffusivities[0],
                "diffusivity_y": diffusivities[1],
                "diffusivity_z": diffusivities[2],
                "mesh_porosity": porosity,
                "num_elements": n_elements,
                "mesh_size": self.mesh_size,
                "mode": "periodic",
            }
            for i in range(3):
                for j in range(3):
                    self._diagnostics[f"D_{'xyz'[i]}{'xyz'[j]}"] = float(D_tensor[i, j])
        else:
            stiffness, coords = _assemble_fluid_p1_stiffness()
            diffusivities = [
                _solve_dirichlet_effective_diffusivity(stiffness, coords, self.packing.box, axis)
                for axis in range(3)
            ]
            self._diagnostics = {
                "diffusivity_x": float(diffusivities[0]),
                "diffusivity_y": float(diffusivities[1]),
                "diffusivity_z": float(diffusivities[2]),
                "mesh_porosity": porosity,
                "num_elements": n_elements,
                "mesh_size": self.mesh_size,
                "mode": "dirichlet",
            }

        self._result = np.array([*diffusivities, porosity, n_elements], dtype=float)
        return self._result

    def diagnostics(self) -> dict[str, float]:
        return dict(self._diagnostics)

    def close(self) -> None:
        if self.mesher is not None:
            self.mesher.finalize()
            self.mesher = None


# ---------------------------------------------------------------------------
# Periodic cell problem
# ---------------------------------------------------------------------------

def _solve_periodic_diffusion(box: tuple[float, float, float]) -> np.ndarray:
    """Solve the periodic FEM cell problem and return the 3×3 D_eff tensor."""
    K, F, node_coords, fluid_vol = _assemble_fluid_p1_system()
    Lx, Ly, Lz = box
    V = Lx * Ly * Lz
    phi_f = fluid_vol / V

    pairs = _find_periodic_node_pairs(node_coords, box)
    n_full = len(node_coords)
    dof_condensed, n_reduced = _condense_dofs(n_full, pairs)

    # Build restriction matrix R: n_full × n_reduced (each row has one 1)
    R = coo_matrix(
        (np.ones(n_full, dtype=float),
         (np.arange(n_full), dof_condensed)),
        shape=(n_full, n_reduced),
    ).tocsr()

    # Reduced stiffness
    K_r = (R.T @ K @ R).tocsr()

    # Regularise null space (constant vector): eps * I
    eps = float(K_r.diagonal().max()) * 1e-10
    K_r = K_r + speye(n_reduced, format="csr") * eps

    # Reduced load vectors (n_reduced × 3)
    F_r = (R.T @ F)

    D_eff = np.zeros((3, 3), dtype=float)
    chi_full = np.zeros((n_full, 3), dtype=float)

    for j in range(3):
        chi_r = spsolve(K_r, -F_r[:, j])
        chi_r -= chi_r.mean()
        chi_full[:, j] = R @ chi_r

    for i in range(3):
        for j in range(3):
            D_eff[i, j] = phi_f * float(i == j) + float(chi_full[:, j] @ F[:, i]) / V

    _log.debug(
        "Periodic FEM: phi_f=%.3f, D_eff diag=[%.4f, %.4f, %.4f]",
        phi_f, D_eff[0, 0], D_eff[1, 1], D_eff[2, 2],
    )
    return D_eff


def _assemble_fluid_p1_system() -> tuple[object, np.ndarray, np.ndarray, float]:
    """Assemble K, load matrix F, node coords, and fluid volume from current Gmsh model.

    Returns
    -------
    K : csr_matrix  (n × n)
    F : ndarray     (n × 3) — load vectors for all three axes
    node_coords : ndarray  (n × 3)
    fluid_volume : float
    """
    node_tags, coord_flat, _ = gmsh.model.mesh.getNodes()
    all_coords = np.asarray(coord_flat, dtype=float).reshape(-1, 3)
    tag_to_global = {int(t): i for i, t in enumerate(node_tags)}

    fluid_entities = _physical_group_tags(3, "fluid")
    if not fluid_entities:
        raise RuntimeError("No 'fluid' physical group found in the Gmsh model.")

    # Collect all fluid linear tetrahedra (element type 4)
    all_tet_nodes_list: list[np.ndarray] = []
    for ent in fluid_entities:
        etypes, _, enode_tags = gmsh.model.mesh.getElements(3, ent)
        for etype, flat in zip(etypes, enode_tags):
            if int(etype) == 4:
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
        x = node_coords[local]  # (4, 3)
        A_mat = np.column_stack((np.ones(4), x))  # (4, 4)
        try:
            Ainv = np.linalg.inv(A_mat)
        except np.linalg.LinAlgError:
            continue
        detA = float(np.linalg.det(A_mat))
        volume = abs(detA) / 6.0
        if volume <= 0.0:
            continue

        # grads[i] = grad(φ_i) = Ainv[1:, i]  →  grads shape (4, 3)
        grads = Ainv[1:, :].T
        ke = volume * (grads @ grads.T)  # (4, 4)
        fe = volume * grads              # (4, 3): fe[i,j] = ∫ ∂φ_i/∂x_j

        fluid_volume += volume
        for a, ga in enumerate(local):
            F[ga] += fe[a]
            for b, gb in enumerate(local):
                rows.append(ga)
                cols.append(gb)
                K_data.append(float(ke[a, b]))

    K = coo_matrix((K_data, (rows, cols)), shape=(n, n)).tocsr()
    return K, F, node_coords, fluid_volume


def _find_periodic_node_pairs(
    node_coords: np.ndarray,
    box: tuple[float, float, float],
    tol_rel: float = 1e-6,
) -> list[tuple[int, int]]:
    """Return (slave_idx, master_idx) pairs for all three periodic face directions.

    Slave nodes are those on the positive faces (x+, y+, z+); their masters
    are on the corresponding negative faces with the same transverse coords.
    Corner/edge nodes are identified with the negative-face master for every
    axis on which they sit on the positive face, and chains are resolved at
    the end so all slaves eventually point to the corner master.
    """
    n = len(node_coords)
    # dof[i] = canonical master index for node i; initially identity
    dof = np.arange(n, dtype=np.intp)
    Lx, Ly, Lz = box

    for axis, L in enumerate((Lx, Ly, Lz)):
        tol = tol_rel * max(L, 1.0)
        lo_mask = node_coords[:, axis] <= tol
        hi_mask = node_coords[:, axis] >= L - tol

        lo_idx = np.where(lo_mask)[0]
        hi_idx = np.where(hi_mask)[0]
        if lo_idx.size == 0 or hi_idx.size == 0:
            continue

        other = [a for a in range(3) if a != axis]
        lo_transverse = node_coords[lo_idx][:, other]  # (n_lo, 2)

        for i in hi_idx:
            y_i = node_coords[i, other[0]]
            z_i = node_coords[i, other[1]]
            dist = np.abs(lo_transverse[:, 0] - y_i) + np.abs(lo_transverse[:, 1] - z_i)
            j_local = int(np.argmin(dist))
            if dist[j_local] < tol:
                dof[i] = lo_idx[j_local]

    # Resolve chains (e.g. edge or corner nodes): repeat 3 times is enough
    for _ in range(3):
        dof = dof[dof]

    pairs = [(int(i), int(dof[i])) for i in range(n) if dof[i] != i]
    return pairs


def _condense_dofs(
    n_full: int,
    pairs: list[tuple[int, int]],
) -> tuple[np.ndarray, int]:
    """Map full DOF indices to reduced DOF indices by collapsing slave→master.

    Returns
    -------
    dof_condensed : ndarray (n_full,) — reduced DOF index for each node
    n_reduced : int
    """
    canon = np.arange(n_full, dtype=np.intp)
    for slave, master in pairs:
        canon[slave] = master

    # Resolve any remaining chains
    for _ in range(3):
        canon = canon[canon]

    unique_vals = sorted(set(canon.tolist()))
    remap = {old: new for new, old in enumerate(unique_vals)}
    dof_condensed = np.array([remap[int(c)] for c in canon], dtype=np.intp)
    return dof_condensed, len(unique_vals)


# ---------------------------------------------------------------------------
# Dirichlet fallback
# ---------------------------------------------------------------------------

def _assemble_fluid_p1_stiffness() -> tuple[object, np.ndarray]:
    """Assemble P1 Laplace stiffness on fluid tets only.  Dirichlet path."""
    node_tags, coord_flat, _ = gmsh.model.mesh.getNodes()
    all_coords = np.asarray(coord_flat, dtype=float).reshape(-1, 3)
    tag_to_global = {int(tag): i for i, tag in enumerate(node_tags)}

    fluid_entities = _physical_group_tags(3, "fluid")
    if not fluid_entities:
        raise RuntimeError("No fluid physical group found in the Gmsh model.")

    tets: list[list[int]] = []
    used_tags: set[int] = set()
    for entity_tag in fluid_entities:
        element_types, _, element_node_tags = gmsh.model.mesh.getElements(3, entity_tag)
        for element_type, flat_tags in zip(element_types, element_node_tags):
            if int(element_type) != 4:
                continue
            tags = np.asarray(flat_tags, dtype=np.int64).reshape(-1, 4)
            for tet_tags in tags:
                tet = [int(tag) for tag in tet_tags]
                tets.append(tet)
                used_tags.update(tet)

    if not tets:
        raise RuntimeError("No linear tetrahedra found in the fluid physical group.")

    used_sorted = sorted(used_tags)
    compact = {tag: i for i, tag in enumerate(used_sorted)}
    coords = np.array([all_coords[tag_to_global[tag]] for tag in used_sorted], dtype=float)

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    for tet_tags in tets:
        local = [compact[tag] for tag in tet_tags]
        x = coords[local]
        A = np.column_stack((np.ones(4), x))
        detA = float(np.linalg.det(A))
        volume = abs(detA) / 6.0
        if volume <= 0.0:
            continue
        coeffs = np.linalg.inv(A)
        grads = coeffs[1:, :].T
        ke = volume * (grads @ grads.T)
        for i, gi in enumerate(local):
            for j, gj in enumerate(local):
                rows.append(gi)
                cols.append(gj)
                data.append(float(ke[i, j]))

    stiffness = coo_matrix((data, (rows, cols)), shape=(len(used_sorted), len(used_sorted))).tocsr()
    return stiffness, coords


def _solve_dirichlet_effective_diffusivity(
    stiffness,
    coords: np.ndarray,
    box: tuple[float, float, float],
    axis: int,
) -> float:
    """Solve one unit-gradient Dirichlet conduction problem; return D_eff."""
    lengths = np.asarray(box, dtype=float)
    length = float(lengths[axis])
    area = float(np.prod(np.delete(lengths, axis)))
    x = coords[:, axis]
    tol = max(1e-9, 1e-7 * length)
    lo = x <= tol
    hi = x >= length - tol
    fixed = lo | hi
    if not np.any(lo) or not np.any(hi):
        raise RuntimeError(f"Missing Dirichlet nodes on axis {axis}.")

    values = np.zeros(coords.shape[0], dtype=float)
    values[hi] = 1.0
    free = ~fixed
    if np.any(free):
        Kff = stiffness[free][:, free]
        Kfd = stiffness[free][:, fixed]
        rhs = -Kfd @ values[fixed]
        values[free] = spsolve(Kff, rhs)

    energy = float(values @ (stiffness @ values))
    diffusivity = energy * length / area
    return max(diffusivity, 0.0)
