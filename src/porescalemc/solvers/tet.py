"""Tetrahedral finite-element diffusion solver adapter.

The solver builds a Gmsh tetrahedral mesh of the pore space and assembles a
linear P1 finite-element Laplace operator on the fluid tetrahedra.  For each
coordinate direction it solves a unit Dirichlet conduction test between the two
opposite box faces, with natural no-flux conditions on the remaining outer
faces and grain walls.  The resulting energy gives a directional effective
diffusivity estimate.

This is a useful body-fitted FEM validation path for the pure-Python MLMC
plumbing.  It is still simpler than a fully periodic homogenization FEM solve:
periodic node constraints and cell-problem fluctuations are left for a later,
more careful implementation.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from porescalemc.geometry.grains import Packing
from porescalemc.solvers.base import SolverProtocol

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
    """Gmsh tetrahedral P1-FEM directional diffusion estimate.

    The returned vector is
    ``[D_eff_x, D_eff_y, D_eff_z, mesh_porosity, num_elements]``.
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
        stiffness, coords = _assemble_fluid_p1_stiffness()
        diffusivities = [
            _solve_dirichlet_effective_diffusivity(stiffness, coords, self.packing.box, axis)
            for axis in range(3)
        ]
        self._result = np.array([*diffusivities, porosity, n_elements], dtype=float)
        self._diagnostics = {
            "diffusivity_x": float(diffusivities[0]),
            "diffusivity_y": float(diffusivities[1]),
            "diffusivity_z": float(diffusivities[2]),
            "mesh_porosity": porosity,
            "num_elements": n_elements,
            "mesh_size": self.mesh_size,
        }
        return self._result

    def diagnostics(self) -> dict[str, float]:
        return dict(self._diagnostics)

    def close(self) -> None:
        if self.mesher is not None:
            self.mesher.finalize()
            self.mesher = None


def _assemble_fluid_p1_stiffness() -> tuple[coo_matrix, np.ndarray]:
    """Assemble the scalar P1 Laplace stiffness matrix on fluid tetrahedra."""
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
            if int(element_type) != 4:  # 4-node linear tetrahedron
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
    """Solve one unit-gradient Dirichlet conduction problem and return D_eff."""
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
