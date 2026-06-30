"""Tetrahedral mesh-backed diffusion solver adapter.

This module intentionally keeps the current tetrahedral path modest: it builds
the Gmsh mesh, reports mesh QoIs, and returns a Bruggeman-style diffusivity
estimate from mesh porosity.  That is useful for MLMC hierarchy plumbing and
mesh convergence studies before a verified P1 finite-element cell-problem
assembly is added.
"""

from __future__ import annotations

import numpy as np

from porescalemc.geometry.grains import Packing
from porescalemc.solvers.base import SolverProtocol

try:
    from porescalemc.geometry.mesher import GmshPackingMesher, _GMSH_AVAILABLE
    _HAS_GMSH = bool(_GMSH_AVAILABLE)
except ImportError:
    _HAS_GMSH = False


TET_DIFFUSION_QOI_NAMES = (
    "diffusivity_estimate",
    "mesh_porosity",
    "num_elements",
)


class TetDiffusionSolver(SolverProtocol):
    """Gmsh tetrahedral mesh-backed diffusivity estimate.

    The returned vector is ``[D_eff_estimate, mesh_porosity, num_elements]``.
    ``D_eff_estimate = porosity**bruggeman_exponent`` by default, a standard
    scalar closure used here as an explicit placeholder until full FEM assembly
    is implemented.
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
        bruggeman_exponent: float = 1.5,
    ):
        if not _HAS_GMSH:
            raise ImportError("TetDiffusionSolver requires the optional gmsh dependency.")
        self.packing = packing
        self.mesh_size = float(mesh_size)
        self.periodic = periodic
        self.order = int(order)
        self.algo_3d = int(algo_3d)
        self.verbosity = int(verbosity)
        self.bruggeman_exponent = float(bruggeman_exponent)
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
        diffusivity = max(porosity, 0.0) ** self.bruggeman_exponent
        self._result = np.array([diffusivity, porosity, n_elements], dtype=float)
        self._diagnostics = {
            "diffusivity_estimate": float(diffusivity),
            "mesh_porosity": porosity,
            "num_elements": n_elements,
            "mesh_size": self.mesh_size,
            "bruggeman_exponent": self.bruggeman_exponent,
        }
        return self._result

    def diagnostics(self) -> dict[str, float]:
        return dict(self._diagnostics)

    def close(self) -> None:
        if self.mesher is not None:
            self.mesher.finalize()
            self.mesher = None
