"""Gmsh-based mesher solver for MLMC (soft dependency).

Provides a SolverProtocol implementation that builds a mesh from the packing
using GmshPackingMesher and returns simple mesh-derived QoIs.

Intended for future use with discrete mesh-based solvers. Currently returns
mesh porosity and element count (scaled) as demo QoI vector.

This allows MLMC to run over different mesh fidelities (via level-dependent
mesh_size) even before full PDE-on-mesh solvers are plugged in.
"""

from __future__ import annotations

import numpy as np

from porescalemc.geometry.grains import Packing
from porescalemc.solvers.base import SolverProtocol

try:
    from porescalemc.geometry.mesher import GmshPackingMesher, _GMSH_AVAILABLE
    _HAS_MESHER = bool(_GMSH_AVAILABLE)
except ImportError:
    _HAS_MESHER = False


class GmshMesherQoISolver(SolverProtocol):
    """Solver that uses Gmsh mesher to produce mesh-based QoIs.

    QoI vector: [mesh_porosity, num_elements / 1000.0 ]  (simple scalars for demo)

    Parameters
    ----------
    mesh_size : float
        Target mesh size. Can be level-dependent via factory.
    periodic : bool
        Pass to mesher.
    verbosity : int
    order : int
        Gmsh element order passed to the mesher.
    algo_3d : int
        Gmsh 3-D meshing algorithm tag.
    """

    def __init__(
        self,
        mesh_size: float = 0.05,
        periodic: bool = True,
        verbosity: int = 0,
        order: int = 1,
        algo_3d: int = 4,
        packing: Packing | None = None,
    ):
        if not _HAS_MESHER:
            raise ImportError("GmshPackingMesher not available (gmsh not installed)")
        self.mesh_size = mesh_size
        self.periodic = periodic
        self.verbosity = verbosity
        self.order = order
        self.algo_3d = algo_3d
        self.packing = packing
        self.mesher: GmshPackingMesher | None = None
        self._result: np.ndarray | None = None

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

    def solve(self) -> np.ndarray:
        if self.mesher is None or self.packing is None:
            raise RuntimeError("setup(packing) must be called first.")
        self.mesher.build(self.packing)
        porosity = self.mesher.get_mesh_porosity()
        n_elems = self.mesher.get_num_elements()
        # Simple vector QoI for MLMC demo
        self._result = np.array([porosity, float(n_elems) / 1000.0], dtype=float)
        return self._result

    def close(self) -> None:
        if self.mesher is not None:
            self.mesher.finalize()
            self.mesher = None
