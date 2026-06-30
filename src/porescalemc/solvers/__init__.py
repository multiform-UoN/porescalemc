"""Pure-Python solver adapters used by the rewritten core package."""

from __future__ import annotations

from porescalemc.solvers.base import SolverProtocol
from porescalemc.solvers.dummy import DummySolver
from porescalemc.solvers.fourier import FourierSolver
from porescalemc.solvers.packing import (
    PACKING_QOI_NAMES,
    PackingStatsSolver,
    packing_statistics,
)
from porescalemc.solvers.spectral import (
    SpectralAdvectionDiffusionSolver,
    SpectralDiffusionSolver,
    SpectralStokesSolver,
)
from porescalemc.solvers.tet import TET_DIFFUSION_QOI_NAMES, TetDiffusionSolver
from porescalemc.solvers.voxel import VOXEL_DIFFUSION_QOI_NAMES, VoxelDiffusionSolver
from porescalemc.solvers.registry import (
    SOLVER_REGISTRY,
    SolverInfo,
    available_qois,
    available_solvers,
    solver_class,
)

__all__ = [
    "SolverProtocol",
    "DummySolver",
    "FourierSolver",
    "PACKING_QOI_NAMES",
    "PackingStatsSolver",
    "packing_statistics",
    "SpectralDiffusionSolver",
    "SpectralStokesSolver",
    "SpectralAdvectionDiffusionSolver",
    "VoxelDiffusionSolver",
    "VOXEL_DIFFUSION_QOI_NAMES",
    "TetDiffusionSolver",
    "TET_DIFFUSION_QOI_NAMES",
    "SOLVER_REGISTRY",
    "SolverInfo",
    "available_solvers",
    "available_qois",
    "solver_class",
]

# Gmsh mesher QoI solver is a soft dependency: importing the class is cheap,
# construction raises ImportError if the gmsh Python module is unavailable.
try:
    from porescalemc.solvers.gmsh_mesher import GmshMesherQoISolver
    __all__.append("GmshMesherQoISolver")
except ImportError:
    pass
