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
]
