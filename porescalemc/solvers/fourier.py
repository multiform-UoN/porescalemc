"""Solver that uses the analytical Fourier field for QoIs.

Provides access to:
- porosity computed from the smoothed solid-fraction field
- upscaled permeability (harmonic mean via chosen law)

This gives a natural discretization level ('g' in hierarchy) that is
independent of the geometric packing resolution.
"""

from __future__ import annotations

import numpy as np

from porescalemc.config import FourierConfig
from porescalemc.geometry.fourier_field import (
    porosity_field,
    upscaled_permeability,
)
from porescalemc.geometry.grains import Packing
from porescalemc.solvers.base import SolverProtocol


class FourierSolver(SolverProtocol):
    """Returns Fourier-derived QoIs for a packing.

    By default returns [mean_solid_fraction_from_field, upscaled_permeability].
    """

    def __init__(
        self,
        packing: Packing | None = None,
        fourier_config: FourierConfig | None = None,
        law: str = "kozeny_carman",
    ):
        self.packing = packing
        self.fourier_config = fourier_config or FourierConfig()
        self.law = law
        self._result: np.ndarray | None = None

    def setup(self, packing: Packing) -> None:
        self.packing = packing
        self._result = None

    def solve(self) -> np.ndarray:
        if self.packing is None:
            raise RuntimeError("FourierSolver.setup must be called first.")

        solid = porosity_field(self.packing, self.fourier_config)
        mean_solid = float(np.mean(solid))

        # Upscaled permeability (scalar for now)
        perm = upscaled_permeability(
            self.packing, self.fourier_config, law=self.law
        )

        self._result = np.array([mean_solid, perm], dtype=float)
        return self._result

    def close(self) -> None:
        pass
