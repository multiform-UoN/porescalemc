"""Common solver protocol.

The MLMC estimator only needs a narrow interface: prepare a packing, return a
numeric result, and release resources.  External PDE solvers can implement this
later without changing the estimator.
"""

from __future__ import annotations

from typing import Protocol

import numpy as np

from porescalemc.geometry.grains import Packing


class SolverProtocol(Protocol):
    """Minimal solver interface consumed by :class:`MLMCEstimator`."""

    def setup(self, packing: Packing) -> None:
        """Prepare the solver for ``packing``."""

    def solve(self) -> np.ndarray:
        """Return a one-dimensional numeric result."""

    def close(self) -> None:
        """Release any resources held by the solver."""
