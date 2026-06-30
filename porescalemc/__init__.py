"""
porescalemc — Multilevel Monte Carlo for pore-scale UQ on random packings.

Pure-Python core for generating hierarchical sphere/ellipsoid packings
and running MLMC on cheap algebraic or Fourier-derived QoIs.

Key modules:
- config: PackingConfig, MLMCConfig, FourierConfig
- geometry: sample_packing_from_name, Packing, Grain, Fourier field
- solvers: PackingStatsSolver, SolverProtocol
- mlmc: MLMCEstimator (with deterministic paired sampling)

Quick start (packing statistics):

    from porescalemc import run_mlmc
    from porescalemc.config import PackingConfig, MLMCConfig

    pcfg = PackingConfig(mu=0.08, n_grains=25)
    mcfg = MLMCConfig(n_levels=3, tolerance=0.05)

    mean, err = run_mlmc(pcfg, mcfg)
    print(mean, err)

Basic logging setup:

    from porescalemc import setup_logging
    setup_logging("INFO")
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from porescalemc.config import (
    PackingConfig,
    MLMCConfig,
    FourierConfig,
    SpectralConfig,
    spectral_hierarchy_level,
)
from porescalemc.geometry.placement import sample_packing_from_name
from porescalemc.solvers.fourier import FourierSolver
from porescalemc.solvers.packing import PackingStatsSolver
from porescalemc.mlmc.estimator import MLMCEstimator

if TYPE_CHECKING:
    from numpy.typing import NDArray

__version__ = "0.1.0"
__author__ = "Matteo Icardi"


def setup_logging(level: str = "INFO") -> None:
    """Convenience one-liner for basic logging during development."""
    import logging
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


def run_mlmc(
    packing_config: PackingConfig,
    mlmc_config: MLMCConfig,
    *,
    qoi_names: tuple[str, ...] | None = None,
    n_workers: int = 1,
    algorithm: str | None = None,
) -> tuple["np.ndarray", "np.ndarray"]:
    """High-level convenience wrapper for the common packing-stats case.

    Uses deterministic name-based sampling so coarse/fine pairs share the
    same realization (essential for MLMC variance reduction).

    Returns (mean, error_bound).
    """
    import numpy as np

    from porescalemc.geometry.placement import sample_packing_from_name
    from porescalemc.solvers.packing import PackingStatsSolver

    if algorithm is not None:
        mlmc_config = MLMCConfig(**{**mlmc_config.__dict__, "algorithm": algorithm})
    if n_workers != 1:
        mlmc_config = MLMCConfig(**{**mlmc_config.__dict__, "n_workers": n_workers})

    qois = qoi_names or (
        "porosity", "solid_fraction", "specific_surface_area", "mean_equivalent_radius"
    )

    def geo(level, name, cfg):
        return sample_packing_from_name(level, name, cfg)

    def sol(packing):
        return PackingStatsSolver(packing, qoi_names=qois)

    est = MLMCEstimator(
        name="run_mlmc",
        packing_config=packing_config,
        mlmc_config=mlmc_config,
        geometry_factory=geo,
        solver_factory=sol,
        qoi_fn=lambda raw, p: raw,
    )
    return est.run()
