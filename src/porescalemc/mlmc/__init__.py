"""MLMC sub-package: statistics, parallel workers, and the top-level estimator."""

from __future__ import annotations

from porescalemc.mlmc.estimator import MLMCEstimator, PairSampler
from porescalemc.mlmc.statistics import (
    MLMCStats,
    compute_mlmc_stats,
    estimate_observed_rates,
    mlmc_diagnostics_report,
    optimal_sample_counts,
    check_convergence,
)

__all__ = [
    "MLMCEstimator",
    "PairSampler",
    "MLMCStats",
    "compute_mlmc_stats",
    "estimate_observed_rates",
    "mlmc_diagnostics_report",
    "optimal_sample_counts",
    "check_convergence",
]
