"""
MLMC statistics: error estimation and sample-count optimisation.

All functions are pure (no I/O, no global state) and operate on plain numpy
arrays.  They port the core math from the original mlmc.py, fixing the
Python 3 integer-division bug (``len(avg)/2`` → ``len(avg)//2``) and
removing all global variable dependencies.

Sample format
-------------
Each MLMC sample at level l is the concatenated vector::

    s = [dQ, Q]  with dQ = Q_l - Q_{l-1},  Q = Q_l

For l = 0 (no coarser level): ``dQ = Q_0``, so ``s = [Q_0, Q_0]``.

This doubled format lets us compute both the level-l increment statistics
(from ``dQ``) and the coarser-level value (from ``Q - dQ = Q_{l-1}``) from a
single sample vector, without storing two separate arrays.

MLMC estimator
--------------
The MLMC estimate of E[Q] is::

    Q_MLMC = Σ_l E[Q_l - Q_{l-1}]  ≈  Σ_l mean(dQ_l)

Statistical error (CLT bound)::

    e_stat = z_{α/2} · sqrt(Σ_l Var[dQ_l] / M_l)

Bias (last-level proxy)::

    e_bias ≈ |E[dQ_L]|  (assumes geometric decay of bias)

Optimal sample counts (Giles 2015)::

    M_l* ∝ sqrt(Var_l / W_l) · (z_{α/2} / ((1-θ) · TOL · E[Q]))^2

where θ = error_split and W_l is the average work per sample at level l.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
from scipy.stats import norm as scipy_norm

from porescalemc.config import MLMCConfig

_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data container
# ---------------------------------------------------------------------------

@dataclass
class MLMCStats:
    """Per-level statistics computed from MLMC paired samples.

    Attributes
    ----------
    means : list of np.ndarray
        E[Q_l - Q_{l-1}] per level, shape (nvar,) each.
    variances : list of np.ndarray
        Var[Q_l - Q_{l-1}] per level, shape (nvar,) each.
    n_samples : list of int
        Effective sample count per level.
    work : list of float
        Average wall-clock seconds per sample per level.
    estimator : np.ndarray, shape (nvar,)
        MLMC estimate = sum of means.
    stat_error : np.ndarray, shape (nvar,)
        Statistical error bound at the requested confidence level.
    bias : np.ndarray, shape (nvar,)
        Last-level mean (proxy for discretisation bias).
    """

    means: list[np.ndarray]
    variances: list[np.ndarray]
    n_samples: list[int]
    work: list[float]
    estimator: np.ndarray
    stat_error: np.ndarray
    bias: np.ndarray

    def summary(self) -> str:
        """Human-readable summary string."""
        lines = [
            f"MLMC estimate:       {self.estimator}",
            f"Statistical error:   {self.stat_error}",
            f"Bias proxy:          {self.bias}",
            f"Samples per level:   {self.n_samples}",
        ]
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Sample variance (unbiased, Bessel-corrected)
# ---------------------------------------------------------------------------

def _sample_var(x: list[np.ndarray], axis: int = 0) -> np.ndarray:
    """Unbiased sample variance along ``axis`` (denominator N-1).

    Equivalent to ``np.var(x, ddof=1, axis=axis)`` but works on a list of
    arrays (consistent with the original code's ``samplevar``).
    """
    arr = np.array(x)
    n = arr.shape[axis]
    if n < 2:
        return np.zeros(arr.shape[1:] if axis == 0 else arr.shape[:axis])
    return np.var(arr, axis=axis, ddof=1)


# ---------------------------------------------------------------------------
# Core statistics
# ---------------------------------------------------------------------------

def compute_mlmc_stats(
    samples: list[list[np.ndarray]],
    work: list[float],
    config: MLMCConfig,
) -> MLMCStats:
    """Compute MLMC statistics from nested per-level paired samples.

    Parameters
    ----------
    samples : list of lists of np.ndarray
        ``samples[l][i]`` is the i-th sample at level l, a 1-D array of
        length 2*nvar: first nvar entries are dQ = Q_l - Q_{l-1}, last nvar
        are Q = Q_l.  For l=0, both halves are Q_0.
    work : list of float
        Total (summed) wall-clock time per level.  Average per-sample work is
        work[l] / len(samples[l]).
    config : MLMCConfig
        Uses ``confidence`` for the z-score.

    Returns
    -------
    MLMCStats
    """
    if not samples or not any(samples):
        raise ValueError("samples must contain at least one non-empty level.")

    calpha = float(scipy_norm.ppf(config.confidence))

    for l, sl in enumerate(samples):
        if len(sl) == 0:
            raise ValueError(f"Level {l} has no samples.")

    # The doubled-vector format [dQ, Q] means each sample has 2*nvar components.
    # Integer division is essential here (Python 3 fix: / gives float, // gives int).
    nvar = len(samples[0][0]) // 2

    means = []
    variances = []
    n_samples_list = []
    work_per_sample = []

    for l, sl in enumerate(samples):
        arr = np.array(sl)          # shape (M_l, 2*nvar)
        dq = arr[:, :nvar]          # extract dQ = Q_l - Q_{l-1} (first half)
        means.append(np.mean(dq, axis=0))
        variances.append(_sample_var(list(dq), axis=0))
        n_samples_list.append(len(sl))
        w_total = work[l] if l < len(work) else 0.0
        work_per_sample.append(w_total / max(1, len(sl)))

    # MLMC estimator = sum of all level-increment means.
    estimator = np.sum(means, axis=0)

    # CLT statistical error bound: e_stat = z_α · sqrt(Σ_l Var_l / M_l).
    # Dividing each level's variance by its sample count gives the variance of
    # the level-l mean; summing gives the variance of the MLMC sum.
    varnorm = np.array([v / max(1, m) for v, m in zip(variances, n_samples_list)])
    stat_error = calpha * np.sqrt(np.sum(varnorm, axis=0))

    # Bias proxy calculation depending on mode.
    bias_mode = getattr(config, "bias_mode", "last_increment")
    if bias_mode == "adaptive" and len(means) >= 2:
        # Adaptive bias estimate: assuming |E[dQ_l]| ~ C · r^{-αl} (geometric decay),
        # the remaining bias |E[Q] - Q_MLMC| ≈ |E[dQ_L]|² / (|E[dQ_{L-1}]| - |E[dQ_L]|).
        abs_l = np.abs(means[-1])
        abs_lm1 = np.abs(means[-2])
        denom = abs_lm1 - abs_l
        bias_adaptive = np.zeros_like(means[-1])
        for idx in range(nvar):
            if denom[idx] > 1e-9:
                val = (abs_l[idx] ** 2) / denom[idx]
                bias_adaptive[idx] = np.sign(means[-1][idx]) * val
            else:
                bias_adaptive[idx] = means[-1][idx]  # fallback: use last increment directly
        bias = bias_adaptive
    else:
        # Default: last-level increment as proxy for remaining bias.
        # If bias decays geometrically, |E[dQ_L]| ≈ |E[Q] - Q_MLMC| / (r^α - 1).
        bias = means[-1].copy()

    # Optional Richardson extrapolation to remove leading-order bias.
    # The MLMC estimator Q_MLMC estimates E[Q_L], not E[Q].  The remaining bias is
    #   E[Q] - E[Q_L] ≈ E[dQ_{L+1}] ≈ E[dQ_L] / r^α
    # so the Richardson correction is E[dQ_L] / (r^α - 1), which is means[-1] / factor.
    if getattr(config, "extrapolate_bias", False) and len(means) >= 2:
        alpha = getattr(config, "alpha", 2.0)
        factor = (2.0 ** alpha - 1.0)   # r = 2 is the standard refinement ratio
        if factor > 1e-9:
            correction = means[-1] / factor
            estimator = estimator + correction
            bias = np.zeros_like(bias)  # bias largely removed by extrapolation

    _log.debug("MLMC stats: estimate=%s, stat_error=%s, bias=%s", estimator, stat_error, bias)

    return MLMCStats(
        means=means,
        variances=variances,
        n_samples=n_samples_list,
        work=work_per_sample,
        estimator=estimator,
        stat_error=stat_error,
        bias=bias,
    )


# ---------------------------------------------------------------------------
# Optimal sample allocation
# ---------------------------------------------------------------------------

def optimal_sample_counts(
    stats: MLMCStats,
    config: MLMCConfig,
) -> list[int]:
    """Compute optimal per-level sample counts using the Giles formula.

    The optimal allocation minimises total work subject to the statistical
    error constraint::

        M_l* = ceil( z_α² / ((1-θ)² TOL² E[Q]²)
                     · sqrt(Var_l / W_l)
                     · Σ_l' sqrt(Var_{l'} W_{l'}) )

    Parameters
    ----------
    stats : MLMCStats
        Current statistics (variances and work per level).
    config : MLMCConfig
        Uses ``confidence``, ``tolerance``, ``error_split``, ``min_samples``.

    Returns
    -------
    list of int
        Optimal M_l for each level, at least ``config.min_samples``.
    """
    calpha = float(scipy_norm.ppf(config.confidence))
    tol = config.tolerance
    theta = config.error_split   # fraction of TOL reserved for bias; (1-θ) for stat error
    n_levels = len(stats.means)
    nvar = len(stats.estimator)

    variances = np.array([v for v in stats.variances])   # (L, nvar)
    work = np.array(stats.work)                           # (L,)
    work = np.maximum(work, 1e-10)                        # avoid division by zero

    # For multi-component QoI, drive sample counts by the most demanding component.
    var_max = variances.max(axis=1)   # (L,)

    # Scale by estimator magnitude for relative (not absolute) tolerance.
    Q_mag = np.abs(stats.estimator).max()
    Q_mag = max(float(Q_mag), 1e-10)

    # Giles optimal allocation: M_l* ∝ sqrt(Var_l/W_l) · (normalisation constant).
    # This minimises Σ M_l W_l subject to z_α² Σ Var_l/M_l ≤ ((1-θ)·TOL·|Q|)².
    coeff = (calpha / ((1.0 - theta) * tol * Q_mag)) ** 2
    sqrt_var_over_w = np.sqrt(var_max / work)
    sum_sqrt_var_w = float(np.sum(np.sqrt(var_max * work)))
    optimal = coeff * sqrt_var_over_w * sum_sqrt_var_w

    counts = [max(config.min_samples, int(np.ceil(m))) for m in optimal]
    _log.debug("Optimal sample counts: %s", counts)
    return counts


# ---------------------------------------------------------------------------
# Convergence check
# ---------------------------------------------------------------------------

def check_convergence(
    stats: MLMCStats,
    config: MLMCConfig,
) -> tuple[bool, str]:
    """Check whether the MLMC estimator has converged.

    Both the statistical error and the bias must be within their respective
    fractions of the total tolerance::

        |bias| / |estimate|   < θ · TOL
        stat_error / |estimate| < (1 - θ) · TOL

    The split θ = error_split separates bias (deterministic, removed by refining
    the hierarchy) from stat error (stochastic, removed by adding samples).
    Both criteria must be satisfied simultaneously for convergence.

    Parameters
    ----------
    stats : MLMCStats
    config : MLMCConfig

    Returns
    -------
    converged : bool
    reason : str
        Human-readable explanation.
    """
    tol = config.tolerance
    theta = config.error_split  # θ fraction of TOL assigned to bias
    Q_mag = np.abs(stats.estimator)
    Q_mag = np.maximum(Q_mag, 1e-15)  # avoid division by zero for near-zero estimators

    rel_bias = np.abs(stats.bias) / Q_mag
    rel_stat = np.abs(stats.stat_error) / Q_mag

    bias_ok = bool((rel_bias < theta * tol).all())
    stat_ok = bool((rel_stat < (1.0 - theta) * tol).all())

    if bias_ok and stat_ok:
        return True, "Converged: both bias and statistical error within tolerance."
    elif not bias_ok:
        reason = (
            f"Not converged: relative bias {rel_bias.max():.3e} "
            f"exceeds {theta * tol:.3e}."
        )
    else:
        reason = (
            f"Not converged: relative stat error {rel_stat.max():.3e} "
            f"exceeds {(1-theta) * tol:.3e}."
        )
    return False, reason
