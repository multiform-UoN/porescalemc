"""
Parallel sample execution using :mod:`concurrent.futures`.

Workers run in separate processes to bypass the GIL and achieve true
parallelism.  Each call to ``sampler(level, name)`` is independent.

NaN results (indicating failed/rejected samples) are silently dropped.
"""

from __future__ import annotations

import logging
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Callable

import numpy as np

_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Worker function (must be top-level for pickling)
# ---------------------------------------------------------------------------

def _call_sampler(sampler_fn, level: int, name: str):
    """Execute sampler_fn(level, name) and return (level, result, elapsed)."""
    t0 = time.perf_counter()
    try:
        result = sampler_fn(level, name)
    except Exception as exc:
        _log.warning("Sampler raised at level %d name %s: %s", level, name, exc)
        result = None
    elapsed = time.perf_counter() - t0
    return level, result, elapsed


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_level_samples(
    sampler: Callable[[int, str], np.ndarray | None],
    level: int,
    n_samples: int,
    base_name: str,
    n_workers: int = 1,
) -> tuple[list[np.ndarray], float]:
    """Run ``n_samples`` calls of ``sampler(level, name)`` in parallel.

    Parameters
    ----------
    sampler : Callable[[int, str], np.ndarray | None]
        Function that draws one MLMC paired sample at the given level.
        Must be picklable (i.e., a module-level function or an object with
        a ``__call__`` method that pickles cleanly).
        Should return ``None`` or a NaN array to signal a failed sample.
    level : int
        MLMC level index.
    n_samples : int
        Number of samples to draw.
    base_name : str
        Prefix for sample names (useful for logging and file naming).
    n_workers : int
        Number of parallel processes.  1 → serial execution in the current
        process (avoids pickle overhead for trivial samplers).

    Returns
    -------
    results : list of np.ndarray
        Valid (non-NaN) samples.  May be shorter than ``n_samples``.
    total_time : float
        Sum of elapsed times across all samples (wall-clock seconds).
    """
    results = []
    total_time = 0.0

    names = [f"{base_name}_l{level:02d}_{i:06d}" for i in range(n_samples)]

    if n_workers <= 1:
        # Serial path: avoids ProcessPoolExecutor startup overhead, which is significant
        # for small samplers (packing generation is fast; forking is not free).
        for name in names:
            _, result, elapsed = _call_sampler(sampler, level, name)
            total_time += elapsed
            if _is_valid(result):
                results.append(np.asarray(result))
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = {
                pool.submit(_call_sampler, sampler, level, name): name
                for name in names
            }
            for future in as_completed(futures):
                _, result, elapsed = future.result()
                total_time += elapsed
                if _is_valid(result):
                    results.append(np.asarray(result))

    _log.debug(
        "Level %d: %d/%d valid samples, total_time=%.2fs",
        level, len(results), n_samples, total_time,
    )
    return results, total_time


def run_all_levels(
    sampler: Callable[[int, str], np.ndarray | None],
    n_samples_per_level: list[int],
    base_name: str,
    n_workers: int = 1,
) -> tuple[list[list[np.ndarray]], list[float]]:
    """Run samples at every level and collect results.

    Parameters
    ----------
    sampler : Callable
        As in :func:`run_level_samples`.
    n_samples_per_level : list of int
        Number of samples to draw at each level.
    base_name : str
    n_workers : int

    Returns
    -------
    all_samples : list[list[np.ndarray]]
        all_samples[l] = list of valid sample arrays at level l.
    work_per_level : list[float]
        Total wall-clock time per level.
    """
    all_samples = []
    work_per_level = []

    for level, n in enumerate(n_samples_per_level):
        samples, work = run_level_samples(sampler, level, n, base_name, n_workers)
        all_samples.append(samples)
        work_per_level.append(work)

    return all_samples, work_per_level


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _is_valid(result) -> bool:
    """Return True if result is a non-None, non-NaN array.

    NaN is the sentinel value used by PairSampler to signal a failed solve
    (solver exception, shape mismatch, etc.).  Such samples are discarded
    rather than propagated as errors.
    """
    if result is None:
        return False
    arr = np.asarray(result)
    if arr.ndim == 0:
        return not bool(np.isnan(arr))
    return bool(~np.isnan(arr).any())
