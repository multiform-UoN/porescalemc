"""
High-level MLMCEstimator class.

Orchestrates the full MLMC workflow:
    1. Generate a random packing at each (level, sample) pair.
    2. Run the solver (or any QoI-returning callable).
    3. Form the pair difference [Q_l - Q_{l-1}, Q_l].
    4. Collect statistics and (optionally) adapt sample counts.
    5. Return the MLMC estimate and error bound.

Solver protocol
---------------
Any object that satisfies :class:`porescalemc.solvers.base.SolverProtocol`
works:

    solver.setup(packing)  → None
    solver.solve()         → np.ndarray  (raw solver output)
    solver.close()         → None

The ``qoi_fn`` then maps (solver_output, packing) → np.ndarray of QoI values.

Geometry factory
----------------
``geometry_factory(level, name) → Packing``

For MLMC pair estimation, the *same* geometric realisation is used at both
level l and level l−1.  The factory uses the level argument to scale the
packing according to the MLMCConfig hierarchy string (see config.hierarchy_level).

Serialisation
-------------
Results are saved/loaded with :mod:`pickle` so that long runs can be
interrupted and resumed.
"""

from __future__ import annotations

import copy
import logging
import pickle
import time
from dataclasses import dataclass
from typing import Any, Callable

import numpy as np

from porescalemc.config import MLMCConfig, PackingConfig, hierarchy_level
from porescalemc.geometry.grains import Packing
from porescalemc.geometry.placement import sample_packing_from_name
from porescalemc.mlmc.statistics import MLMCStats, check_convergence, compute_mlmc_stats, optimal_sample_counts
from porescalemc.solvers.packing import PackingStatsSolver
from porescalemc.mlmc.workers import run_all_levels

_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Sampler helpers
# ---------------------------------------------------------------------------

def identity_qoi(solver_output: Any, packing: Packing) -> np.ndarray:
    """Default QoI function: return the solver output as a numeric array."""
    del packing
    return np.asarray(solver_output, dtype=float)


@dataclass
class PairSampler:
    """Callable paired-sample generator used by the MLMC workers.

    For level l:
        - Builds a packing at level l−1 (coarse) using the geometry factory.
        - Runs the solver to get Q_{l-1}.
        - Builds the level-l packing using the same realisation name.
        - Runs the solver again to get Q_l.
        - Returns np.concatenate([Q_l - Q_{l-1}, Q_l]).

    For level 0:
        - Returns np.concatenate([Q_0, Q_0]).

    The same ``name`` is passed to the coarse and fine geometry calls.  A
    deterministic geometry factory can therefore couple the random stream for
    the two resolutions, which is essential for MLMC variance reduction.
    """
    packing_config: PackingConfig
    mlmc_config: MLMCConfig
    geometry_factory: Callable[[int, str, PackingConfig], Packing]
    solver_factory: Callable[[Packing], Any]
    qoi_fn: Callable[[Any, Packing], np.ndarray]

    def __call__(self, level: int, name: str) -> np.ndarray | None:
        """Return one paired sample or None if either solve fails."""

        # --- coarse level (l-1) ---
        q_coarse = np.zeros(0)
        if level > 0:
            p_cfg_coarse, _ = hierarchy_level(self.mlmc_config, self.packing_config, level - 1)
            try:
                # The key MLMC symmetry: pass the same `name` to both coarse and fine.
                # A deterministic geometry_factory will therefore produce the same
                # random realisation at both resolutions, maximising the covariance
                # Cov(Q_l, Q_{l-1}) and minimising Var[Q_l - Q_{l-1}].
                packing_coarse = self.geometry_factory(level - 1, name, p_cfg_coarse)
                q_coarse = self._solve_qoi(packing_coarse)
            except Exception as exc:
                _log.warning("Coarse solve failed at level %d: %s", level, exc)
                return None
            if np.isnan(q_coarse).any():
                return None

        # --- fine level (l) ---
        p_cfg_fine, _ = hierarchy_level(self.mlmc_config, self.packing_config, level)
        try:
            packing_fine = self.geometry_factory(level, name, p_cfg_fine)
            q_fine = self._solve_qoi(packing_fine)
        except Exception as exc:
            _log.warning("Fine solve failed at level %d: %s", level, exc)
            return None
        if np.isnan(q_fine).any():
            return None

        if level == 0:
            # Level 0: no coarser level, so dQ = Q_0 and the sample is [Q_0, Q_0].
            return np.concatenate([q_fine, q_fine])
        else:
            if q_fine.shape != q_coarse.shape:
                _log.warning(
                    "QoI shape mismatch at level %d: fine=%s coarse=%s",
                    level,
                    q_fine.shape,
                    q_coarse.shape,
                )
                return None
            # Return [dQ, Q] = [Q_l - Q_{l-1}, Q_l] packed into one array.
            return np.concatenate([q_fine - q_coarse, q_fine])

    def _solve_qoi(self, packing: Packing) -> np.ndarray:
        solver = self.solver_factory(packing)
        try:
            solver.setup(packing)
            raw = solver.solve()
            qoi = np.asarray(self.qoi_fn(raw, packing), dtype=float)
        finally:
            solver.close()
        return np.atleast_1d(qoi)


# ---------------------------------------------------------------------------
# Main estimator class
# ---------------------------------------------------------------------------

class MLMCEstimator:
    """Multilevel Monte Carlo estimator for expected values of QoI.

    Parameters
    ----------
    name : str
        Identifier for this run (used for file names when saving).
    packing_config : PackingConfig
        Base (level-0) geometry configuration.
    mlmc_config : MLMCConfig
        MLMC algorithm parameters.
    geometry_factory : Callable[[int, str, PackingConfig], Packing]
        Called as ``geometry_factory(level, name, packing_cfg) → Packing``.
        The factory should use the provided PackingConfig (already scaled for
        the given level via :func:`~porescalemc.config.hierarchy_level`).
    solver_factory : Callable[[Packing], SolverProtocol]
        Called as ``solver_factory(packing) → solver``.
    qoi_fn : Callable[[Any, Packing], np.ndarray]
        Maps ``(solver_output, packing) → 1-D array of QoI values``.

    Examples
    --------
    >>> from porescalemc.config import PackingConfig, MLMCConfig
    >>> from porescalemc.geometry.placement import sample_packing
    >>> from porescalemc.solvers.dummy import DummySolver
    >>>
    >>> def geo_factory(level, name, cfg): return sample_packing_from_name(level, name, cfg)
    >>> def sol_factory(packing): return DummySolver(packing)
    >>> def qoi(out, packing): return out
    >>>
    >>> est = MLMCEstimator("test", PackingConfig(), MLMCConfig(n_levels=2, m0=8),
    ...                     geo_factory, sol_factory, qoi)
    >>> mean, error = est.run()
    """

    def __init__(
        self,
        name: str,
        packing_config: PackingConfig,
        mlmc_config: MLMCConfig,
        geometry_factory: Callable | None = None,
        solver_factory: Callable | None = None,
        qoi_fn: Callable | None = None,
    ):
        self.name = name
        # Deep-copy both configs so that _add_level can safely mutate self.mlmc_config
        # without affecting the caller's original objects.
        self.packing_config = copy.deepcopy(packing_config)
        self.mlmc_config = copy.deepcopy(mlmc_config)
        self.geometry_factory = geometry_factory or sample_packing_from_name
        self.solver_factory = solver_factory or PackingStatsSolver
        self.qoi_fn = qoi_fn or identity_qoi

        self._samples: list[list[np.ndarray]] = [[] for _ in range(mlmc_config.n_levels)]
        self._work: list[float] = [0.0] * mlmc_config.n_levels
        self._stats: MLMCStats | None = None
        self._sampler = PairSampler(
            packing_config,
            mlmc_config,
            self.geometry_factory,
            self.solver_factory,
            self.qoi_fn,
        )

    # ------------------------------------------------------------------ #
    # Public interface                                                     #
    # ------------------------------------------------------------------ #

    def run(self) -> tuple[np.ndarray, np.ndarray]:
        """Run the MLMC estimator.

        Dispatches to the adaptive (Giles) or fixed algorithm depending on
        ``mlmc_config.algorithm``.

        Returns
        -------
        mean : np.ndarray
            MLMC estimate of E[QoI].
        error : np.ndarray
            Estimated total error (|bias| + stat_error).
        """
        cfg = self.mlmc_config
        t0 = time.perf_counter()

        if cfg.algorithm == "adaptive":
            self._run_adaptive()
        else:
            self._run_fixed()

        elapsed = time.perf_counter() - t0
        _log.info("MLMC run complete in %.2fs", elapsed)

        if self._stats is None:
            raise RuntimeError("MLMC run did not produce statistics.")
        stats = self._stats
        mean = stats.estimator
        error = np.abs(stats.bias) + stats.stat_error
        return mean, error

    @property
    def stats(self) -> MLMCStats | None:
        """Most recently computed statistics, or None before ``run``."""
        return self._stats

    def save(self, path: str) -> None:
        """Save samples, work, stats, and configs to ``path`` (pickle)."""
        data = {
            "name": self.name,
            "samples": self._samples,
            "work": self._work,
            "stats": self._stats,
            "packing_config": self.packing_config,
            "mlmc_config": self.mlmc_config,
        }
        with open(path, "wb") as f:
            pickle.dump(data, f)
        _log.info("Saved estimator state to %s", path)

    @classmethod
    def load_data(cls, path: str) -> dict:
        """Load saved estimator data from ``path``."""
        with open(path, "rb") as f:
            return pickle.load(f)

    @classmethod
    def load(cls, path: str) -> "MLMCEstimator":
        """Reconstruct an estimator from a saved file.

        Configs are restored from the pickle.  Supply geometry/solver/qoi
        factories before calling ``run()`` again to continue sampling.
        """
        data = cls.load_data(path)
        est = cls(
            name=data.get("name", "loaded"),
            packing_config=data.get("packing_config", PackingConfig()),
            mlmc_config=data.get("mlmc_config", MLMCConfig()),
        )
        est._samples = data.get("samples", [])
        est._work = data.get("work", [0.0] * len(est._samples))
        est._stats = data.get("stats")
        return est

    def plot_convergence(self, filename: str | None = None) -> None:
        """Plot MLMC convergence diagnostics.

        Shows per-level mean and variance of the pair differences.
        Requires matplotlib (optional dependency).

        Parameters
        ----------
        filename : str, optional
            If given, save the figure to this path; otherwise display it.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            _log.error("matplotlib is required for plotting. Install with: pip install matplotlib")
            return

        if self._stats is None:
            _log.warning("No statistics available. Run the estimator first.")
            return

        stats = self._stats
        levels = list(range(len(stats.means)))
        nvar = len(stats.means[0])

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        for v in range(nvar):
            means_v = [abs(m[v]) for m in stats.means]
            vars_v = [v_arr[v] for v_arr in stats.variances]

            axes[0].semilogy(levels, means_v, "o-", label=f"QoI {v}")
            axes[1].semilogy(levels, vars_v, "s--", label=f"QoI {v}")

        axes[0].set_xlabel("Level ℓ")
        axes[0].set_ylabel("|E[Q_ℓ − Q_{ℓ−1}]|")
        axes[0].set_title("Mean of level differences")
        axes[0].grid(True)
        axes[0].legend()

        axes[1].set_xlabel("Level ℓ")
        axes[1].set_ylabel("Var[Q_ℓ − Q_{ℓ−1}]")
        axes[1].set_title("Variance of level differences")
        axes[1].grid(True)
        axes[1].legend()

        fig.tight_layout()
        if filename:
            fig.savefig(filename, bbox_inches="tight")
            _log.info("Saved convergence plot to %s", filename)
        else:
            plt.show()
        plt.close(fig)

    # ------------------------------------------------------------------ #
    # Internal algorithms                                                  #
    # ------------------------------------------------------------------ #

    def _run_fixed(self) -> None:
        """Non-adaptive MLMC: run fixed number of samples per level."""
        cfg = self.mlmc_config
        n_per_level = cfg.initial_sample_counts
        _log.info("Fixed MLMC: levels=%d, samples=%s", cfg.n_levels, n_per_level)

        new_samples, new_work = run_all_levels(
            self._sampler, n_per_level, self.name, cfg.n_workers
        )
        self._merge(new_samples, new_work)
        self._compute_stats()

    def _run_adaptive(self) -> None:
        """Adaptive MLMC: iteratively update sample counts until convergence."""
        cfg = self.mlmc_config
        n_per_level = cfg.initial_sample_counts
        _log.info("Adaptive MLMC: initial levels=%d, samples=%s", cfg.n_levels, n_per_level)

        # max_outer caps iterations in case the convergence check is never satisfied
        # (e.g. if the bias never decays because the QoI has no level dependence).
        max_outer = 20
        for iteration in range(max_outer):
            new_samples, new_work = run_all_levels(
                self._sampler, n_per_level, self.name, cfg.n_workers
            )
            self._merge(new_samples, new_work)
            self._compute_stats()

            converged, reason = check_convergence(self._stats, cfg)
            _log.info("Iteration %d: %s", iteration, reason)
            if converged:
                break

            added_level = False
            if self._bias_too_large() and len(self._samples) < cfg.max_levels:
                self._add_level()
                added_level = True

            optimal = optimal_sample_counts(self._stats, cfg)
            if added_level:
                # New finest level gets the minimum initial allocation.
                optimal = list(optimal) + [cfg.min_samples]
            # Draw only the *additional* samples needed (don't redo existing ones).
            n_per_level = [
                max(0, optimal[l] - len(self._samples[l]))
                for l in range(len(optimal))
            ]
            n_per_level = [max(cfg.min_samples, n) if len(self._samples[l]) == 0 else n
                           for l, n in enumerate(n_per_level)]
        else:
            _log.warning("Adaptive MLMC: reached max iterations (%d).", max_outer)

    def _merge(self, new_samples: list[list[np.ndarray]], new_work: list[float]) -> None:
        """Merge newly collected samples into ``self._samples``."""
        for l, (samps, work) in enumerate(zip(new_samples, new_work)):
            if l >= len(self._samples):
                self._samples.append([])
                self._work.append(0.0)
            if self.mlmc_config.reuse_samples:
                self._samples[l].extend(samps)
            else:
                self._samples[l] = samps
            self._work[l] += work

    def _compute_stats(self) -> None:
        """Recompute statistics from all accumulated samples."""
        self._stats = compute_mlmc_stats(
            self._samples, self._work, self.mlmc_config
        )

    def _bias_too_large(self) -> bool:
        """Return True if the bias exceeds its tolerance fraction."""
        if self._stats is None:
            return True
        tol = self.mlmc_config.tolerance
        theta = self.mlmc_config.error_split
        Q_mag = np.maximum(np.abs(self._stats.estimator), 1e-15)
        rel_bias = np.abs(self._stats.bias) / Q_mag
        return bool((rel_bias >= theta * tol).any())

    def _add_level(self) -> None:
        """Add a new finest level to the hierarchy.

        We mutate self.mlmc_config here — safe because __init__ made a deep copy,
        so the caller's original MLMCConfig object is untouched.
        """
        self._samples.append([])
        self._work.append(0.0)
        self.mlmc_config.n_levels = max(self.mlmc_config.n_levels, len(self._samples))
        _log.info("Added level %d (total levels: %d)", len(self._samples) - 1, len(self._samples))
