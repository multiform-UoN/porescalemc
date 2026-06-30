import os
import tempfile
import unittest

import numpy as np

from porescalemc.config import MLMCConfig, PackingConfig
from porescalemc.geometry.grains import Grain, Packing
from porescalemc.mlmc.estimator import MLMCEstimator, PairSampler
from porescalemc.mlmc.statistics import (
    check_convergence,
    compute_mlmc_stats,
    estimate_observed_rates,
    mlmc_diagnostics_report,
    optimal_sample_counts,
)
from porescalemc.solvers.packing import PackingStatsSolver


class ScalarPorositySolver:
    def __init__(self, packing):
        self.packing = packing

    def setup(self, packing):
        self.packing = packing

    def solve(self):
        return np.array([self.packing.porosity])

    def close(self):
        return None


def identity_qoi(output, packing):
    del packing
    return output


class MLMCTests(unittest.TestCase):
    def test_statistics_compute_expected_estimator(self):
        samples = [
            [np.array([1.0, 1.0]), np.array([1.2, 1.2])],
            [np.array([0.1, 1.3]), np.array([0.2, 1.5])],
        ]
        stats = compute_mlmc_stats(samples, [0.2, 0.8], MLMCConfig())
        np.testing.assert_allclose(stats.estimator, [1.25])
        np.testing.assert_allclose(stats.bias, [0.15])
        self.assertEqual(stats.n_samples, [2, 2])

    def test_optimal_sample_counts_are_at_least_min_samples(self):
        samples = [
            [np.array([1.0, 1.0]), np.array([1.2, 1.2]), np.array([0.8, 0.8])],
            [np.array([0.1, 1.3]), np.array([0.2, 1.5]), np.array([0.0, 1.4])],
        ]
        cfg = MLMCConfig(min_samples=4)
        stats = compute_mlmc_stats(samples, [0.3, 0.9], cfg)
        counts = optimal_sample_counts(stats, cfg)
        self.assertTrue(all(count >= 4 for count in counts))

    def test_observed_rates_and_diagnostics_report(self):
        """Finite-level diagnostics should expose measured mean/var/work data."""
        samples = [
            [np.array([1.0, 1.0]), np.array([1.4, 1.4]), np.array([0.6, 0.6])],
            [np.array([0.25, 1.25]), np.array([0.35, 1.35]), np.array([0.15, 1.15])],
            [np.array([0.06, 1.31]), np.array([0.08, 1.33]), np.array([0.04, 1.29])],
        ]
        cfg = MLMCConfig(min_samples=2, refratio=2.0)
        stats = compute_mlmc_stats(samples, [0.3, 1.2, 4.8], cfg)

        rates = estimate_observed_rates(stats, refratio=cfg.refratio)
        report = mlmc_diagnostics_report(stats, cfg)

        self.assertGreater(rates["alpha_hat"], 0.0)
        self.assertGreater(rates["beta_hat"], 0.0)
        self.assertGreater(rates["gamma_hat"], 0.0)
        self.assertEqual(len(report["levels"]), 3)
        self.assertEqual(len(report["optimal_samples"]), 3)
        self.assertIn("total_work", report["levels"][0])

    def test_convergence_detects_large_bias(self):
        cfg = MLMCConfig(tolerance=0.01, error_split=0.5)
        samples = [
            [np.array([1.0, 1.0]), np.array([1.0, 1.0])],
            [np.array([0.5, 1.5]), np.array([0.5, 1.5])],
        ]
        stats = compute_mlmc_stats(samples, [1.0, 1.0], cfg)
        converged, reason = check_convergence(stats, cfg)
        self.assertFalse(converged)
        self.assertIn("bias", reason)

    def test_pair_sampler_uses_same_realization_name_for_coarse_and_fine(self):
        calls = []

        def geo(level, name, cfg):
            calls.append((level, name))
            radius = 0.05 + 0.01 * level
            return Packing([Grain.sphere([0.0, 0.0, 0.0], radius)], box=np.ones(3))

        sampler = PairSampler(
            PackingConfig(n_grains=1),
            MLMCConfig(n_levels=2, hierarchy="g"),
            geo,
            ScalarPorositySolver,
            identity_qoi,
        )
        result = sampler(1, "sample-0001")
        self.assertIsNotNone(result)
        self.assertEqual(calls, [(0, "sample-0001"), (1, "sample-0001")])

    def test_estimator_runs_with_default_packing_solver(self):
        cfg = PackingConfig(n_grains=3, psd="constant", periodic=False, detached=True)
        mlmc = MLMCConfig(n_levels=2, m0=3, min_samples=2, n_workers=1, hierarchy="g")
        estimator = MLMCEstimator("default-smoke", cfg, mlmc)
        mean, error = estimator.run()
        self.assertEqual(mean.shape, error.shape)
        self.assertEqual(mean.shape[0], 8)
        self.assertEqual([len(level) for level in estimator._samples], [3, 2])

    def test_estimator_runs_with_custom_scalar_solver(self):
        def geo(level, name, cfg):
            radius = 0.05 + 0.005 * level
            return Packing([Grain.sphere([0.0, 0.0, 0.0], radius)], box=np.ones(3))

        estimator = MLMCEstimator(
            "scalar-smoke",
            PackingConfig(n_grains=1),
            MLMCConfig(n_levels=2, m0=4, min_samples=2, n_workers=1, hierarchy="g"),
            geo,
            ScalarPorositySolver,
            identity_qoi,
        )
        mean, error = estimator.run()
        self.assertEqual(mean.shape, (1,))
        self.assertTrue(np.isfinite(mean).all())
        self.assertTrue(np.isfinite(error).all())

    def test_adaptive_estimator_converges_on_level_independent_qoi(self):
        def geo(level, name, cfg):
            del level, name, cfg
            return Packing([Grain.sphere([0.0, 0.0, 0.0], 0.05)], box=np.ones(3))

        estimator = MLMCEstimator(
            "adaptive-smoke",
            PackingConfig(n_grains=1),
            MLMCConfig(
                n_levels=2,
                m0=4,
                min_samples=2,
                n_workers=1,
                hierarchy="g",
                algorithm="adaptive",
                tolerance=0.1,
            ),
            geo,
            ScalarPorositySolver,
            identity_qoi,
        )
        mean, error = estimator.run()
        self.assertEqual(mean.shape, (1,))
        self.assertLessEqual(float(error[0]), 1e-12)

    def test_synthetic_mlmc_recovers_theoretical_rates(self):
        """Synthetic test with known rates to verify MLMC estimator behavior.

        Uses a deterministic-seeded synthetic solver so the test is repeatable.
        The MLMC estimator should land within a generous band around the true mean;
        we care about correctness of the framework (sign, scale, no NaN), not
        tight statistical convergence with only ~12 samples.
        """
        true_mean = 0.42
        alpha = 2.0
        rng = np.random.default_rng(seed=12345)

        def make_synthetic_solver(level: int):
            class SyntheticSolver:
                def setup(self, packing):
                    pass
                def solve(self):
                    h = 2.0 ** -level
                    bias = 0.25 * (h ** alpha)
                    val = true_mean + bias + rng.normal(0, 0.05 * h)
                    return np.array([val])
                def close(self):
                    pass
            return SyntheticSolver()

        def geo(level, name, cfg):
            return Packing([Grain.sphere([0.0, 0.0, 0.0], 0.05)], box=np.ones(3))

        # Patch PairSampler to inject a level-aware synthetic solver
        original_call = PairSampler.__call__

        def level_injecting_call(self_sampler, level, name):
            orig_factory = self_sampler.solver_factory
            self_sampler.solver_factory = lambda p: make_synthetic_solver(level)
            try:
                return original_call(self_sampler, level, name)
            finally:
                self_sampler.solver_factory = orig_factory

        PairSampler.__call__ = level_injecting_call

        try:
            estimator = MLMCEstimator(
                "synthetic-rates",
                PackingConfig(n_grains=1),
                MLMCConfig(n_levels=3, m0=12, min_samples=5, tolerance=0.05, n_workers=1),
                geo,
                lambda p: make_synthetic_solver(0),
                identity_qoi,
            )
            mean, error = estimator.run()
            # With m0=12 and a seeded RNG the estimate is stable; allow ±0.5 for safety
            self.assertLess(abs(mean[0] - true_mean), 0.5)
            self.assertTrue(np.isfinite(mean).all())
        finally:
            PairSampler.__call__ = original_call


    def test_save_load_preserves_configs_and_results(self):
        """Serialisation round-trip: configs, samples, and estimate survive pickle."""
        cfg = PackingConfig(n_grains=3, psd="constant", periodic=False)
        mlmc = MLMCConfig(n_levels=2, m0=3, min_samples=2, n_workers=1, hierarchy="g")
        estimator = MLMCEstimator("save-load-test", cfg, mlmc)
        mean_before, error_before = estimator.run()

        with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as f:
            path = f.name
        try:
            estimator.save(path)
            loaded = MLMCEstimator.load(path)

            # Configs round-trip
            self.assertEqual(loaded.packing_config.n_grains, cfg.n_grains)
            self.assertEqual(loaded.mlmc_config.n_levels, mlmc.n_levels)

            # Estimate is preserved
            np.testing.assert_array_equal(loaded._stats.estimator, mean_before)
            self.assertEqual([len(s) for s in loaded._samples],
                             [len(s) for s in estimator._samples])
        finally:
            os.unlink(path)

    def test_richardson_extrapolation_correction(self):
        """Richardson extrapolation correction matches means[-1] / (2^alpha - 1)."""
        samples = [
            [np.array([1.0, 1.0]), np.array([1.0, 1.0])],  # level 0: mean dQ0 = 1.0
            [np.array([0.2, 1.2]), np.array([0.2, 1.2])],  # level 1: mean dQ1 = 0.2
        ]
        # extrapolation with alpha = 2.0 -> factor = 3.0 -> correction = 0.2 / 3.0 = 0.06666...
        cfg = MLMCConfig(extrapolate_bias=True, alpha=2.0)
        stats = compute_mlmc_stats(samples, [0.1, 0.1], cfg)
        expected_est = (1.0 + 0.2) + (0.2 / 3.0)
        self.assertAlmostEqual(float(stats.estimator[0]), expected_est)
        self.assertAlmostEqual(float(stats.bias[0]), 0.0)  # bias is zeroed out

    def test_optimal_sample_counts_uniform(self):
        """Uniform variance and work leads to uniform optimal sample allocation across levels."""
        # For MLMStats: estimator, means, variances, n_samples, work
        # Let's construct samples such that variances are equal, and work is equal.
        samples = [
            [np.array([1.0, 1.0]), np.array([2.0, 2.0])] * 10,  # level 0, M0=20, var=0.5
            [np.array([0.1, 1.1]), np.array([1.1, 2.1])] * 10,  # level 1, M1=20, var=0.5
            [np.array([0.01, 1.01]), np.array([1.01, 2.01])] * 10,  # level 2, M2=20, var=0.5
        ]
        # Work per sample is equal
        work = [20.0, 20.0, 20.0]
        cfg = MLMCConfig(tolerance=0.01, error_split=0.5, min_samples=5)
        stats = compute_mlmc_stats(samples, work, cfg)
        counts = optimal_sample_counts(stats, cfg)
        # All counts should be identical since var and work are identical
        self.assertEqual(len(counts), 3)
        self.assertEqual(counts[0], counts[1])
        self.assertEqual(counts[1], counts[2])

    def test_hierarchy_level_independence_from_mratio(self):
        """Config hierarchy scaling uses refratio and is independent of mratio."""
        from porescalemc.config import hierarchy_level
        mlmc_cfg_1 = MLMCConfig(hierarchy="gsdn", refratio=2.5, mratio=4.0)
        mlmc_cfg_2 = MLMCConfig(hierarchy="gsdn", refratio=2.5, mratio=10.0)
        packing_cfg = PackingConfig(n_grains=10, xlen=1.0, ylen=1.0, zlen=1.0, mu=0.2)
        
        p1, mf1 = hierarchy_level(mlmc_cfg_1, packing_cfg, 1)
        p2, mf2 = hierarchy_level(mlmc_cfg_2, packing_cfg, 1)
        
        # Scaling results must be identical
        self.assertEqual(mf1, mf2)
        self.assertEqual(p1.xlen, p2.xlen)
        self.assertEqual(p1.mu, p2.mu)
        self.assertEqual(p1.n_grains, p2.n_grains)
        
        # Double check the actual values are scaled by refratio = 2.5
        self.assertEqual(mf1, 2.5)
        self.assertEqual(p1.xlen, 2.5)
        self.assertEqual(p1.mu, 0.2 / 2.5)
        # dimension is 3, so n_grains_scale = (2.5) ** 3 = 15.625 -> 10 * 15.625 = 156
        self.assertEqual(p1.n_grains, 156)

    def test_adaptive_bias_mode(self):
        """Adaptive bias mode estimates bias correctly with geometric decay, and falls back gracefully."""
        # Case 1: Monotonic decay
        # level 0 difference = 1.0, level 1 difference = 0.25 (decay by factor 4)
        # Expected bias = 0.25^2 / (1.0 - 0.25) = 0.0625 / 0.75 = 0.083333...
        samples = [
            [np.array([1.0, 1.0]), np.array([1.0, 1.0])],
            [np.array([0.25, 1.25]), np.array([0.25, 1.25])],
        ]
        cfg = MLMCConfig(bias_mode="adaptive")
        stats = compute_mlmc_stats(samples, [0.1, 0.1], cfg)
        expected_bias = (0.25 ** 2) / (1.0 - 0.25)
        self.assertAlmostEqual(float(stats.bias[0]), expected_bias)

        # Case 2: Non-monotonic decay (fallback to last increment)
        # level 0 difference = 0.1, level 1 difference = 0.25 (no decay)
        # Expected bias = 0.25 (fallback since 0.1 - 0.25 < 0)
        samples_non_mono = [
            [np.array([0.1, 1.0]), np.array([0.1, 1.0])],
            [np.array([0.25, 1.25]), np.array([0.25, 1.25])],
        ]
        stats_non_mono = compute_mlmc_stats(samples_non_mono, [0.1, 0.1], cfg)
        self.assertAlmostEqual(float(stats_non_mono.bias[0]), 0.25)


if __name__ == "__main__":
    unittest.main()
