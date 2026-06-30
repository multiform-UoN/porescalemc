"""Example demonstrating MLMCEstimator using the new SpectralDiffusionSolver."""

from __future__ import annotations

import numpy as np

from porescalemc.config import (
    FourierConfig,
    MLMCConfig,
    PackingConfig,
    SpectralConfig,
    spectral_hierarchy_level,
)
from porescalemc.mlmc.estimator import MLMCEstimator
from porescalemc.solvers.spectral import SpectralDiffusionSolver


def main() -> None:
    print("=== porescalemc Spectral MLMC Demo ===")

    # 1. Setup configurations
    # We use a small packing size for demonstration speed
    pcfg = PackingConfig(
        n_grains=10,
        psd="constant",
        mu=0.08,
        xlen=1.0,
        ylen=1.0,
        zlen=1.0,
        periodic=True,
    )

    # MLMC config: 2 levels, 8 coarse samples, 3 minimum samples
    mlmc_cfg = MLMCConfig(
        n_levels=2,
        m0=8,
        min_samples=3,
        tolerance=0.05,
        n_workers=1,
        hierarchy="g",
    )

    # 2. Spectral/Fourier resolution parameters
    spectral_cfg = SpectralConfig(
        resolution=10,
        eta=1e-4,
        max_iter=100,
        tol=1e-4,
        n_directions=1,
        resolution_level_factor=1.5,
    )

    # 3. Solver factory returning the VPM Diffusion solver
    def solver_factory(packing, level: int = 0) -> SpectralDiffusionSolver:
        cfg_level = spectral_hierarchy_level(spectral_cfg, level)
        fourier_cfg = FourierConfig(
            resolution=cfg_level.resolution,
            smoothing=cfg_level.smoothing,
            smoothing_length=cfg_level.smoothing_length,
        )
        return SpectralDiffusionSolver(
            packing=packing,
            spectral_config=cfg_level,
            fourier_config=fourier_cfg,
        )

    def identity_qoi(output, packing) -> np.ndarray:
        del packing
        return output

    from porescalemc.geometry.placement import sample_packing_from_name
    def geo_factory(level, name, cfg):
        return sample_packing_from_name(level, name, cfg)

    print("\nRunning MLMCEstimator with SpectralDiffusionSolver...")
    estimator = MLMCEstimator(
        name="spectral-diffusion-uq",
        packing_config=pcfg,
        mlmc_config=mlmc_cfg,
        geometry_factory=geo_factory,
        solver_factory=solver_factory,
        qoi_fn=identity_qoi,
    )

    mean, error = estimator.run()

    print("\nMLMC Result:")
    print(f"  Mean effective diffusivity:       {mean[0]:.6f}")
    print(f"  Error Bound (Bias + Stat):        {error[0]:.6f}")
    print(f"  Samples drawn per level:          {estimator.stats.n_samples}")
    print(f"  Average work per sample (s):      {estimator.stats.work}")

    diagnostics = estimator.diagnostics()
    print("  Observed rates:")
    for name, value in diagnostics["observed_rates"].items():
        print(f"    {name}: {value:.3g}")
    print("  Level table:")
    for row in diagnostics["levels"]:
        print(
            "    level {level:.0f}: mean={mean_abs_max:.3e}, "
            "var={variance_max:.3e}, work={work_per_sample:.3e}, samples={n_samples:.0f}"
            .format(**row)
        )

    estimator.plot_convergence("spectral_mlmc_diagnostics.png")
    print("  Wrote spectral_mlmc_diagnostics.png if matplotlib is available.")

    print("\n=== Done. All pure Python spectral solvers under VPM. ===")


if __name__ == "__main__":
    main()
