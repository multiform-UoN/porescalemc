#!/usr/bin/env python3
"""
Pure-Python MLMC example using the new porescalemc package.

Demonstrates:
- Deterministic packing generation via name
- PackingStatsSolver (algebraic QoIs)
- Fixed and adaptive MLMC
- Same realization used for coarse/fine levels (key for variance reduction)

Run with:
    PYTHONPATH=. python examples/pure_python_mlmc.py
"""

from __future__ import annotations

import numpy as np

from porescalemc.config import PackingConfig, MLMCConfig, FourierConfig
from porescalemc.geometry.placement import sample_packing_from_name
from porescalemc.solvers.dummy import DummySolver
from porescalemc.mlmc.estimator import MLMCEstimator
from porescalemc.solvers.packing import PackingStatsSolver


def main():
    print("=== porescalemc Pure-Python MLMC Demo ===\n")

    # Base packing config (level 0)
    pcfg = PackingConfig(
        mu=0.08,
        n_grains=25,
        coeff_var=0.15,
        psd="lognormal",
        periodic=True,
        detached=True,
        min_porosity=0.35,
        use_jodrey_tory=True,
        jt_min_dist=1.0,
    )

    # MLMC config
    mcfg = MLMCConfig(
        n_levels=3,
        m0=8,
        mratio=4.0,
        min_samples=4,
        tolerance=0.08,
        error_split=0.5,
        n_workers=1,           # set >1 for parallel (requires picklable objects)
        algorithm="fixed",     # or "adaptive"
        hierarchy="g",         # 'g' = grid refinement (affects future Fourier levels)
        reuse_samples=True,
    )

    # Common geometry factory (name-based for correlated pairs)
    def geo_factory(level, name, cfg):
        return sample_packing_from_name(level, name, cfg)

    # Use the built-in algebraic packing solver (8 QoIs by default)
    def solver_factory(packing):
        return PackingStatsSolver(packing, qoi_names=(
            "porosity", "solid_fraction", "specific_surface_area",
            "mean_equivalent_radius"
        ))

    # Identity QoI (solver already returns the vector we want)
    def qoi_fn(raw, packing):
        return raw

    est = MLMCEstimator(
        name="demo_packing",
        packing_config=pcfg,
        mlmc_config=mcfg,
        geometry_factory=geo_factory,
        solver_factory=solver_factory,
        qoi_fn=qoi_fn,
    )

    mean, error = est.run()
    print("\nFixed MLMC result:")
    print("  Mean (QoIs):", mean)
    print("  Error bound:", error)
    print("  Stats summary:\n" + est.stats.summary())

    # Quick demonstration of adaptive
    mcfg2 = MLMCConfig(
        n_levels=2,
        m0=6,
        tolerance=0.15,
        algorithm="adaptive",
        max_levels=4,
    )
    est2 = MLMCEstimator(
        name="demo_adaptive",
        packing_config=pcfg,
        mlmc_config=mcfg2,
        geometry_factory=geo_factory,
        solver_factory=solver_factory,
    )
    mean2, err2 = est2.run()
    print("\nAdaptive MLMC result:")
    print("  Mean:", mean2)
    print("  Error:", err2)
    print("  Final samples per level:", [len(s) for s in est2._samples])

    print("\n=== Done. All pure Python, no external solvers. ===")


if __name__ == "__main__":
    main()
