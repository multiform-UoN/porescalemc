# porescalemc

Multilevel Monte Carlo (MLMC) for uncertainty quantification of pore-scale effective
properties in random sphere/ellipsoid packings — pure Python, no external PDE solvers
required for the core workflow.

## Citation

If you use this software, please cite:

> Icardi et al. (2016), *Pore-scale simulation of fluid flow and solute dispersion in
> three-dimensional porous media*, Advances in Water Resources.
> DOI: [10.1016/j.advwatres.2016.01.004](https://doi.org/10.1016/j.advwatres.2016.01.004)

The accompanying paper describing the methods is `porescalemc.tex` / `porescalemc.pdf`.

## Features

- **Random packings**: Random Sequential Addition of spheres and oriented ellipsoids,
  with optional Jodrey-Tory compaction for higher solid fractions.
- **Analytical Fourier field** (novel): exact k-space representation of each grain
  (sphere or ellipsoid) without voxelisation. Delivers the smooth porosity field,
  structure factor, two-point autocorrelation, and correlation length purely from
  algebra.
- **Porosity-permeability laws**: Kozeny-Carman, Ergun (with Reynolds correction),
  Wen-Yu; upscaled via harmonic mean of the Fourier-derived porosity field.
- **MLMC estimator**: Giles (2008/2015) algorithm with fixed or adaptive sample
  allocation, deterministic name-based pairing for variance reduction, and Richardson
  extrapolation for bias correction.
- **Typed configuration**: all parameters in `PackingConfig`, `MLMCConfig`,
  `FourierConfig` dataclasses — no global state, no `exec()`.

## Quick start

```python
from porescalemc import run_mlmc
from porescalemc.config import PackingConfig, MLMCConfig

mean, error = run_mlmc(
    PackingConfig(mu=0.08, n_grains=25),
    MLMCConfig(n_levels=3, m0=20, tolerance=0.05),
)
print(f"Porosity: {mean[0]:.4f} ± {error[0]:.4f}")
```

Or run the bundled example:

```bash
PYTHONPATH=. python3 examples/pure_python_mlmc.py
```

## Installation

```bash
pip install .
# with optional plotting:
pip install ".[plot]"
# with test dependencies:
pip install ".[test]"
```

Requires Python ≥ 3.10, NumPy ≥ 1.24, SciPy ≥ 1.10.

## Package layout

```
porescalemc/
├── config.py              # PackingConfig, MLMCConfig, FourierConfig, hierarchy_level
├── geometry/
│   ├── grains.py          # Grain, Packing, stretching_length, ellipsoid_distance
│   ├── placement.py       # Random Sequential Addition, seed_from_name
│   ├── fourier_field.py   # Analytical FT, porosity field, structure factor, autocorrelation
│   ├── transforms.py      # Kozeny-Carman, Ergun, Wen-Yu permeability laws
│   └── jodrey_tory.py     # Jodrey-Tory compaction algorithm
├── mlmc/
│   ├── statistics.py      # compute_mlmc_stats, optimal_sample_counts, check_convergence
│   ├── estimator.py       # MLMCEstimator, PairSampler
│   └── workers.py         # Parallel sample execution (concurrent.futures)
└── solvers/
    ├── base.py            # SolverProtocol (typing.Protocol)
    ├── dummy.py           # DummySolver (returns zeros)
    ├── packing.py         # PackingStatsSolver (8 algebraic QoIs)
    └── fourier.py         # FourierSolver (porosity + permeability from Fourier field)
```

## Fourier field — mathematical background

The key innovation is the analytical Fourier transform of each grain, avoiding
voxelisation:

**Sphere** (radius *r*, centre **x**₀):

    F(k) = V · exp(-2πi k·x₀) · G(2π|k|r)

**Ellipsoid** (transformation matrix **M**, centre **x**₀):

    F(k) = |det M| · (4π/3) · exp(-2πi k·x₀) · G(2π|M^T k|)

where `G(ρ) = 3[sin(ρ) - ρ cos(ρ)] / ρ³` with `G(0) = 1`.

The porosity field is reconstructed as:

    φ(x) = Re[ IFFT( Σᵢ Fᵢ(k) · K(k) ) ] / V_domain

where `K(k)` is an optional smoothing kernel (`'cell'`, `'gaussian'`, or `'none'`).

## MLMC hierarchy levels

Each MLMC level *l* uses a scaled `PackingConfig` determined by `MLMCConfig.hierarchy`:

| Flag  | Effect at level *l*                                          |
|-------|--------------------------------------------------------------|
| `'g'` | Fourier grid resolution × `refratio^l`                       |
| `'d'` | Domain length × `refratio^l`                                 |
| `'s'` | Mean grain radius ÷ `refratio^l`                             |
| `'n'` | Grain count × `refratio^(dim·l)` (constant packing fraction) |

`refratio` (spatial refinement, default 2.0) and `mratio` (sample-count ratio,
default 4.0) are **independent** parameters.

## Running tests

```bash
python3 -m pytest tests/ -v
```

All 27 tests pass. No external solvers (OpenFOAM, gmsh, GetDP) are needed.

## Legacy code

The original 2016 code is preserved in `code/` for historical reference. It
requires OpenFOAM ≥ 2.3, GMSH, and GetDP. The new `porescalemc/` package
replaces it for pure-Python workflows; the external-solver interface is retained
as `SolverProtocol` for future extension.

## License

GPL — see `LICENSE`.

## Authors

Matteo Icardi (with G. Boccardo, H. Hoel, N. Quadrio)
November 2013 – June 2026
