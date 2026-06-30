# porescalemc

Pure-Python implementation of Multilevel Monte Carlo (MLMC) for uncertainty
quantification on random pore-scale geometries.

This is a modern reimplementation focused on a clean, testable Python core. It
generates hierarchical sphere and oriented ellipsoid packings, computes cheap
algebraic and spectral (FFT + Volume Penalisation) quantities of interest, runs
MLMC (fixed and adaptive), estimates computational speedups versus standard
Monte Carlo, and exports results as structured VTK files for ParaView
visualization.

The original 2013-2016 research code and paper are preserved for reference
(see `awr.tex` and `legacy/`). External PDE solvers (OpenFOAM/GetDP) are left
for future work; a Gmsh mesh-only adapter is present as an optional,
experimental bridge for future discrete solvers.

## Key Features

- **Geometry**: Deterministic, name-based generation of sphere/ellipsoid
  packings with Jodrey-Tory compaction and hierarchical refinement.
- **Cheap QoIs**: Algebraic packing statistics, analytical Fourier
  solid-fraction fields, structure/correlation statistics, and upscaled
  permeability.
- **Diffusion/permeability solvers**: Spectral VPM diffusion/Stokes/advection-
  diffusion on the Fourier grid, plus a voxel finite-volume diffusion solver
  for independent checks.
- **MLMC**: Giles-style estimator with paired sampling, adaptive convergence,
  optimal sample allocation, bias diagnostics, and level-aware solver factories.
- **Diagnostics**: Observed finite-level mean, variance, work, sample counts,
  and alpha/beta/gamma fits without relying on asymptotic assumptions.
- **Post-processing**: XML VTK ImageData (`.vti`) export and quick slice plots
  for solid fraction, porosity, correctors, and velocity fields.
- **Optional Gmsh meshing**: Soft-dependency OCC mesher, mesh-QoI solver, and
  tetrahedral mesh-backed diffusion estimate for conformal tetrahedral meshes,
  including rotated ellipsoids.

## Installation

```bash
pip install -e .
pip install -e ".[plot,test]"
# optional mesh-only Gmsh bridge
pip install -e ".[mesh]"
```

Requires Python >= 3.10, NumPy, and SciPy.

The package uses the standard `src/porescalemc/` layout.

## Quick Start

```python
from porescalemc import run_mlmc
from porescalemc.config import PackingConfig, MLMCConfig

pcfg = PackingConfig(mu=0.08, n_grains=25, periodic=True)
mcfg = MLMCConfig(n_levels=3, tolerance=0.05, algorithm="fixed")

mean, err = run_mlmc(pcfg, mcfg)
print(mean, err)
```

## Spectral Solver Example

```python
from porescalemc.config import PackingConfig, MLMCConfig, SpectralConfig, spectral_hierarchy_level
from porescalemc.geometry.placement import sample_packing_from_name
from porescalemc.mlmc.estimator import MLMCEstimator
from porescalemc.solvers.spectral import SpectralDiffusionSolver

def geo_factory(level, name, pcfg):
    return sample_packing_from_name(level, name, pcfg)

def solver_factory(packing, level=0):
    cfg = spectral_hierarchy_level(SpectralConfig(resolution=12, eta=1e-4), level)
    return SpectralDiffusionSolver(spectral_config=cfg)

est = MLMCEstimator(
    name="diffusion-uq",
    packing_config=PackingConfig(n_grains=12),
    mlmc_config=MLMCConfig(n_levels=3, m0=8),
    geometry_factory=geo_factory,
    solver_factory=solver_factory,
)
mean, err = est.run()
diagnostics = est.diagnostics()
```

## Solver and QoI Discovery

```python
from porescalemc.solvers import available_solvers, available_qois, solver_class

print(available_solvers())
print(available_qois("packing"))
print(available_qois("voxel_diffusion"))

VoxelDiffusionSolver = solver_class("voxel_diffusion")
```

Built-in solver keys currently include:

- `packing`: algebraic packing statistics.
- `fourier`: Fourier-field porosity/permeability closures.
- `spectral_diffusion`, `spectral_stokes`, `spectral_advection_diffusion`.
- `voxel_diffusion`: periodic voxel finite-volume diffusion cell solve.
- `tet_diffusion`: optional Gmsh tetrahedral P1-FEM diffusion estimate.

## VTK Export

```python
from porescalemc.config import FourierConfig
from porescalemc.geometry.fourier_field import porosity_field, save_structured_vti

fcfg = FourierConfig(resolution=24, smoothing="cell")
phi_s = porosity_field(packing, fcfg)
origin = (-0.5, -0.5, -0.5)
spacing = (1.0 / phi_s.shape[0],) * 3

save_structured_vti(phi_s, origin, spacing, "solid_fraction.vti", "solid_fraction")
```

## Examples

- `examples/pure_python_mlmc.py`: basic packing-statistics MLMC.
- `examples/spectral_mlmc.py`: spectral diffusion inside MLMC with diagnostic
  plots.
- `examples/spectral_solvers.py`: Hasimoto validation, diffusion/Stokes,
  advection-diffusion, VTI export, and speedup reporting.

Run after `pip install -e .`:

```bash
python examples/spectral_solvers.py
```

Or from a clean source tree without editable install:

```bash
PYTHONPATH=/path/to/porescalemc/src python examples/spectral_solvers.py
```

## Testing

```bash
python -m unittest discover -v
# or, if pytest is installed:
python -m pytest tests/ -q
```

Tests that actually invoke Gmsh mesh generation are skipped by default even
when the `gmsh` Python package is installed. To run them explicitly:

```bash
PORESCALEMC_RUN_GMSH_TESTS=1 python -m unittest tests.test_mesher -v
```

No OpenFOAM, GetDP, or Gmsh execution is required for the default test suite.

## References

- Icardi et al., "On the predictivity of pore-scale simulations: an
  Uncertainty Quantification approach", *Advances in Water Resources* (2016).
- MLMC foundations: Giles (2008, 2015), Cliffe et al., Haji-Ali et al.

See `porescalemc.tex` for detailed mathematical documentation.

## License

GPL (see `LICENSE`). Please cite the 2016 paper when using or building on this
work.

## Status

Pure-Python core is complete and tested: geometry, Fourier fields, spectral VPM
solvers, voxel diffusion, MLMC engine, post-processing, and speedup analysis.
The Gmsh mesher/tet diffusion path is an optional experimental bridge; legacy
OpenFOAM/GetDP solver paths remain reference material for future development.
