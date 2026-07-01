# porescalemc

Pure-Python Multilevel Monte Carlo (MLMC) for uncertainty quantification on
random pore-scale geometries.

Generates sphere/ellipsoid packings, solves effective diffusivity and
permeability via spectral (FFT/VPM), voxel FV, or tetrahedral FEM methods,
and runs MLMC with Giles-optimal sample allocation and adaptive level
management.

---

## Installation

```bash
pip install -e .
pip install -e ".[plot,test]"    # matplotlib + pytest
pip install -e ".[mesh]"         # optional: Gmsh tetrahedral FEM
```

Requires Python ≥ 3.10, NumPy, SciPy.

---

## Quick start — one line

```python
from porescalemc import run_study

mean, err = run_study()
print(f"D_eff = {mean[0]:.4f} ± {err[0]:.4f}")
```

`run_study()` uses defaults: `spectral_diffusion` solver, 20 grains, 3 MLMC
levels.  Adjust any parameter as a keyword argument.

---

## Choosing solver, geometry, and MLMC settings

```python
from porescalemc import run_study, list_solvers, list_qois

# See all available solvers and their QoI names
list_solvers()

# Run a study with explicit settings
mean, err = run_study(
    # --- solver ---
    solver="spectral_diffusion",   # see list_solvers()
    resolution=20,                 # FFT grid points per unit length
    n_directions=1,                # 1 = x-only, 3 = full 3×3 tensor
    eta=1e-5,                      # VPM penalisation (smaller = sharper interface)

    # --- geometry ---
    n_grains=20,                   # grains per unit-volume box
    mu=0.08,                       # mean grain radius
    coeff_var=0.1,                 # radius spread (0 = monodisperse)
    box_size=1.0,                  # cubic box side length

    # --- MLMC ---
    n_levels=3,                    # coarse-to-fine levels
    m0=20,                         # samples at coarsest level
    tolerance=0.05,                # target relative RMSE
    hierarchy="g",                 # 'g'=grid, 'd'=domain, 'n'=grains, or combo
    algorithm="fixed",             # 'fixed' or 'adaptive'
    n_workers=1,                   # parallel workers (1 = serial)
)
```

### Available solvers

| Name | What it computes | Key params |
|------|-----------------|------------|
| `packing` | Porosity, surface area, packing statistics | — |
| `fourier` | Fourier-field porosity + permeability closure | `resolution` |
| `spectral_diffusion` | Effective diffusivity D_eff (VPM) | `resolution`, `eta`, `n_directions` |
| `spectral_stokes` | Effective permeability K_eff (VPM Stokes) | `resolution`, `eta`, `n_directions` |
| `spectral_advection_diffusion` | D_eff + K_eff at finite Pe | `resolution`, `eta`, `n_directions` |
| `voxel_diffusion` | D_eff from periodic voxel FV cell solve | `resolution`, `n_directions` |
| `tet_diffusion` | D_eff from periodic P1 FEM (requires gmsh) | mesh_size |

### QoI names

```python
from porescalemc import list_qois

list_qois("spectral_diffusion")
# spectral_diffusion: D_eff_x, D_eff_y, D_eff_z, D_tensor_3x3_if_n_directions_3

list_qois("voxel_diffusion")
# voxel_diffusion: diffusivity_x, diffusivity_y, diffusivity_z

list_qois()   # all solvers
```

---

## Common recipes

### Stokes permeability

```python
mean, err = run_study(
    solver="spectral_stokes",
    n_grains=15, mu=0.10,
    resolution=24, n_directions=1,
    n_levels=3, m0=16, tolerance=0.05,
)
print(f"K_eff = {mean[0]:.4e} ± {err[0]:.4e}")
```

### Full 3×3 diffusivity tensor

```python
mean, err = run_study(
    solver="spectral_diffusion",
    n_grains=20, mu=0.08,
    resolution=20, n_directions=3,   # returns 9 components
    n_levels=3, m0=20, tolerance=0.05,
)
D_tensor = mean.reshape(3, 3)
print(D_tensor)
```

### Voxel FV (independent check)

```python
mean, err = run_study(
    solver="voxel_diffusion",
    n_grains=20, mu=0.08,
    resolution=16, n_directions=1,
    n_levels=3, m0=16, tolerance=0.05,
)
```

### Adaptive MLMC (auto-level and sample allocation)

```python
mean, err = run_study(
    solver="spectral_diffusion",
    algorithm="adaptive",
    n_levels=2,       # starting levels; may add more
    m0=10, tolerance=0.03,
)
```

### Polydisperse packing

```python
mean, err = run_study(
    solver="spectral_diffusion",
    n_grains=25, mu=0.07,
    coeff_var=0.3,     # wide lognormal PSD
    psd="lognormal",
)
```

---

## Full config objects (advanced)

For fine-grained control, pass complete config objects.  These override
the keyword shortcuts above:

```python
from porescalemc import run_study
from porescalemc.config import PackingConfig, MLMCConfig, SpectralConfig

mean, err = run_study(
    solver="spectral_diffusion",
    packing_config=PackingConfig(
        n_grains=30, mu=0.07,
        coeff_var=0.2, psd="lognormal",
        ellipsoid=False,
        periodic=True,
    ),
    mlmc_config=MLMCConfig(
        n_levels=4, m0=32,
        tolerance=0.02, algorithm="fixed",
        hierarchy="g", refratio=2.0,
    ),
    solver_config=SpectralConfig(
        resolution=20, eta=1e-5,
        n_directions=3,
        max_iter=2000, tol=1e-7,
    ),
)
```

---

## MLMCEstimator (lowest-level interface)

For custom solver factories or non-standard geometries:

```python
from porescalemc import MLMCEstimator
from porescalemc.config import PackingConfig, MLMCConfig, SpectralConfig, spectral_hierarchy_level
from porescalemc.geometry.placement import sample_packing_from_name
from porescalemc.solvers.spectral import SpectralDiffusionSolver
import numpy as np

def geo_factory(level, name, pcfg):
    return sample_packing_from_name(level, name, pcfg)

def solver_factory(packing, level=0):
    cfg = spectral_hierarchy_level(SpectralConfig(resolution=12, eta=1e-4), level)
    return SpectralDiffusionSolver(spectral_config=cfg)

est = MLMCEstimator(
    name="diffusion-uq",
    packing_config=PackingConfig(n_grains=20, mu=0.08),
    mlmc_config=MLMCConfig(n_levels=3, m0=16, tolerance=0.05),
    geometry_factory=geo_factory,
    solver_factory=solver_factory,
)
mean, err = est.run()
print(est.diagnostics())         # per-level samples, means, variances, rates
```

The `solver_factory` can accept an optional `level` keyword argument
(detected via `inspect.signature`).  Use this to refine the mesh or lower
`eta` at finer MLMC levels.

---

## Geometry: PackingConfig options

| Field | Default | Meaning |
|-------|---------|---------|
| `n_grains` | 20 | Target number of grains |
| `mu` | 0.08 | Mean grain radius |
| `coeff_var` | 0.1 | CoV of radius PSD |
| `psd` | `"lognormal"` | `"lognormal"` / `"uniform"` / `"constant"` |
| `xlen`, `ylen`, `zlen` | 1.0 | Box side lengths |
| `periodic` | `True` | Wrap grains across box faces |
| `ellipsoid` | `False` | Use random-orientation ellipsoids |
| `detached` | `True` | Reject overlapping grain placements |
| `use_jodrey_tory` | `False` | Post-process to reduce overlaps |

---

## MLMC: MLMCConfig options

| Field | Default | Meaning |
|-------|---------|---------|
| `n_levels` | 3 | Initial number of levels |
| `m0` | 32 | Samples at coarsest level |
| `tolerance` | 0.01 | Target relative RMSE |
| `hierarchy` | `"g"` | Level scaling: `g`=grid, `d`=domain, `n`=grain count |
| `refratio` | 2.0 | Resolution doubling per level |
| `algorithm` | `"fixed"` | `"fixed"` or `"adaptive"` |
| `n_workers` | 1 | Parallel workers |
| `mratio` | 4.0 | Sample-count ratio between levels |
| `beta` | 1.2 | Assumed variance decay rate |
| `alpha` | 2.5 | Assumed bias decay rate |
| `extrapolate_bias` | `False` | Richardson bias correction |

---

## Diagnostics

```python
mean, err = run_study(solver="spectral_diffusion", ...)
# For detailed diagnostics, use MLMCEstimator directly:
est = MLMCEstimator(...)
mean, err = est.run()
diag = est.diagnostics()
# diag keys: level_means, level_variances, level_costs, sample_counts,
#            observed_alpha, observed_beta, observed_gamma, speedup_vs_mc
```

---

## VTK export (ParaView)

```python
from porescalemc.config import FourierConfig, SpectralConfig
from porescalemc.geometry.fourier_field import porosity_field, save_structured_vti
from porescalemc.geometry.fourier_field import save_spectral_fields_to_vti
from porescalemc.geometry.grains import Grain, Packing
from porescalemc.solvers.spectral import SpectralDiffusionSolver

packing = Packing([Grain.sphere([0.5, 0.5, 0.5], 0.15)], box=[1.0, 1.0, 1.0])
fcfg = FourierConfig(resolution=24)

# Save solid-fraction field
phi_s = porosity_field(packing, fcfg)
n = phi_s.shape[0]
save_structured_vti(phi_s, (-0.5,)*3, (1/n,)*3, "solid_fraction.vti", "solid_fraction")

# Save solver fields (solid fraction + correctors)
cfg = SpectralConfig(resolution=24, n_directions=1)
solver = SpectralDiffusionSolver(spectral_config=cfg)
solver.setup(packing)
solver.solve()
save_spectral_fields_to_vti(solver, packing, cfg, base_name="debug")
```

---

## Tetrahedral FEM (optional, requires gmsh)

```bash
pip install -e ".[mesh]"
```

```python
from porescalemc import run_study

mean, err = run_study(
    solver="tet_diffusion",
    n_grains=10, mu=0.10,
    n_levels=2, m0=5, tolerance=0.2,
)
```

Or directly:

```python
from porescalemc.geometry.grains import Grain, Packing
from porescalemc.solvers.tet import TetDiffusionSolver

packing = Packing([Grain.sphere([0.5, 0.5, 0.5], 0.2)], box=[1.0, 1.0, 1.0])
solver = TetDiffusionSolver(mesh_size=0.08, periodic=True)
solver.setup(packing)
result = solver.solve()   # [D_xx, D_yy, D_zz, porosity, n_elements]
solver.close()
```

`periodic=True` (default) solves the full homogenisation cell problem with
DOF condensation on periodic boundary nodes.  `periodic=False` uses a
simpler unit-gradient Dirichlet test.

---

## Testing

```bash
pip install -e ".[test]"
python -m pytest tests/ -q                   # 78 tests, ~5 s, no heavy deps
PORESCALEMC_RUN_GMSH_TESTS=1 python -m pytest tests/test_mesher.py -v
```

---

## References

- Icardi et al., *Advances in Water Resources* (2016) — original research.
- Giles (2008, 2015) — MLMC foundations.

See `awr.tex` for mathematical documentation of the spectral VPM and
homogenisation formulations.

## License

GPL.  Cite the 2016 paper when building on this work.
