# Reviewer prompt for porescalemc

You are reviewing the Python package `porescalemc` — a Multilevel Monte Carlo
(MLMC) framework for uncertainty quantification of pore-scale effective
properties (permeability, diffusivity) via random sphere/ellipsoid packings.

The code is at the repo root.  **Do not run OpenFOAM, gmsh, or GetDP.**
Python-only tests are in `tests/`.

---

## 1. Run the test suite

```bash
python3 -m pytest tests/ -v
```

All 27 tests must pass.

---

## 2. Run the example

```bash
PYTHONPATH=. python3 examples/pure_python_mlmc.py
```

Expected: two MLMC runs (fixed and adaptive) print a 4-component QoI estimate
and error bound, then "=== Done. All pure Python, no external solvers. ===".

---

## 3. Things to exercise specifically

The novel scientific contribution is the **analytical Fourier-space grain
field** (`porescalemc/geometry/fourier_field.py`).  Check these properties:

| Property | How to test |
|----------|-------------|
| DC component = solid volume | `packing_fourier_transform(packing, N, N, N)[0][0,0,0].real ≈ packing.solid_volume` |
| Parseval theorem | `(1/(N³·V)) Σ|F_k|² ≈ solid_volume/V` (±0.01 for N=10) |
| Autocorrelation peak at zero lag | `autocorrelation(packing, cfg)[0,0,0] == 1.0` |
| Ellipsoid reduces to sphere for isotropic M | `grain_fourier_transform(sphere, …) == grain_fourier_transform(ellipsoid_with_M=r·I, …)` |
| Porosity field mean = 1 − porosity | `porosity_field(packing, cfg).mean() ≈ 1 − packing.porosity` (±0.001) |

---

## 4. MLMC math to verify

- **Estimator = sum of level means**: check `MLMCStats.estimator ≈ sum(MLMCStats.means)`.
- **Richardson extrapolation**: with `extrapolate_bias=True`, the correction added
  to the estimator should equal `means[-1] / (2^alpha − 1)` (not the *difference*
  of the last two level means).
- **Optimal sample counts**: with all variances equal and all work equal, the
  optimal allocation should be uniform across levels (`M_l* ≈ M_0` for all l).
- **Serialisation**: `estimator.save(path); est2 = MLMCEstimator.load(path)` should
  restore the same estimate, error, packing config, and mlmc config.

---

## 5. Config hierarchy independence

`MLMCConfig.mratio` (sample count ratio) and `MLMCConfig.refratio` (spatial
refinement ratio) are now independent fields (defaults 4.0 and 2.0).  Verify
that `hierarchy_level(mlmc_cfg, packing_cfg, level)` uses `refratio` for domain
/ grain-size / mesh scaling, and that changing `mratio` alone does not alter the
geometry at any level.

---

## 6. What the code does NOT do (by design)

- No mesh generation (gmsh) — `FourierSolver` computes permeability purely from
  the analytical FT without any mesh.
- No PDE solves (OpenFOAM/GetDP) — all QoIs come from algebraic laws
  (Kozeny-Carman, Ergun, Wen-Yu) or packing statistics.
- No global state / `exec(open(...))` — all configuration is typed dataclasses.
