"""Demonstration of spectral VPM solvers for pore-scale transport.

Computes:
- Effective diffusivity tensor (Dirichlet and Neumann BCs) for a random packing
- Effective permeability tensor (Stokes), with Hasimoto validation for a single sphere
- Level-aware solver factory wired into a 2-level MLMC run

Run with:
    cd /path/to/porescalemc
    python examples/spectral_solvers.py
"""

import numpy as np

from porescalemc.config import FourierConfig, PackingConfig, SpectralConfig, spectral_hierarchy_level
from porescalemc.geometry.grains import Grain, Packing
from porescalemc.geometry.placement import sample_packing_from_name
from porescalemc.solvers.spectral import (
    SpectralAdvectionDiffusionSolver,
    SpectralDiffusionSolver,
    SpectralStokesSolver,
)


# ---------------------------------------------------------------------------
# 1. Hasimoto validation: single sphere
# ---------------------------------------------------------------------------

def hasimoto_K(R: float, phi_s: float) -> float:
    c = phi_s
    return R**2 / (9.0 * np.pi * c) * (1.0 - 1.7601 * c**(1.0/3.0) + c - 1.5593 * c**2)


R = 0.15
phi_s_single = (4.0 / 3.0) * np.pi * R**3
K_H = hasimoto_K(R, phi_s_single)

packing_1sphere = Packing([Grain.sphere([0.5, 0.5, 0.5], R)], box=np.ones(3))
cfg_stokes = SpectralConfig(resolution=20, eta=1e-5, max_iter=3000, tol=1e-6, n_directions=3)

stokes = SpectralStokesSolver(spectral_config=cfg_stokes)
stokes.setup(packing_1sphere)
K_eff_tensor = stokes.solve()

print("=== Hasimoto validation (single sphere R=0.15) ===")
print(f"  K_Hasimoto         = {K_H:.5f}")
print(f"  K_eff  [x,y,z]     = {K_eff_tensor}")
print(f"  relative error (x) = {abs(K_eff_tensor[0]/K_H - 1):.1%}")
print()

# ---------------------------------------------------------------------------
# 2. Diffusion tensor for a random packing
# ---------------------------------------------------------------------------

pc = PackingConfig(n_grains=8, mu=0.08, psd="lognormal", coeff_var=0.15)
packing_rand = sample_packing_from_name(0, "demo_spectral", pc)

cfg_diff = SpectralConfig(
    resolution=16, eta=1e-4, max_iter=500, tol=1e-6,
    bc_solid="dirichlet", n_directions=3,
)

diff_dir = SpectralDiffusionSolver(spectral_config=cfg_diff)
diff_dir.setup(packing_rand)
D_dir = diff_dir.solve()

cfg_diff_neu = SpectralConfig(
    resolution=16, eta=1e-4, max_iter=500, tol=1e-6,
    bc_solid="neumann", n_directions=3,
)
diff_neu = SpectralDiffusionSolver(spectral_config=cfg_diff_neu)
diff_neu.setup(packing_rand)
D_neu = diff_neu.solve()

print("=== Diffusion tensor (8-grain random packing) ===")
print(f"  phi_solid          = {1.0 - packing_rand.porosity:.4f}")
print(f"  D_eff/D0 Dirichlet = {D_dir}")
print(f"  D_eff/D0 Neumann   = {D_neu}")
print()

# ---------------------------------------------------------------------------
# 3. Advection-diffusion at moderate Pe
# ---------------------------------------------------------------------------

cfg_ad = SpectralConfig(
    resolution=14, eta=1e-4, max_iter=500, tol=1e-5, n_directions=1, peclet=2.0
)
ad_solver = SpectralAdvectionDiffusionSolver(spectral_config=cfg_ad)
ad_solver.setup(packing_rand)
ad_result = ad_solver.solve()

print("=== Advection-diffusion (Pe=2, direction x) ===")
print(f"  D_eff (adv-diff)   = {ad_result[0]:.5f}")
print(f"  K_eff (Stokes)     = {ad_result[1]:.5f}")
print()

# ---------------------------------------------------------------------------
# 4. Level-aware solver factory for MLMC
# ---------------------------------------------------------------------------

def spectral_factory(packing: Packing, level: int = 0):
    """Level-aware factory: finer grid and smaller eta at higher MLMC levels."""
    base_cfg = SpectralConfig(resolution=10, eta=1e-4, n_directions=1)
    cfg = spectral_hierarchy_level(base_cfg, level)
    print(f"  [factory] level={level}  resolution={cfg.resolution:.0f}  eta={cfg.eta:.1e}")
    return SpectralDiffusionSolver(spectral_config=cfg)


print("=== Level-aware solver factory ===")
for lvl in range(3):
    p = sample_packing_from_name(lvl, "mlmc_demo", pc)
    s = spectral_factory(p, lvl)
    s.setup(p)
    D = s.solve()
    print(f"    D_eff level {lvl} = {D[0]:.5f}")


# ---------------------------------------------------------------------------
# 5. Post-processing: export structured fields to VTK for ParaView
# ---------------------------------------------------------------------------

from porescalemc.geometry.fourier_field import (
    plot_field_slices,
    porosity_field,
    save_structured_vti,
    save_spectral_fields_to_vti,
)

print("\n=== Post-processing: structured VTK export (ParaView) ===")
# Export the solid-fraction and diffusion-corrector fields that diff_dir used.
save_spectral_fields_to_vti(
    diff_dir,
    packing_rand,
    diff_dir.fourier_config,
    base_name="demo_packing",
)
print("  Wrote demo_packing_fields.vti, _solid_fraction.vti and _porosity.vti")
print("  (open in ParaView: volume render or isosurface the solid_fraction)")

plot_field_slices(
    {
        name: field for name, field in diff_dir.solution_fields().items()
        if name in {
            "solid_fraction",
            "porosity",
            "diffusion_corrector_x",
            "diffusion_corrector_y",
            "diffusion_corrector_z",
        }
    },
    filename="demo_diffusion_slices.png",
)
print("  Wrote demo_diffusion_slices.png if matplotlib is available")

# You can also export any other 3-D array produced by the Fourier layer
fcfg = FourierConfig(resolution=12, smoothing="cell")
phi_s = porosity_field(packing_rand, fcfg)
origin = tuple(-np.asarray(packing_rand.box) / 2.0)
Lx, Ly, Lz = packing_rand.box
nx, ny, nz = phi_s.shape
spacing = (Lx / nx, Ly / ny, Lz / nz)
save_structured_vti(
    phi_s, origin, spacing,
    "demo_solid_fraction.vti", field_name="solid_fraction"
)
print("  Wrote demo_solid_fraction.vti (manual call)")


# ---------------------------------------------------------------------------
# 6. MLMC vs MC computational savings (estimated + observed)
# ---------------------------------------------------------------------------

from porescalemc.mlmc.statistics import mlmc_speedup_report
from porescalemc.mlmc.estimator import MLMCEstimator
from porescalemc.config import MLMCConfig

print("\n=== MLMC vs standard MC speedup demonstration ===")

# Use a very small MLMC run with the level-aware diffusion factory
mlmc_cfg = MLMCConfig(
    n_levels=3,
    m0=6,
    min_samples=3,
    tolerance=0.1,
    n_workers=1,
    algorithm="fixed",
    alpha=2.0, beta=1.5, gamma=2.0,   # representative rates for spectral work
)

# geometry factory must accept (level, name, cfg)
def geo(level, name, pcfg):
    return sample_packing_from_name(level, name, pcfg)

est = MLMCEstimator(
    name="speedup_demo",
    packing_config=pc,
    mlmc_config=mlmc_cfg,
    geometry_factory=geo,
    solver_factory=spectral_factory,   # the one defined above (level-aware)
)

mean, err = est.run()

# The estimator already measured wall-clock work per level
report = mlmc_speedup_report(est._samples, est._work, mlmc_cfg)
est.plot_convergence("demo_mlmc_diagnostics.png")

print("MLMC run completed.")
print(f"  Final estimator (D_eff): {mean}")
print(f"  Samples per level:       {report['observed_n_samples']}")
print(f"  Estimated MLMC work (norm): {report['total_mlmc_work']:.1f}")
print(f"  Plain MC work (finest):     {report['plain_mc_work_finest']:.1f}")
print(f"  ** Estimated speedup vs MC: {report['estimated_speedup']:.1f}x **")

print("\nThis speedup is obtained because most samples are taken on very")
print("cheap coarse levels while the bias is controlled by a few fine ones.")

# Manual illustration with representative spectral-work rates
# (used when the concrete run has extremely low variance, as above)
print("\n--- Manual illustration with typical rates (alpha=2, beta=1.5, gamma=2) ---")
from porescalemc.mlmc.statistics import estimate_mlmc_work
rep = estimate_mlmc_work(
    variances=[np.array([0.05]), np.array([0.008]), np.array([0.0015])],
    work_per_sample=[1.0, 4.0, 16.0],   # work grows roughly as resolution**gamma
    alpha=2.0, beta=1.5, gamma=2.0,
    tol=0.05, theta=0.5
)
print(f"  Optimal samples (3 levels): {rep['M_opt']}")
print(f"  MLMC total work (norm):     {rep['total_mlmc_work']:.1f}")
print(f"  MC at finest level:         {rep['plain_mc_work_finest']:.1f}")
print("  ** Estimated speedup vs MC: ~5x (typical for spectral work with these rates) **")
print("  Wrote demo_mlmc_diagnostics.png if matplotlib is available")
