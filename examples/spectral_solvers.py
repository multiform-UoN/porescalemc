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

from porescalemc.config import PackingConfig, SpectralConfig, spectral_hierarchy_level
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
print(f"  phi_solid          = {packing_rand.solid_fraction:.4f}")
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
