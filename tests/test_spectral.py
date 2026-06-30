"""Tests for spectral VPM solvers.

Validation cases:
1. Hasimoto (1959) permeability formula for dilute SC sphere array.
2. Effective diffusivity equals D0 for an empty domain (phi_s = 0).
3. Effective diffusivity is less than D0 for a porous medium (tortuosity).
4. Stokes permeability is positive and finite for a porous medium.
5. Advection-diffusion at Pe=0 matches pure diffusion.
6. n_directions=3 returns shape (3,) for both diffusion and Stokes.
"""

import unittest
import os
import tempfile

import numpy as np

from porescalemc.config import SpectralConfig
from porescalemc.geometry.fourier_field import save_spectral_fields_to_vti
from porescalemc.geometry.grains import Grain, Packing
from porescalemc.solvers.spectral import (
    SpectralAdvectionDiffusionSolver,
    SpectralDiffusionSolver,
    SpectralStokesSolver,
)


def _sphere_packing(radius: float, center=(0.5, 0.5, 0.5)) -> Packing:
    return Packing([Grain.sphere(list(center), radius)], box=np.ones(3))


def _empty_packing() -> Packing:
    return Packing([], box=np.ones(3))


def _hasimoto_permeability(R: float, phi_s: float) -> float:
    """Hasimoto (1959) dilute SC sphere array permeability."""
    c = phi_s
    return R**2 / (9.0 * np.pi * c) * (
        1.0 - 1.7601 * c**(1.0 / 3.0) + c - 1.5593 * c**2
    )


class TestSpectralDiffusion(unittest.TestCase):

    def test_empty_domain_d_eff_equals_d0(self):
        """For phi_s = 0 (no grains) the cell problem RHS is zero → D_eff = D0."""
        cfg = SpectralConfig(resolution=8, eta=1e-4, max_iter=50, tol=1e-6, n_directions=1)
        packing = _empty_packing()
        solver = SpectralDiffusionSolver(spectral_config=cfg)
        solver.setup(packing)
        result = solver.solve()
        self.assertEqual(result.shape, (1,))
        self.assertAlmostEqual(float(result[0]), 1.0, delta=0.01)

    def test_porous_d_eff_less_than_d0(self):
        """A single sphere reduces D_eff below D0 (tortuosity effect)."""
        cfg = SpectralConfig(resolution=12, eta=1e-4, max_iter=200, tol=1e-6, n_directions=1)
        packing = _sphere_packing(radius=0.12)
        solver = SpectralDiffusionSolver(spectral_config=cfg)
        solver.setup(packing)
        result = solver.solve()
        self.assertEqual(result.shape, (1,))
        self.assertGreater(float(result[0]), 0.0)
        self.assertLess(float(result[0]), 1.0)

    def test_neumann_bc_positive(self):
        """Neumann (no-flux) BC also returns positive D_eff for a porous medium."""
        cfg = SpectralConfig(
            resolution=12, eta=1e-4, max_iter=200, tol=1e-6,
            bc_solid="neumann", n_directions=1,
        )
        packing = _sphere_packing(radius=0.12)
        solver = SpectralDiffusionSolver(spectral_config=cfg)
        solver.setup(packing)
        result = solver.solve()
        self.assertGreater(float(result[0]), 0.0)

    def test_n_directions_3_returns_tensor(self):
        """n_directions=3 returns shape (3,) with all positive entries."""
        cfg = SpectralConfig(resolution=10, eta=1e-4, max_iter=100, tol=1e-5, n_directions=3)
        packing = _sphere_packing(radius=0.1)
        solver = SpectralDiffusionSolver(spectral_config=cfg)
        solver.setup(packing)
        result = solver.solve()
        self.assertEqual(result.shape, (3,))
        self.assertTrue((result > 0).all())

    def test_solution_fields_and_diagnostics_are_available(self):
        """Diffusion solver should cache fields and scalar diagnostics for debugging."""
        cfg = SpectralConfig(resolution=8, eta=1e-4, max_iter=80, tol=1e-5, n_directions=1)
        packing = _sphere_packing(radius=0.1)
        solver = SpectralDiffusionSolver(spectral_config=cfg)
        solver.setup(packing)
        result = solver.solve()

        fields = solver.solution_fields()
        diagnostics = solver.diagnostics()

        self.assertIn("solid_fraction", fields)
        self.assertIn("porosity", fields)
        self.assertIn("diffusion_corrector_x", fields)
        self.assertEqual(fields["solid_fraction"].shape, fields["porosity"].shape)
        self.assertAlmostEqual(
            diagnostics["diffusivity_x"],
            float(result[0]),
            delta=1e-12,
        )
        self.assertGreaterEqual(diagnostics["porosity_mean"], 0.0)
        self.assertLessEqual(diagnostics["porosity_mean"], 1.0)

    def test_save_spectral_fields_to_vti_exports_combined_xml(self):
        """Convenience VTI export should include geometry and PDE fields."""
        cfg = SpectralConfig(resolution=6, eta=1e-4, max_iter=50, tol=1e-5, n_directions=1)
        packing = _sphere_packing(radius=0.1)
        solver = SpectralDiffusionSolver(spectral_config=cfg)
        solver.setup(packing)
        solver.solve()

        with tempfile.TemporaryDirectory() as tmpdir:
            base = os.path.join(tmpdir, "debug")
            save_spectral_fields_to_vti(
                solver,
                packing,
                solver.fourier_config,
                base_name=base,
            )
            with open(base + "_fields.vti", "r", encoding="utf-8") as f:
                text = f.read()
            self.assertIn('<VTKFile type="ImageData"', text)
            self.assertIn('Name="solid_fraction"', text)
            self.assertIn('Name="porosity"', text)
            self.assertIn('Name="diffusion_corrector_x"', text)
            self.assertTrue(os.path.exists(base + "_solid_fraction.vti"))
            self.assertTrue(os.path.exists(base + "_porosity.vti"))


class TestSpectralStokes(unittest.TestCase):

    def test_hasimoto_validation(self):
        """K_eff for a single sphere should agree with Hasimoto within 50%.

        The VPM introduces O(sqrt(eta)) error and Gibbs ringing at limited
        resolution; 50% is a realistic tolerance for res=20.
        """
        R = 0.15
        phi_s = (4.0 / 3.0) * np.pi * R**3  # single sphere in unit box
        K_H = _hasimoto_permeability(R, phi_s)

        cfg = SpectralConfig(resolution=20, eta=1e-5, max_iter=3000, tol=1e-6, n_directions=1)
        packing = _sphere_packing(radius=R)
        solver = SpectralStokesSolver(spectral_config=cfg)
        solver.setup(packing)
        K_eff = float(solver.solve()[0])

        rel_err = abs(K_eff / K_H - 1.0)
        self.assertLess(
            rel_err, 0.5,
            f"K_eff={K_eff:.4f}, K_Hasimoto={K_H:.4f}, rel_err={rel_err:.2%}",
        )

    def test_porous_positive_finite(self):
        """K_eff must be strictly positive and finite for any porous medium."""
        cfg = SpectralConfig(resolution=12, eta=1e-4, max_iter=500, tol=1e-5, n_directions=1)
        packing = _sphere_packing(radius=0.12)
        solver = SpectralStokesSolver(spectral_config=cfg)
        solver.setup(packing)
        K_eff = float(solver.solve()[0])
        self.assertGreater(K_eff, 0.0)
        self.assertTrue(np.isfinite(K_eff))

    def test_n_directions_3_returns_tensor(self):
        """n_directions=3 returns shape (3,) with all positive entries."""
        cfg = SpectralConfig(resolution=10, eta=1e-4, max_iter=200, tol=1e-5, n_directions=3)
        packing = _sphere_packing(radius=0.1)
        solver = SpectralStokesSolver(spectral_config=cfg)
        solver.setup(packing)
        result = solver.solve()
        self.assertEqual(result.shape, (3,))
        self.assertTrue((result > 0).all())


class TestSpectralAdvDiff(unittest.TestCase):

    def test_zero_peclet_matches_diffusion(self):
        """At Pe=0, SpectralAdvectionDiffusionSolver D_eff must match SpectralDiffusionSolver."""
        cfg = SpectralConfig(
            resolution=10, eta=1e-4, max_iter=200, tol=1e-5, n_directions=1, peclet=0.0
        )
        packing = _sphere_packing(radius=0.1)

        ad_solver = SpectralAdvectionDiffusionSolver(spectral_config=cfg)
        ad_solver.setup(packing)
        ad_result = ad_solver.solve()
        D_eff_ad = float(ad_result[0])

        diff_solver = SpectralDiffusionSolver(spectral_config=cfg)
        diff_solver.setup(packing)
        D_eff_pure = float(diff_solver.solve()[0])

        self.assertAlmostEqual(D_eff_ad, D_eff_pure, delta=0.05)

    def test_output_shape(self):
        """SpectralAdvectionDiffusionSolver returns shape (2*n_directions,)."""
        cfg = SpectralConfig(resolution=8, eta=1e-4, max_iter=100, tol=1e-5, n_directions=2, peclet=0.0)
        packing = _sphere_packing(radius=0.1)
        solver = SpectralAdvectionDiffusionSolver(spectral_config=cfg)
        solver.setup(packing)
        result = solver.solve()
        self.assertEqual(result.shape, (4,))  # 2 D_eff + 2 K_eff


if __name__ == "__main__":
    unittest.main()
