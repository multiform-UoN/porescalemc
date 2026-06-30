import unittest

import numpy as np

from porescalemc.geometry.grains import Grain, Packing
from porescalemc.solvers.dummy import DummySolver
from porescalemc.solvers.packing import PACKING_QOI_NAMES, PackingStatsSolver, packing_statistics


class SolverTests(unittest.TestCase):
    def test_packing_statistics_contains_expected_descriptors(self):
        packing = Packing(
            [Grain.sphere([0.0, 0.0, 0.0], 0.1), Grain.sphere([0.4, 0.0, 0.0], 0.1)],
            box=np.ones(3),
        )
        stats = packing_statistics(packing)
        self.assertEqual(set(PACKING_QOI_NAMES), set(stats))
        self.assertGreater(stats["specific_surface_area"], 0.0)
        self.assertGreater(stats["number_density"], 0.0)

    def test_dummy_solver_returns_numeric_vector(self):
        packing = Packing([Grain.sphere([0.0, 0.0, 0.0], 0.1)], box=np.ones(3))
        solver = DummySolver(packing)
        solver.setup(packing)
        result = solver.solve()
        self.assertEqual(result.shape, (len(PACKING_QOI_NAMES),))
        self.assertFalse(np.isnan(result[:6]).any())

    def test_packing_solver_can_select_qois(self):
        packing = Packing([Grain.sphere([0.0, 0.0, 0.0], 0.1)], box=np.ones(3))
        solver = PackingStatsSolver(packing, qoi_names=("porosity", "number_density"))
        solver.setup(packing)
        result = solver.solve()
        self.assertEqual(result.shape, (2,))
        self.assertAlmostEqual(result[0], packing.porosity)

    def test_spectral_solvers(self):
        from porescalemc.solvers.spectral import SpectralDiffusionSolver, SpectralStokesSolver
        from porescalemc.config import SpectralConfig

        packing = Packing([Grain.sphere([0.0, 0.0, 0.0], 0.08)], box=np.ones(3))
        cfg = SpectralConfig(resolution=8, eta=1e-3, max_iter=50, tol=1e-4, n_directions=1)

        diff_solver = SpectralDiffusionSolver(spectral_config=cfg)
        diff_solver.setup(packing)
        qoi_diff = diff_solver.solve()
        self.assertEqual(qoi_diff.shape, (1,))
        self.assertTrue(np.isfinite(qoi_diff[0]))
        self.assertGreater(qoi_diff[0], 0.0)

        stokes_solver = SpectralStokesSolver(spectral_config=cfg)
        stokes_solver.setup(packing)
        qoi_stokes = stokes_solver.solve()
        self.assertEqual(qoi_stokes.shape, (1,))
        self.assertTrue(np.isfinite(qoi_stokes[0]))
        self.assertGreater(qoi_stokes[0], 0.0)

    def test_solver_registry_exposes_qois(self):
        from porescalemc.solvers import available_qois, available_solvers, solver_class
        from porescalemc.solvers.voxel import VoxelDiffusionSolver

        names = available_solvers()
        self.assertIn("packing", names)
        self.assertIn("voxel_diffusion", names)
        self.assertIn("tet_diffusion", names)
        self.assertIn("porosity", available_qois("packing"))
        self.assertIn("diffusivity_x", available_qois("voxel_diffusion"))
        self.assertIn("diffusivity_x", available_qois("tet_diffusion"))
        self.assertIs(solver_class("voxel_diffusion"), VoxelDiffusionSolver)

    def test_voxel_diffusion_empty_domain_is_free_diffusion(self):
        from porescalemc.config import FourierConfig
        from porescalemc.solvers.voxel import VoxelDiffusionSolver

        packing = Packing([], box=np.ones(3))
        solver = VoxelDiffusionSolver(
            fourier_config=FourierConfig(resolution=5, smoothing="cell"),
            n_directions=2,
            max_iter=80,
            tol=1e-8,
        )
        solver.setup(packing)
        result = solver.solve()
        self.assertEqual(result.shape, (2,))
        np.testing.assert_allclose(result, np.ones(2), atol=1e-10)
        self.assertIn("voxel_diffusivity", solver.solution_fields())

    def test_voxel_diffusion_porous_domain_is_finite(self):
        from porescalemc.config import FourierConfig
        from porescalemc.solvers.voxel import VoxelDiffusionSolver

        packing = Packing([Grain.sphere([0.5, 0.5, 0.5], 0.12)], box=np.ones(3))
        solver = VoxelDiffusionSolver(
            fourier_config=FourierConfig(resolution=6, smoothing="cell"),
            n_directions=1,
            max_iter=100,
            tol=1e-6,
        )
        solver.setup(packing)
        result = solver.solve()
        self.assertEqual(result.shape, (1,))
        self.assertTrue(np.isfinite(result[0]))
        self.assertGreater(result[0], 0.0)
        self.assertLessEqual(result[0], 1.05)

    def test_tet_diffusion_soft_dependency_guard(self):
        import porescalemc.solvers.tet as tet_module

        old_has_gmsh = tet_module._HAS_GMSH
        try:
            tet_module._HAS_GMSH = False
            with self.assertRaises(ImportError):
                tet_module.TetDiffusionSolver()
        finally:
            tet_module._HAS_GMSH = old_has_gmsh


if __name__ == "__main__":
    unittest.main()
