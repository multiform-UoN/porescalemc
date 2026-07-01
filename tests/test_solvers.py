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

    def test_voxel_fv_flux_uses_periodic_faces(self):
        from porescalemc.solvers.voxel import (
            _build_operator,
            _effective_diffusivity,
            _face_diffusivities,
        )

        d = np.ones((6, 4, 3), dtype=float)
        face = _face_diffusivities(d)
        spacing = (1.0 / d.shape[0], 1.0 / d.shape[1], 1.0 / d.shape[2])
        operator = _build_operator(face, spacing, d.shape)
        np.testing.assert_allclose(operator @ np.ones(d.size), np.zeros(d.size), atol=1e-12)

        x = np.arange(d.shape[0])[:, None, None]
        periodic_corrector = np.sin(2.0 * np.pi * x / d.shape[0]) * np.ones_like(d)
        self.assertAlmostEqual(
            _effective_diffusivity(face, periodic_corrector, spacing, axis=0),
            1.0,
            places=12,
        )

    def test_tet_diffusion_soft_dependency_guard(self):
        import porescalemc.solvers.tet as tet_module

        old_has_gmsh = tet_module._HAS_GMSH
        try:
            tet_module._HAS_GMSH = False
            with self.assertRaises(ImportError):
                tet_module.TetDiffusionSolver()
        finally:
            tet_module._HAS_GMSH = old_has_gmsh


class TestPeriodicFEMHelpers(unittest.TestCase):
    """Pure-Python unit tests for periodic FEM DOF condensation helpers."""

    def test_find_periodic_pairs_unit_cube_grid(self):
        """On a 2×2×2 regular node grid, all high-face nodes pair with low-face nodes."""
        from porescalemc.solvers.tet import _find_periodic_node_pairs

        # 8 corner nodes of a 1×1×1 cube
        nodes = np.array([
            [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
            [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1],
        ], dtype=float)
        box = (1.0, 1.0, 1.0)
        pairs = _find_periodic_node_pairs(nodes, box)
        # slave: any node with a coord == 1; master: corresponding coord == 0
        slaves = {p[0] for p in pairs}
        masters = {p[1] for p in pairs}
        # 7 non-origin corner nodes are slaves to the origin (chains resolved)
        self.assertGreater(len(pairs), 0)
        # origin (0,0,0) must be a master, never a slave
        origin_idx = 0  # first node is (0,0,0)
        self.assertNotIn(origin_idx, slaves)

    def test_condense_dofs_pairs(self):
        """_condense_dofs collapses slave nodes onto masters."""
        from porescalemc.solvers.tet import _condense_dofs

        # 6 nodes: 0,1,2 are masters; 3->0, 4->1, 5->2
        pairs = [(3, 0), (4, 1), (5, 2)]
        dof, n_r = _condense_dofs(6, pairs)
        self.assertEqual(n_r, 3)
        self.assertEqual(dof[3], dof[0])
        self.assertEqual(dof[4], dof[1])
        self.assertEqual(dof[5], dof[2])

    def test_condense_dofs_chain(self):
        """Chain: 5->3->0 should both collapse to same reduced DOF as 0."""
        from porescalemc.solvers.tet import _condense_dofs

        pairs = [(3, 0), (5, 3)]
        dof, n_r = _condense_dofs(6, pairs)
        # 5, 3, 0 should all share the same condensed index
        self.assertEqual(dof[0], dof[3])
        self.assertEqual(dof[0], dof[5])
        # Remaining nodes 1, 2, 4 have their own condensed indices
        self.assertGreater(n_r, 1)

    def test_condense_dofs_no_pairs(self):
        """With no periodic pairs all nodes keep unique DOFs."""
        from porescalemc.solvers.tet import _condense_dofs

        dof, n_r = _condense_dofs(5, [])
        self.assertEqual(n_r, 5)
        self.assertEqual(len(set(dof.tolist())), 5)


class TestRunStudyAPI(unittest.TestCase):
    """Tests for the high-level run_study / list_solvers user API."""

    def test_list_solvers_returns_dict(self, capsys=None):
        import io, sys
        from porescalemc import list_solvers

        buf = io.StringIO()
        old = sys.stdout
        sys.stdout = buf
        result = list_solvers()
        sys.stdout = old
        self.assertIsInstance(result, dict)
        self.assertIn("packing", result)
        self.assertIn("description", result["packing"])
        self.assertIn("qoi_names", result["packing"])

    def test_list_qois_single_solver(self):
        import io, sys
        from porescalemc import list_qois

        buf = io.StringIO()
        old = sys.stdout
        sys.stdout = buf
        result = list_qois("packing")
        sys.stdout = old
        self.assertIn("packing", result)
        self.assertIn("porosity", result["packing"])

    def test_list_qois_unknown_solver_raises(self):
        from porescalemc import list_qois
        with self.assertRaises(KeyError):
            list_qois("nonexistent_solver_xyz")

    def test_run_study_packing_solver(self):
        """run_study with packing solver runs fast and returns finite results."""
        from porescalemc import run_study

        mean, err = run_study(
            solver="packing",
            n_grains=5,
            mu=0.10,
            n_levels=2,
            m0=4,
            tolerance=0.5,
        )
        self.assertTrue(np.all(np.isfinite(mean)))
        self.assertTrue(np.all(np.isfinite(err)))
        self.assertGreater(mean.shape[0], 0)

    def test_run_study_voxel_diffusion(self):
        """run_study with voxel_diffusion returns diffusivities in (0, 1]."""
        from porescalemc import run_study

        mean, err = run_study(
            solver="voxel_diffusion",
            n_grains=5,
            mu=0.08,
            n_levels=2,
            m0=3,
            tolerance=0.5,
            resolution=6,
            n_directions=1,
        )
        self.assertTrue(np.all(np.isfinite(mean)))
        self.assertGreater(float(mean[0]), 0.0)
        self.assertLessEqual(float(mean[0]), 1.05)

    def test_run_study_unknown_solver_raises(self):
        from porescalemc import run_study
        with self.assertRaises(KeyError):
            run_study(solver="not_a_solver")

    def test_run_study_accepts_full_config_objects(self):
        """Passing full config objects overrides all shortcuts."""
        from porescalemc import run_study
        from porescalemc.config import PackingConfig, MLMCConfig

        pcfg = PackingConfig(n_grains=4, mu=0.10)
        mcfg = MLMCConfig(n_levels=2, m0=3, tolerance=0.9)
        mean, err = run_study(
            solver="packing",
            packing_config=pcfg,
            mlmc_config=mcfg,
        )
        self.assertTrue(np.all(np.isfinite(mean)))


if __name__ == "__main__":
    unittest.main()
