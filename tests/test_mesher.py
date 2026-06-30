"""Tests for GmshPackingMesher — skipped when gmsh is not installed."""

from __future__ import annotations

import unittest
import os

from porescalemc.geometry.grains import Grain, Packing

try:
    import gmsh  # type: ignore
    _GMSH_AVAILABLE = True
except ImportError:
    _GMSH_AVAILABLE = False

_RUN_GMSH_TESTS = _GMSH_AVAILABLE and os.environ.get("PORESCALEMC_RUN_GMSH_TESTS") == "1"


class TestGmshMesherPurePython(unittest.TestCase):
    """Pure-Python checks for mesher helper logic and soft dependency behavior."""

    def test_point_in_grain_helper_supports_spheres_and_rotated_ellipsoids(self):
        from porescalemc.geometry.mesher import _point_in_any_grain
        import numpy as np

        sphere = Grain.sphere(center=[0.5, 0.5, 0.5], radius=0.2)
        self.assertTrue(_point_in_any_grain(0.5, 0.5, 0.5, [sphere]))
        self.assertFalse(_point_in_any_grain(0.0, 0.0, 0.0, [sphere]))

        angle = np.pi / 4
        R = np.array([
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle),  np.cos(angle), 0.0],
            [0.0,            0.0,           1.0],
        ])
        M = R @ np.diag([0.25, 0.12, 0.08])
        ellipsoid = Grain.ellipsoid(center=[0.5, 0.5, 0.5], M=M)
        self.assertTrue(_point_in_any_grain(0.5, 0.5, 0.5, [ellipsoid]))
        self.assertFalse(_point_in_any_grain(0.9, 0.9, 0.9, [ellipsoid]))

    def test_soft_dependency_constructor_guard(self):
        import porescalemc.geometry.mesher as mesher_module
        import porescalemc.solvers.gmsh_mesher as solver_module

        old_gmsh_available = mesher_module._GMSH_AVAILABLE
        old_has_mesher = solver_module._HAS_MESHER
        try:
            mesher_module._GMSH_AVAILABLE = False
            solver_module._HAS_MESHER = False
            with self.assertRaises(ImportError):
                mesher_module.GmshPackingMesher()
            with self.assertRaises(ImportError):
                solver_module.GmshMesherQoISolver()
        finally:
            mesher_module._GMSH_AVAILABLE = old_gmsh_available
            solver_module._HAS_MESHER = old_has_mesher

    def test_periodic_surface_pairing_matches_translated_boxes(self):
        import porescalemc.geometry.mesher as mesher_module

        class FakeModel:
            def __init__(self):
                self.boxes = {
                    1: (0.0, 0.0, 0.0, 0.0, 0.4, 0.4),
                    2: (0.0, 0.5, 0.5, 0.0, 1.0, 1.0),
                    11: (1.0, 0.5, 0.5, 1.0, 1.0, 1.0),
                    12: (1.0, 0.0, 0.0, 1.0, 0.4, 0.4),
                }

            def getBoundingBox(self, dim, tag):
                del dim
                return self.boxes[tag]

        class FakeGmsh:
            model = FakeModel()

        old_gmsh = getattr(mesher_module, "gmsh", None)
        try:
            mesher_module.gmsh = FakeGmsh()
            pairs = mesher_module._matching_periodic_surface_pairs(
                slave_tags=[1, 2],
                master_tags=[11, 12],
                translation=[1.0, 0.0, 0.0],
            )
            self.assertEqual(set(pairs), {(1, 12), (2, 11)})
        finally:
            if old_gmsh is None:
                delattr(mesher_module, "gmsh")
            else:
                mesher_module.gmsh = old_gmsh


@unittest.skipUnless(
    _RUN_GMSH_TESTS,
    "gmsh mesh-generation tests require gmsh and PORESCALEMC_RUN_GMSH_TESTS=1",
)
class TestGmshPackingMesher(unittest.TestCase):
    """Smoke tests: build a tiny mesh and verify it is non-empty."""

    def _single_sphere_packing(self, radius: float = 0.2) -> Packing:
        return Packing(
            box=(1.0, 1.0, 1.0),
            grains=[Grain.sphere(center=[0.5, 0.5, 0.5], radius=radius)],
        )

    def test_import_succeeds(self):
        from porescalemc.geometry.mesher import GmshPackingMesher
        mesher = GmshPackingMesher(mesh_size=0.15, periodic=False, verbosity=0)
        self.assertIsNotNone(mesher)

    def test_build_non_periodic_runs(self):
        from porescalemc.geometry.mesher import GmshPackingMesher
        import tempfile
        packing = self._single_sphere_packing()
        mesher = GmshPackingMesher(mesh_size=0.15, periodic=False, verbosity=0)
        mesher.build(packing)
        with tempfile.NamedTemporaryFile(suffix=".msh", delete=False) as f:
            path = f.name
        try:
            mesher.write(path)
            self.assertGreater(os.path.getsize(path), 0)
        finally:
            os.unlink(path)
            mesher.finalize()

    def test_build_periodic_runs(self):
        from porescalemc.geometry.mesher import GmshPackingMesher
        import tempfile
        packing = self._single_sphere_packing(radius=0.15)
        mesher = GmshPackingMesher(mesh_size=0.18, periodic=True, verbosity=0)
        mesher.build(packing)
        with tempfile.NamedTemporaryFile(suffix=".msh", delete=False) as f:
            path = f.name
        try:
            mesher.write(path)
            self.assertGreater(os.path.getsize(path), 0)
        finally:
            os.unlink(path)
            mesher.finalize()

    def test_ellipsoid_grain(self):
        from porescalemc.geometry.mesher import GmshPackingMesher
        import tempfile
        import numpy as np
        # Axis-aligned ellipsoid: diagonal M with semi-axes 0.3, 0.2, 0.15
        M = np.diag([0.3, 0.2, 0.15])
        packing = Packing(
            box=(1.0, 1.0, 1.0),
            grains=[Grain.ellipsoid(center=[0.5, 0.5, 0.5], M=M)],
        )
        mesher = GmshPackingMesher(mesh_size=0.18, periodic=False, verbosity=0)
        mesher.build(packing)
        with tempfile.NamedTemporaryFile(suffix=".msh", delete=False) as f:
            path = f.name
        try:
            mesher.write(path)
            self.assertGreater(os.path.getsize(path), 0)
        finally:
            os.unlink(path)
            mesher.finalize()

    def test_constructor_available_when_gmsh_enabled(self):
        from porescalemc.geometry.mesher import GmshPackingMesher
        m = GmshPackingMesher(mesh_size=0.2, verbosity=0)
        m.finalize()  # should be a no-op before build()

    def test_mlmc_integration_with_gmsh_mesher(self):
        """Test that Gmsh mesher integrates with MLMC machinery via custom solver.

        Uses level-dependent mesh_size via factory. QoI is mesh-based porosity
        + scaled element count from GmshMesherQoISolver.
        Uses a rotated ellipsoid to verify full orientation support.
        """
        from porescalemc.geometry.grains import Grain, Packing
        from porescalemc.solvers.gmsh_mesher import GmshMesherQoISolver
        from porescalemc.mlmc.estimator import MLMCEstimator
        from porescalemc.config import MLMCConfig, PackingConfig
        import numpy as np

        # Rotated ellipsoid to test the fix
        angle = np.pi / 6
        R = np.array([
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle), np.cos(angle), 0],
            [0, 0, 1]
        ])
        S = np.diag([0.18, 0.12, 0.09])
        M = R @ S
        packing0 = Packing(
            box=(1.0, 1.0, 1.0),
            grains=[Grain.ellipsoid(center=[0.5, 0.5, 0.5], M=M)]
        )

        def geo_factory(level, name, pcfg):
            # For demo, same simple packing, level affects solver not geo here
            return packing0

        def solver_factory(packing, level=0):
            # Level dependent mesh size: coarser for low levels
            mesh_size = 0.25 / (2 ** level)   # 0.25, 0.125, 0.0625
            return GmshMesherQoISolver(mesh_size=mesh_size, periodic=False, verbosity=0)

        mlmc_cfg = MLMCConfig(
            n_levels=2,
            m0=2,
            min_samples=1,
            tolerance=0.2,
            n_workers=1,
            algorithm="fixed",
        )

        est = MLMCEstimator(
            name="gmsh_mlmc_test",
            packing_config=PackingConfig(),
            mlmc_config=mlmc_cfg,
            geometry_factory=geo_factory,
            solver_factory=solver_factory,
        )
        mean, err = est.run()

        # Basic sanity: should return 2-element vector, finite
        self.assertEqual(mean.shape, (2,))
        self.assertTrue(np.all(np.isfinite(mean)))
        self.assertTrue(np.all(np.isfinite(err)))

        self.assertGreater(mean[0], 0.0)
        self.assertLess(mean[0], 1.0)


if __name__ == "__main__":
    unittest.main()
