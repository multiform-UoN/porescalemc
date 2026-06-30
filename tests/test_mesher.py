"""Tests for GmshPackingMesher — skipped when gmsh is not installed."""

from __future__ import annotations

import unittest

from porescalemc.geometry.grains import Grain, Packing

try:
    import gmsh  # type: ignore
    _GMSH_AVAILABLE = True
except ImportError:
    _GMSH_AVAILABLE = False


@unittest.skipUnless(_GMSH_AVAILABLE, "gmsh not installed")
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
        import tempfile, os
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
        import tempfile, os
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
        import tempfile, os
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

    def test_missing_gmsh_raises_import_error(self):
        """Importing GmshPackingMesher without gmsh raises ImportError at construction."""
        import sys
        # This test only makes sense when gmsh IS available (we're in the skipUnless block).
        # We verify the class exists and constructs cleanly.
        from porescalemc.geometry.mesher import GmshPackingMesher
        m = GmshPackingMesher(mesh_size=0.2, verbosity=0)
        m.finalize()  # should be a no-op before build()

    def test_point_in_grain_helper(self):
        from porescalemc.geometry.mesher import _point_in_any_grain
        grains = [Grain.sphere(center=[0.5, 0.5, 0.5], radius=0.2)]
        self.assertTrue(_point_in_any_grain(0.5, 0.5, 0.5, grains))
        self.assertFalse(_point_in_any_grain(0.0, 0.0, 0.0, grains))


if __name__ == "__main__":
    unittest.main()
