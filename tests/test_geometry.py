import math
import unittest

import numpy as np

from porescalemc.config import PackingConfig
from porescalemc.geometry.grains import (
    Grain,
    Packing,
    ellipsoid_distance,
    sample_random_orientation,
)
from porescalemc.geometry.jodrey_tory import jodrey_tory
from porescalemc.geometry.placement import sample_packing, seed_from_name


class GeometryTests(unittest.TestCase):
    def test_sphere_distance_touching(self):
        g1 = Grain.sphere([0.0, 0.0, 0.0], 0.1)
        g2 = Grain.sphere([0.2, 0.0, 0.0], 0.1)
        self.assertAlmostEqual(ellipsoid_distance(g1, g2), 1.0)

    def test_periodic_distance_uses_minimum_image(self):
        g1 = Grain.sphere([-0.49, 0.0, 0.0], 0.02)
        g2 = Grain.sphere([0.49, 0.0, 0.0], 0.02)
        self.assertAlmostEqual(ellipsoid_distance(g1, g2, box=np.ones(3)), 0.5)

    def test_random_orientation_preserves_radii(self):
        M = sample_random_orientation([0.1, 0.2, 0.3], np.random.default_rng(1)).reshape(3, 3)
        radii = np.sqrt((M**2).sum(axis=0))
        np.testing.assert_allclose(np.sort(radii), [0.1, 0.2, 0.3])
        self.assertGreater(np.linalg.det(M), 0.0)

    def test_sample_packing_is_seed_reproducible(self):
        cfg = PackingConfig(n_grains=5, psd="constant", periodic=False, detached=True)
        p1 = sample_packing(cfg, np.random.default_rng(seed_from_name("same")))
        p2 = sample_packing(cfg, np.random.default_rng(seed_from_name("same")))
        np.testing.assert_allclose(p1.centers, p2.centers)
        self.assertEqual(p1.n_grains, p2.n_grains)

    def test_jodrey_tory_reduces_simple_overlap_without_mutating_input(self):
        g1 = Grain.sphere([-0.05, 0.0, 0.0], 0.1)
        g2 = Grain.sphere([0.05, 0.0, 0.0], 0.1)
        packing = Packing([g1, g2], box=np.ones(3))
        cfg = PackingConfig(
            periodic=False,
            detached_bc=0.0,
            max_tries=50,
            use_jodrey_tory=True,
            jt_min_dist=0.99,
            jt_eps=0.75,
            jt_n_moves=1,
        )
        compacted = jodrey_tory(packing, cfg, np.random.default_rng(2))
        self.assertLess(ellipsoid_distance(g1, g2), 0.99)
        self.assertGreaterEqual(
            ellipsoid_distance(compacted.grains[0], compacted.grains[1]),
            0.98,
        )
        np.testing.assert_allclose(packing.centers, np.array([[-0.05, 0.0, 0.0], [0.05, 0.0, 0.0]]))

    def test_packing_porosity_matches_solid_volume_estimate(self):
        grains = [Grain.sphere([0, 0, 0], 0.1), Grain.sphere([0.5, 0, 0], 0.1)]
        packing = Packing(grains, box=np.ones(3))
        expected = 1.0 - 2.0 * 4.0 / 3.0 * math.pi * 0.1**3
        self.assertAlmostEqual(packing.porosity, expected)


if __name__ == "__main__":
    unittest.main()
