import math
import os
import tempfile
import unittest

import numpy as np

from porescalemc.config import FourierConfig
from porescalemc.geometry.fourier_field import (
    autocorrelation,
    correlation_length,
    grain_fourier_transform,
    packing_fourier_transform,
    porosity_field,
    save_structured_vti_fields,
    smoothing_kernel,
    structure_factor,
    upscaled_permeability,
)
from porescalemc.geometry.grains import Grain, Packing


class FourierTests(unittest.TestCase):
    def test_packing_fourier_dc_equals_solid_volume(self):
        packing = Packing([Grain.sphere([0.0, 0.0, 0.0], 0.1)], box=np.ones(3))
        F_k, *_ = packing_fourier_transform(packing, 4, 4, 4)
        self.assertAlmostEqual(F_k[0, 0, 0].real, packing.solid_volume)

    def test_ellipsoid_transform_reduces_to_sphere_for_isotropic_matrix(self):
        kx = np.array([[[0.0]], [[1.0]]])
        ky = np.zeros_like(kx)
        kz = np.zeros_like(kx)
        sphere = Grain.sphere([0.0, 0.0, 0.0], 0.2)
        ellipsoid = Grain.ellipsoid([0.0, 0.0, 0.0], np.eye(3) * 0.2)
        np.testing.assert_allclose(
            grain_fourier_transform(sphere, kx, ky, kz),
            grain_fourier_transform(ellipsoid, kx, ky, kz),
            rtol=1e-12,
            atol=1e-12,
        )

    def test_porosity_field_mean_matches_solid_fraction_reasonably(self):
        packing = Packing([Grain.sphere([0.0, 0.0, 0.0], 0.08)], box=np.ones(3))
        field = porosity_field(packing, FourierConfig(resolution=12, smoothing="cell"))
        self.assertAlmostEqual(float(field.mean()), 1.0 - packing.porosity, delta=1e-3)

    def test_gaussian_smoothing_has_unit_dc_gain(self):
        kx = np.array([[[0.0]], [[1.0]]])
        ky = np.zeros_like(kx)
        kz = np.zeros_like(kx)
        kernel = smoothing_kernel(kx, ky, kz, "gaussian", 0.1, 0.1, 0.1, sigma=2.0)
        self.assertAlmostEqual(float(kernel[0, 0, 0]), 1.0)
        self.assertLess(float(kernel[1, 0, 0]), 1.0)

    def test_gaussian_porosity_field_conserves_mean(self):
        packing = Packing([Grain.sphere([0.0, 0.0, 0.0], 0.07)], box=np.ones(3))
        field = porosity_field(packing, FourierConfig(resolution=12, smoothing="gaussian"))
        self.assertAlmostEqual(float(field.mean()), 1.0 - packing.porosity, delta=1e-3)

    def test_upscaled_permeability_returns_positive_value(self):
        packing = Packing([Grain.sphere([0.0, 0.0, 0.0], 0.05)], box=np.ones(3))
        value = upscaled_permeability(
            packing,
            FourierConfig(resolution=8, smoothing="cell"),
            law="kozeny_carman",
            grain_size=0.1,
        )
        self.assertTrue(math.isfinite(value))
        self.assertGreater(value, 0.0)

    def test_parseval_theorem_holds(self):
        """Parseval: sum |f(x)|^2 dx == sum |F(k)|^2 dk (discrete form).

        For the solid indicator field f ∈ {0,1}, the L2 norm equals the solid
        volume.  The Fourier-side counterpart is (1/V) Σ_k |F_k|^2 = solid_volume.
        """
        packing = Packing([Grain.sphere([0.0, 0.0, 0.0], 0.1)], box=np.ones(3))
        cfg = FourierConfig(resolution=10, smoothing="none")
        F_k, *_ = packing_fourier_transform(packing, 10, 10, 10)
        V = packing.box[0] * packing.box[1] * packing.box[2]
        N = 10 * 10 * 10
        # Discrete Parseval: (1/N) Σ |F_k|^2 == (V/N) * Σ_x |f(x)|^2
        # i.e. (1/(N*V)) Σ |F_k|^2 == mean(f²) == mean(f) (since f ∈ {0,1})
        # which equals solid_fraction = solid_volume / V
        parseval_lhs = float(np.sum(np.abs(F_k) ** 2)) / (N * V)
        parseval_rhs = packing.solid_volume / V
        self.assertAlmostEqual(parseval_lhs, parseval_rhs, delta=0.01)

    def test_autocorrelation_peaks_at_zero_lag(self):
        """Autocorrelation must have its maximum at zero lag (C[0,0,0] = 1)."""
        packing = Packing([Grain.sphere([0.0, 0.0, 0.0], 0.08)], box=np.ones(3))
        cfg = FourierConfig(resolution=8, smoothing="cell")
        C = autocorrelation(packing, cfg)
        self.assertAlmostEqual(float(C[0, 0, 0]), 1.0, places=10)
        # All other values should be <= 1
        self.assertLessEqual(float(C.max()), 1.0 + 1e-10)

    def test_correlation_length_is_positive_and_finite(self):
        """Correlation length must be a positive, finite scalar."""
        packing = Packing(
            [Grain.sphere([0.5, 0.5, 0.5], 0.1)],
            box=np.array([1.0, 1.0, 1.0]),
        )
        cfg = FourierConfig(resolution=8, smoothing="cell")
        lc = correlation_length(packing, cfg)
        self.assertTrue(math.isfinite(lc))
        self.assertGreater(lc, 0.0)

    def test_save_structured_vti_fields_writes_xml_image_data(self):
        """VTI export should be real XML ImageData, not legacy VTK text."""
        fields = {
            "solid_fraction": np.ones((2, 3, 4)),
            "porosity": np.zeros((2, 3, 4)),
        }
        with tempfile.NamedTemporaryFile(suffix=".vti", delete=False) as tmp:
            path = tmp.name
        try:
            save_structured_vti_fields(
                fields,
                origin=(0.0, 0.0, 0.0),
                spacing=(0.5, 0.25, 0.125),
                filename=path,
                active_scalar="solid_fraction",
            )
            with open(path, "r", encoding="utf-8") as f:
                text = f.read()
            self.assertIn('<VTKFile type="ImageData"', text)
            self.assertIn('WholeExtent="0 1 0 2 0 3"', text)
            self.assertIn('Name="solid_fraction"', text)
            self.assertIn('Name="porosity"', text)
        finally:
            os.unlink(path)


if __name__ == "__main__":
    unittest.main()
