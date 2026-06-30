"""Tests for OpenFOAM template structure and case-preparation utilities.

These tests are pure Python — no OpenFOAM binary is required.  They verify:
  1. All 7 templates in templates/openfoam/ have the required file structure.
  2. write_porosity_3d() produces a correctly formatted file.
  3. write_domainsize() produces a valid #include fragment.
  4. prepare_case() copies a template and injects the generated files.
"""

from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
OPENFOAM_TEMPLATES = REPO_ROOT / "templates" / "openfoam"

EXPECTED_TEMPLATES = [
    "darcy",
    "laplacian",
    "NavierStokes",
    "snappy",
    "comp2phase",
    "incomp2phase",
    "discDarcy",
]


def _small_packing():
    """Tiny periodic packing for fast tests."""
    from porescalemc.geometry.placement import sample_packing
    from porescalemc.config import PackingConfig
    return sample_packing(PackingConfig(mu=0.08, n_grains=5))


class TestTemplateStructure(unittest.TestCase):
    """Structural integrity checks — no solver calls needed."""

    def test_all_expected_templates_present(self):
        """All 7 expected template directories exist."""
        for name in EXPECTED_TEMPLATES:
            path = OPENFOAM_TEMPLATES / name
            self.assertTrue(path.is_dir(), f"Missing template: {name}")

    def test_templates_have_required_files(self):
        """Each template has system/controlDict, fvSchemes, fvSolution, constant/."""
        required = [
            "system/controlDict",
            "system/fvSchemes",
            "system/fvSolution",
            "constant",
        ]
        for name in EXPECTED_TEMPLATES:
            tmpl = OPENFOAM_TEMPLATES / name
            if not tmpl.is_dir():
                continue
            for rel in required:
                self.assertTrue(
                    (tmpl / rel).exists(),
                    f"Template '{name}' missing: {rel}",
                )

    def test_templates_have_openfoam13_version(self):
        """controlDict should reference OpenFOAM version 13."""
        for name in EXPECTED_TEMPLATES:
            ctrl = OPENFOAM_TEMPLATES / name / "system" / "controlDict"
            if not ctrl.exists():
                continue
            content = ctrl.read_text()
            self.assertIn(
                "13",
                content,
                f"Template '{name}' controlDict does not mention version 13",
            )

    def test_templates_have_allrun(self):
        """Each template must have an Allrun script."""
        for name in EXPECTED_TEMPLATES:
            tmpl = OPENFOAM_TEMPLATES / name
            if not tmpl.is_dir():
                continue
            self.assertTrue(
                (tmpl / "Allrun").exists(),
                f"Template '{name}' missing Allrun",
            )

    def test_darcy_has_porosity_properties(self):
        """Darcy template must reference porosity.3d in porosityProperties."""
        prop = OPENFOAM_TEMPLATES / "darcy" / "constant" / "porosityProperties"
        if not prop.exists():
            self.skipTest("darcy/constant/porosityProperties not found")
        content = prop.read_text()
        self.assertIn("porosity.3d", content)

    def test_validate_template_helper(self):
        """validate_template() returns empty list for structurally complete templates."""
        from porescalemc.openfoam import validate_template
        for name in EXPECTED_TEMPLATES:
            tmpl = OPENFOAM_TEMPLATES / name
            if not tmpl.is_dir():
                continue
            missing = validate_template(tmpl)
            self.assertEqual(
                missing, [],
                f"Template '{name}' structurally incomplete: {missing}",
            )


class TestWritePorosity3d(unittest.TestCase):

    def test_output_format(self):
        """porosity.3d first line is 'nx ny nz'; subsequent lines have 6 floats."""
        from porescalemc.openfoam import write_porosity_3d
        from porescalemc.config import FourierConfig
        packing = _small_packing()
        cfg = FourierConfig(resolution=5.0)  # coarse grid for speed
        with tempfile.NamedTemporaryFile(suffix=".3d", delete=False, mode="r") as tmp:
            tmppath = tmp.name
        try:
            write_porosity_3d(packing, tmppath, fourier_config=cfg)
            with open(tmppath) as fh:
                lines = fh.readlines()
            header = lines[0].split()
            self.assertEqual(len(header), 3, "Header should be 'nx ny nz'")
            nx, ny, nz = int(header[0]), int(header[1]), int(header[2])
            self.assertEqual(len(lines) - 1, nx * ny * nz)
            # Spot-check a data line
            row = lines[1].split()
            self.assertEqual(len(row), 6, "Each data row must have 6 fields")
            # Porosity values in [0, 1]
            val = float(row[3])
            self.assertGreaterEqual(val, 0.0)
            self.assertLessEqual(val, 1.0)
        finally:
            os.unlink(tmppath)

    def test_porosity_values_sum_to_mean_porosity(self):
        """Mean of written porosity values should match packing's volume fraction."""
        from porescalemc.openfoam import write_porosity_3d
        from porescalemc.config import FourierConfig
        packing = _small_packing()
        cfg = FourierConfig(resolution=8.0)
        with tempfile.NamedTemporaryFile(suffix=".3d", delete=False) as tmp:
            tmppath = tmp.name
        try:
            write_porosity_3d(packing, tmppath, fourier_config=cfg)
            with open(tmppath) as fh:
                lines = fh.readlines()[1:]  # skip header
            vals = [float(ln.split()[3]) for ln in lines]
            mean_porosity = sum(vals) / len(vals)
            # Should be roughly equal to phi_f — allow 5% tolerance
            expected = packing.porosity
            self.assertAlmostEqual(mean_porosity, expected, delta=0.05)
        finally:
            os.unlink(tmppath)


class TestWriteDomainsize(unittest.TestCase):

    def test_output_contains_grid_tokens(self):
        """domainsize must contain xgrid, ygrid, zgrid, scalegrid."""
        from porescalemc.openfoam import write_domainsize
        packing = _small_packing()
        with tempfile.NamedTemporaryFile(delete=False, mode="r", suffix=".txt") as tmp:
            tmppath = tmp.name
        try:
            write_domainsize(packing, tmppath, grid_res=20.0)
            content = open(tmppath).read()
            for token in ("xgrid", "ygrid", "zgrid", "scalegrid", "x1", "x2"):
                self.assertIn(token, content, f"Missing token '{token}' in domainsize")
        finally:
            os.unlink(tmppath)

    def test_grid_counts_match_res(self):
        """xgrid = round(grid_res * Lx)."""
        from porescalemc.openfoam import write_domainsize
        packing = _small_packing()
        with tempfile.NamedTemporaryFile(delete=False, mode="r", suffix=".txt") as tmp:
            tmppath = tmp.name
        try:
            write_domainsize(packing, tmppath, grid_res=30.0)
            content = open(tmppath).read()
            Lx = packing.box[0]
            expected_xgrid = max(1, int(30.0 * Lx))
            self.assertIn(f"xgrid {expected_xgrid};", content)
        finally:
            os.unlink(tmppath)


class TestPrepareCase(unittest.TestCase):

    def test_prepare_darcy_case(self):
        """prepare_case(darcy) creates case dir with domainsize and porosity.3d."""
        from porescalemc.openfoam import prepare_case
        from porescalemc.config import FourierConfig
        tmpl = OPENFOAM_TEMPLATES / "darcy"
        if not tmpl.is_dir():
            self.skipTest("darcy template not found")
        packing = _small_packing()
        cfg = FourierConfig(resolution=5.0)
        with tempfile.TemporaryDirectory() as tmpdir:
            case_dir = Path(tmpdir) / "test_darcy"
            result = prepare_case(
                tmpl, case_dir, packing,
                fourier_config=cfg, grid_res=10.0, overwrite=True,
            )
            self.assertTrue((result / "system" / "controlDict").exists())
            self.assertTrue((result / "constant" / "polyMesh" / "domainsize").exists())
            self.assertTrue((result / "constant" / "porosity.3d").exists())
            self.assertTrue((result / "system" / "tolerance").exists())

    def test_prepare_laplacian_case(self):
        """prepare_case(laplacian) creates a structurally valid case."""
        from porescalemc.openfoam import prepare_case
        from porescalemc.config import FourierConfig
        tmpl = OPENFOAM_TEMPLATES / "laplacian"
        if not tmpl.is_dir():
            self.skipTest("laplacian template not found")
        packing = _small_packing()
        cfg = FourierConfig(resolution=5.0)
        with tempfile.TemporaryDirectory() as tmpdir:
            case_dir = Path(tmpdir) / "test_laplacian"
            result = prepare_case(tmpl, case_dir, packing, fourier_config=cfg, grid_res=10.0)
            self.assertTrue(result.is_dir())
            self.assertTrue((result / "constant" / "porosity.3d").exists())

    def test_prepare_navier_stokes_case(self):
        """prepare_case(NavierStokes) creates a case with domainsize."""
        from porescalemc.openfoam import prepare_case
        from porescalemc.config import FourierConfig
        tmpl = OPENFOAM_TEMPLATES / "NavierStokes"
        if not tmpl.is_dir():
            self.skipTest("NavierStokes template not found")
        packing = _small_packing()
        cfg = FourierConfig(resolution=5.0)
        with tempfile.TemporaryDirectory() as tmpdir:
            case_dir = Path(tmpdir) / "test_ns"
            result = prepare_case(tmpl, case_dir, packing, fourier_config=cfg, grid_res=10.0)
            self.assertTrue((result / "constant" / "polyMesh" / "domainsize").exists())

    def test_prepare_all_templates(self):
        """prepare_case() succeeds for every available template."""
        from porescalemc.openfoam import prepare_case
        from porescalemc.config import FourierConfig
        packing = _small_packing()
        cfg = FourierConfig(resolution=4.0)
        with tempfile.TemporaryDirectory() as tmpdir:
            for name in EXPECTED_TEMPLATES:
                tmpl = OPENFOAM_TEMPLATES / name
                if not tmpl.is_dir():
                    continue
                case_dir = Path(tmpdir) / name
                result = prepare_case(tmpl, case_dir, packing, fourier_config=cfg, grid_res=8.0)
                self.assertTrue(
                    (result / "system" / "controlDict").exists(),
                    f"prepare_case({name}) did not copy controlDict",
                )
                self.assertTrue(
                    (result / "constant" / "porosity.3d").exists(),
                    f"prepare_case({name}) did not write porosity.3d",
                )


if __name__ == "__main__":
    unittest.main()
