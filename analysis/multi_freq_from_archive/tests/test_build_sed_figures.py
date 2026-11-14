from __future__ import annotations

import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import matplotlib.image as mpimg


PACKAGE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_DIR))

from build_sed_figures import (  # noqa: E402
    DEFAULT_BHS_PARAMS,
    DEFAULT_INPUT,
    EXPECTED_SOURCE_ROWS,
    TARGET_SOURCES,
    build_figures,
    load_angular_scales,
    load_continuum_points,
    shared_axes,
)
import build_sed_figures as sed_builder  # noqa: E402


EXPECTED_HASHES = {
    "NGC4258": "bbf6e4231c4945afca6b2b6ae538775a0ebc07442fb216364ed0eabae9098f1d",
    "NGC1194": "410c05bf8d0b61a487442827385dd618ac5c96cd332407093f04672d8d66141a",
    "NGC1068": "347906bf2cf19fe7c5d980d0cc6937efd52d8008f4ded3103d3f2a92493c4e8f",
    "NGC3079": "73efddf740ce211f3c0c2cf0e38e1e75203258cefef350da147948efc027a80b",
    "NGC4945": "9c45f9a5edcdd3dacf23e09fc0e6fd48364583ef23cd935e2565d77022daaef2",
    "Circinus": "3c797c93bdf996134174853c9964f50d1c6bc13648b8586d1cc4204eb5fdad6e",
}


class SedFigureDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.points = load_continuum_points(DEFAULT_INPUT)

    def test_all_manuscript_rows_are_present(self) -> None:
        self.assertEqual(
            EXPECTED_SOURCE_ROWS,
            {
                "NGC4258": 82,
                "NGC1194": 9,
                "NGC1068": 88,
                "NGC3079": 79,
                "NGC4945": 45,
                "Circinus": 53,
            },
        )
        counts = {source: 0 for source in TARGET_SOURCES}
        for point in self.points:
            if point.source in counts:
                counts[point.source] += 1
        self.assertEqual(counts, EXPECTED_SOURCE_ROWS)

    def test_ngc4945_uses_all_corrected_nuclear_refits(self) -> None:
        points = [point for point in self.points if point.source == "NGC4945"]
        self.assertEqual(len(points), 45)
        self.assertTrue(all(not point.upper_limit for point in points))

        by_frequency = {
            round(10 ** point.log_frequency_ghz, 2): 10 ** point.log_flux_density_mjy
            for point in points
        }
        self.assertTrue(math.isclose(by_frequency[85.91], 55.19, abs_tol=1e-10))
        self.assertTrue(math.isclose(by_frequency[349.40], 13.08, abs_tol=1e-10))
        self.assertTrue(math.isclose(by_frequency[351.37], 55.06, abs_tol=1e-10))

    def test_shared_notebook_axes_are_preserved(self) -> None:
        axes = shared_axes(self.points)
        self.assertEqual(axes.frequency_limits, (-0.15000000000000002, 3.15))
        self.assertEqual(axes.flux_limits, (-1.25, 4.25))
        self.assertEqual(axes.beam_color_limits, (-2.0, 2.0))
        self.assertEqual(axes.point_color_limits, (-2.2, 2.2))

    def test_angular_scale_insets_follow_canonical_distances(self) -> None:
        scales = load_angular_scales(DEFAULT_BHS_PARAMS)
        expected = {
            "NGC4258": ("0.03673", "0.00054"),
            "NGC1194": ("0.255", "0.019"),
            "NGC1068": ("0.0520", "0.0025"),
            "NGC3079": ("0.0839", "0.0073"),
            "NGC4945": ("0.0179", None),
            "Circinus": ("0.0204", "0.0039"),
        }
        self.assertEqual(
            {
                source: (scale.display_value, scale.display_uncertainty)
                for source, scale in scales.items()
            },
            expected,
        )
        self.assertEqual(scales["NGC1068"].distance_mpc, 10.72)
        self.assertEqual(scales["NGC1068"].distance_unc_mpc, 0.52)

    def test_angular_scale_rounding_is_decimal_half_up(self) -> None:
        self.assertEqual(
            sed_builder._format_angular_scale(0.05205, 0.0025),
            ("0.0521", "0.0025"),
        )
        self.assertEqual(
            sed_builder._format_angular_scale(0.01245, None),
            ("0.0125", None),
        )

    def test_builder_has_no_network_dependency(self) -> None:
        source = (PACKAGE_DIR / "build_sed_figures.py").read_text(encoding="utf-8")
        for forbidden in ("requests", "urllib", "astroquery", "urlopen"):
            self.assertNotIn(forbidden, source)

    def test_historical_notebooks_cannot_overwrite_publication_seds(self) -> None:
        diagnostic_directories = {
            "Beam_Date_Flux_Graphs.ipynb": "SpectraGraphs_notebook_diagnostic",
            "Variability_Search.ipynb": "VariabilityGraphs_notebook_diagnostic",
            "get_contams.ipynb": "ContaminationGraphs_notebook_diagnostic",
            "flux_size_graph.ipynb": "FluxSizeGraphs_notebook_diagnostic",
        }
        for filename, diagnostic_directory in diagnostic_directories.items():
            notebook = (PACKAGE_DIR / filename).read_text(encoding="utf-8")
            json.loads(notebook)
            self.assertIn(diagnostic_directory, notebook, filename)
            self.assertNotIn("path = f'SpectraGraphs'", notebook, filename)
            self.assertNotIn("pylab.savefig(f'SpectraGraphs/", notebook, filename)
            self.assertNotIn("os.makedirs(f'SpectraGraphs')", notebook, filename)
            self.assertNotIn("os.makedirs('SpectraGraphs')", notebook, filename)

        sed_notebook = (PACKAGE_DIR / "Beam_Date_Flux_Graphs.ipynb").read_text(
            encoding="utf-8"
        )
        self.assertIn("HISTORICAL EXPLORATORY NOTEBOOK", sed_notebook)
        for stale in ("0.068', '0.010", "0.03673', '0.00040"):
            self.assertNotIn(stale, sed_notebook)
        for current in ("0.0520', '0.0025", "0.03673', '0.00054"):
            self.assertIn(current, sed_notebook)


class SedFigureArtifactTests(unittest.TestCase):
    def test_checked_in_pngs_are_exact_deterministic_artifacts(self) -> None:
        repository_dir = PACKAGE_DIR / "SpectraGraphs"
        with tempfile.TemporaryDirectory(prefix="ae8a42-sed-test-") as temp_dir:
            generated_dir = Path(temp_dir)
            generated_manifest = build_figures(DEFAULT_INPUT, generated_dir)

            expected_manifest = json.loads(
                (repository_dir / "sed_figure_manifest.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(generated_manifest, expected_manifest)
            for source in TARGET_SOURCES:
                filename = f"{source}_Spectra.png"
                self.assertEqual(
                    (generated_dir / filename).read_bytes(),
                    (repository_dir / filename).read_bytes(),
                )

    def test_release_hashes_and_dimensions_are_fixed(self) -> None:
        repository_dir = PACKAGE_DIR / "SpectraGraphs"
        manifest = json.loads(
            (repository_dir / "sed_figure_manifest.json").read_text(encoding="utf-8")
        )
        hashes = {
            source: manifest["artifacts"][source]["sha256"]
            for source in TARGET_SOURCES
        }
        self.assertEqual(hashes, EXPECTED_HASHES)
        for source in TARGET_SOURCES:
            image = mpimg.imread(repository_dir / f"{source}_Spectra.png")
            self.assertEqual(image.shape, (1361, 1729, 4))

    def test_revision_working_pngs_are_exact_repository_copies(self) -> None:
        repository_dir = PACKAGE_DIR / "SpectraGraphs"
        audit_root = PACKAGE_DIR.parents[2]
        revision_dir = audit_root / "revision_package/working_clean"
        if not revision_dir.is_dir():
            self.skipTest("adjacent manuscript revision package is not present")
        for source in TARGET_SOURCES:
            filename = f"{source}_Spectra.png"
            self.assertEqual(
                (revision_dir / filename).read_bytes(),
                (repository_dir / filename).read_bytes(),
            )


if __name__ == "__main__":
    unittest.main()
