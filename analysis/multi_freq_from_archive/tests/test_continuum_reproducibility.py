from __future__ import annotations

import csv
import math
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PACKAGE_DIR.parents[1]
sys.path.insert(0, str(PACKAGE_DIR))

import s230_reference  # noqa: E402
from selected_clean_manifest import (  # noqa: E402
    load_selected_clean_manifest,
    resolve_manifest_entry_path,
    selected_image_path,
)


class SelectedCleanManifestTests(unittest.TestCase):
    def test_only_clean_runner_enables_selected_manifest(self):
        clean_runner = (PACKAGE_DIR / "run_continuum_fits_clean.sh").read_text()
        dirty_runner = (PACKAGE_DIR / "run_continuum_fits_dirty.sh").read_text()
        self.assertIn("--use_selected_clean_manifest", clean_runner)
        self.assertNotIn("--use_selected_clean_manifest", dirty_runner)

    def test_readme_generates_every_manifest_selected_iteration(self):
        readme = (REPOSITORY_ROOT / "README.md").read_text()
        for iteration in (6, 7, 10):
            self.assertIn(f"--sciniter {iteration}", readme)

    def test_released_manifest_resolves_selected_flag2_images_from_repo_root(self):
        manifest = REPOSITORY_ROOT / "analysis/sma_calibration/selected_clean_iterations.csv"
        selected, sources = load_selected_clean_manifest(
            manifest,
            date_code="230316",
            flag_number=2,
            repository_root=REPOSITORY_ROOT,
        )

        self.assertEqual(set(selected), {"NGC1068", "NGC3079", "NGC4258"})
        self.assertEqual(set(sources), {"NGC1068", "NGC3079", "NGC4258"})
        expected_iterations = {"NGC1068": 10, "NGC3079": 7, "NGC4258": 6}
        for source, iteration in expected_iterations.items():
            path = selected[source].replace("\\", "/")
            self.assertTrue(path.startswith(str(REPOSITORY_ROOT).replace("\\", "/") + "/"))
            self.assertIn(f"/images/clean/nitercal_{iteration}/", path)
            self.assertNotIn("/nitercal_0/", path)

    def test_manifest_local_paths_resolve_beside_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest_dir = root / "portable"
            image = manifest_dir / "images/clean/nitercal_6/NGC4258_test.image"
            image.parent.mkdir(parents=True)
            image.touch()
            manifest = manifest_dir / "selected.csv"
            manifest.write_text(
                "source,date_code,flagnum,image_path\n"
                "NGC4258,230316,2,images/clean/nitercal_6/NGC4258_test.image\n"
            )

            selected, _ = load_selected_clean_manifest(
                manifest,
                "230316",
                2,
                repository_root=root,
            )
            self.assertEqual(selected_image_path(selected, "NGC4258"), str(image.resolve()))

    def test_dot_analysis_entry_is_repository_relative_even_when_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "analysis/sma_calibration/selected.csv"
            manifest.parent.mkdir(parents=True)
            manifest.touch()
            resolved = resolve_manifest_entry_path(
                "./analysis/sma_calibration/missing.image",
                manifest,
                repository_root=root,
            )
            self.assertEqual(
                resolved,
                (root / "analysis/sma_calibration/missing.image").resolve(),
            )

    def test_missing_selected_image_fails_instead_of_falling_back(self):
        selected = {"NGC4258": "/missing/clean/nitercal_6/NGC4258_test.image"}
        with self.assertRaisesRegex(FileNotFoundError, "refusing to substitute"):
            selected_image_path(selected, "NGC4258", exists=lambda _: False)

    def test_missing_manifest_fails_instead_of_enabling_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaisesRegex(FileNotFoundError, "selection record"):
                load_selected_clean_manifest("missing.csv", repository_root=root)


class S230ReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source_table = PACKAGE_DIR / "fitsumfiles/final/fitsummary_final_withextras.txt"
        cls.rows = s230_reference.build_s230_reference_rows(cls.source_table)
        cls.by_source = {str(row["Source"]): row for row in cls.rows}

    def test_corrected_scientific_central_values_are_reproduced(self):
        expected = {
            "Circinus": 15.499750827004508,
            "ESO558-G009": 0.8190724285773601,
            "IC2560": 1.9742190595784141,
            "MRK1029": 2.8112008668238873,
            "NGC1068": 10.723826310128143,
            "NGC1194": 1.5159610255163574,
            "NGC1320": 0.7804450642546338,
            "NGC2273": 8.176888377947352,
            "NGC2960": 0.2250387857150855,
            "NGC3079": 14.904557617862627,
            "NGC3393": 0.3093864788422666,
            "NGC4258": 12.239913813400241,
            "NGC4388": 3.3817968527055853,
            "NGC4945": 17.527674085288034,
            "NGC5765B": 0.39481163015431925,
            "NGC6264": 2.5839309669058155,
            "NGC6323": 6.154095127287134,
            "UGC3789": 6.3324033035490945,
            "UGC6093": 5.957403241893575,
        }
        self.assertEqual(set(self.by_source), set(expected))
        for source, value in expected.items():
            self.assertTrue(
                math.isclose(
                    self.by_source[source]["S_230_mJy"],
                    value,
                    rel_tol=0.0,
                    abs_tol=1e-12,
                )
            )

    def test_selected_metadata_status_and_stable_ids_are_present(self):
        self.assertEqual(len(self.rows), 19)
        ids = set()
        for row in self.rows:
            self.assertIn(row["flux_status"], {"detection", "upper_limit"})
            for key in (
                "nu0_telescope",
                "nu0_date",
                "nu0_beam_arcsec",
                "nu0_reference",
                "nu0_observation_id",
            ):
                self.assertNotIn(row[key], (None, ""))
            ids.add(row["nu0_observation_id"])
            self.assertNotIn("/home/", str(row["nu0_reference"]))
        self.assertEqual(len(ids), len(self.rows))

        detections = [row for row in self.rows if row["flux_status"] == "detection"]
        upper_limits = [row for row in self.rows if row["flux_status"] == "upper_limit"]
        self.assertEqual(len(detections), 12)
        self.assertEqual(len(upper_limits), 7)
        self.assertEqual(
            Counter(str(row["nu0_telescope"]) for row in self.rows),
            {"ALMA": 11, "SMA": 7, "NOEMA": 1},
        )

        eso558 = self.by_source["ESO558-G009"]
        self.assertEqual(eso558["flux_status"], "detection")
        self.assertEqual(eso558["nu0_telescope"], "ALMA")
        self.assertEqual(eso558["nu0_date"], "2014/12/13_1")
        self.assertEqual(eso558["nu0_GHz"], 313.11)
        self.assertEqual(eso558["S_nu0_mJy"], 0.66)
        self.assertEqual(eso558["sigma_S_nu0_mJy"], 0.07352941176470588)
        self.assertEqual(eso558["nu0_beam_arcsec"], 0.71817)
        self.assertEqual(eso558["S_230_mJy"], 0.8190724285773601)
        self.assertEqual(eso558["sigma_S_230_mJy"], 0.09125138464542781)
        self.assertEqual(
            eso558["nu0_reference"],
            "ALMA/ESO558-G009/"
            "member.uid___A001_X122_X600.calibrated_final_cont_image-ESO558-G009.pbcor.fits",
        )
        self.assertEqual(
            eso558["nu0_observation_id"],
            "ESO558-G009|ALMA|2014-12-13|313.11GHz|"
            "member.uid-A001-X122-X600.calibrated-final-cont-image-ESO558-G009.pbcor.fits",
        )

        changed_references = {
            "NGC1068": ("ALMA", 348.32, 8.02, "2023/07/28_1"),
            "NGC1194": ("ALMA", 310.04, 1.23, "2025/08/09_1"),
            "NGC2960": ("ALMA", 316.43, 0.18, "2025/07/14_1"),
            "NGC3079": ("NOEMA", 229.66, 14.92, "2018/02/07_1"),
        }
        for source, (telescope, frequency, flux, date) in changed_references.items():
            with self.subTest(source=source):
                row = self.by_source[source]
                self.assertEqual(row["nu0_telescope"], telescope)
                self.assertEqual(row["nu0_GHz"], frequency)
                self.assertEqual(row["S_nu0_mJy"], flux)
                self.assertEqual(row["nu0_date"], date)

        self.assertEqual(self.by_source["NGC2960"]["flux_status"], "upper_limit")
        self.assertEqual(self.by_source["NGC3079"]["sigma_S_230_mJy"], None)

        ngc4258 = self.by_source["NGC4258"]
        self.assertEqual(ngc4258["flux_status"], "detection")
        self.assertEqual(ngc4258["nu0_telescope"], "SMA")
        self.assertEqual(ngc4258["nu0_date"], "2023/03/16_1")
        self.assertNotIn("_1", ngc4258["nu0_observation_id"])
        self.assertIn("230316-128-flag2-cal1", ngc4258["nu0_observation_id"])
        self.assertTrue(str(ngc4258["nu0_reference"]).startswith("analysis/sma_calibration/"))

    def test_detection_preferred_and_upper_limit_is_explicit_fallback(self):
        upper = {
            "Source": "TEST",
            "freq_ghz": 230.0,
            "beam_arcsec": 0.1,
            "is_upper_limit": True,
        }
        detection = {
            "Source": "TEST",
            "freq_ghz": 230.0,
            "beam_arcsec": 1.0,
            "is_upper_limit": False,
        }
        chosen = s230_reference.select_chosen_frequency_records([upper, detection])
        self.assertEqual(chosen["TEST"]["flux_status"], "detection")
        self.assertIs(chosen["TEST"]["is_upper_limit"], False)

        chosen = s230_reference.select_chosen_frequency_records([upper])
        self.assertEqual(chosen["TEST"]["flux_status"], "upper_limit")
        self.assertIs(chosen["TEST"]["is_upper_limit"], True)

        eligible_ancillary = {
            **detection,
            "Source": "NGC5495",
        }
        chosen = s230_reference.select_chosen_frequency_records(
            [detection, eligible_ancillary],
            allowed_sources={"TEST"},
        )
        self.assertEqual(set(chosen), {"TEST"})

    def test_formal_sigma_uses_rounded_flux_over_rounded_snr(self):
        row = self.by_source["NGC4258"]
        expected_sigma_nu0 = 12.41 / 6.725
        expected_sigma_230 = expected_sigma_nu0 * (230.0 / 225.51) ** -0.7
        self.assertTrue(
            math.isclose(row["sigma_S_nu0_mJy"], expected_sigma_nu0, rel_tol=0.0, abs_tol=1e-12)
        )
        self.assertTrue(
            math.isclose(row["sigma_S_230_mJy"], expected_sigma_230, rel_tol=0.0, abs_tol=1e-12)
        )

    def test_csv_retains_anchor_metadata_and_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "s230.csv"
            s230_reference.write_s230_reference_csv(self.rows, output)
            with output.open(newline="") as handle:
                written = list(csv.DictReader(handle))

        required = {
            "flux_status",
            "nu0_telescope",
            "nu0_date",
            "nu0_beam_arcsec",
            "nu0_reference",
            "nu0_observation_id",
        }
        self.assertTrue(required.issubset(written[0]))
        self.assertEqual({row["Source"] for row in written}, set(self.by_source))
        self.assertNotIn("NGC5495", self.by_source)
        self.assertFalse(any("/home/" in row["nu0_reference"] for row in written))


if __name__ == "__main__":
    unittest.main()
