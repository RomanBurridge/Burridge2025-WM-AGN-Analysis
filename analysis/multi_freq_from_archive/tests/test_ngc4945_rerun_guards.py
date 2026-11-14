from __future__ import annotations

import csv
import math
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_DIR))

from rerun_ngc4945_nuclear import (  # noqa: E402
    CORRECTION_FIELDS,
    CORRECTION_SCHEMA_VERSION,
    LEGACY_CORRECTION_FIELDS,
    LEGACY_PHOTOMETRY_FIELDS,
    PHOTOMETRY_FIELDS,
    merge_correction_manifest,
    preserve_pre_component_fields,
)
from chosen_frequency_anchors import (  # noqa: E402
    angular_separation_arcsec,
    sexagesimal_dec_to_deg,
    sexagesimal_ra_to_deg,
)

SCRIPT = PACKAGE_DIR / "rerun_ngc4945_nuclear.py"
PHOTOMETRY = PACKAGE_DIR / "ngc4945_nuclear_photometry.csv"
COMPONENTS = PACKAGE_DIR / "ngc4945_component_audit.csv"
CORRECTIONS = PACKAGE_DIR / "validated_continuum_corrections.csv"


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or ()), list(reader)


def write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class NGC4945RecordedRerunGuardTests(unittest.TestCase):
    def test_response_separation_is_derived_from_canonical_icrs_coordinates(
        self,
    ) -> None:
        _, rows = read_csv(COMPONENTS)
        source17 = next(
            row for row in rows if row["record_id"] == "rejected-source17-archived-anchor"
        )
        source18 = next(
            row for row in rows if row["record_id"] == "accepted-source18-fixed-canonical"
        )
        separation = angular_separation_arcsec(
            sexagesimal_ra_to_deg(source17["ra"]),
            sexagesimal_dec_to_deg(source17["dec"]),
            sexagesimal_ra_to_deg(source18["ra"]),
            sexagesimal_dec_to_deg(source18["dec"]),
        )

        _, photometry_rows = read_csv(PHOTOMETRY)
        anchor = next(row for row in photometry_rows if row["frequency_ghz"] == "349.40")
        geometric_beam = math.sqrt(
            float(anchor["beam_major_arcsec"]) * float(anchor["beam_minor_arcsec"])
        )
        self.assertAlmostEqual(separation, 0.2930145029022723, places=12)
        self.assertAlmostEqual(geometric_beam, 0.0783, places=4)
        self.assertAlmostEqual(
            separation / geometric_beam, 3.7422469339703524, places=12
        )

        response = (
            PACKAGE_DIR.parents[2]
            / "revision_package"
            / "response"
            / "CONSOLIDATED_CORRECTION_SUMMARY.tex"
        )
        if not response.exists():
            self.skipTest("adjacent correction response is not present")
        response_text = response.read_text(encoding="utf-8")
        self.assertIn(
            "separated by approximately 0.293 arcsec, or 3.74 times the "
            "0.0783-arcsec geometric beam",
            response_text,
        )
        self.assertNotIn("0.326 arcsec, or 4.16 times", response_text)

    def test_legacy_manifest_migration_is_exact_and_estimator_specific(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            manifest = Path(raw) / "corrections.csv"
            legacy = [
                {
                    "source": "OTHER",
                    "telescope": "ALMA",
                    "frequency_ghz": "100.00",
                    "date_label": "2020/01/01_1",
                    "product": "ALMA/OTHER/a.fits",
                    "old_flux_field": "1.00*",
                    "new_flux_field": "1.20*",
                    "old_snr": "NA",
                    "new_snr": "NA",
                    "correction_reasons": "header_aware_target_frame",
                    "validated_unrounded_flux_mjy": "1.199",
                    "validated_background_only_snr": "NA",
                },
                {
                    "source": "NGC4945",
                    "telescope": "ALMA",
                    "frequency_ghz": "230.00",
                    "date_label": "2020/01/02_1",
                    "product": "ALMA/NGC4945/b.fits",
                    "old_flux_field": "2.00",
                    "new_flux_field": "2.10",
                    "old_snr": "4.000",
                    "new_snr": "4.200",
                    "correction_reasons": (
                        "fixed_nuclear_component_source18+single_background_rms"
                    ),
                    "validated_unrounded_flux_mjy": "2.1",
                    "validated_background_only_snr": "4.200",
                },
            ]
            write_csv(manifest, list(LEGACY_CORRECTION_FIELDS), legacy)
            photometry = [
                {
                    "source": "NGC4945",
                    "telescope": "ALMA",
                    "frequency_ghz": "230.00",
                    "date_label": "2020/01/02_1",
                    "product": "ALMA/NGC4945/b.fits",
                    "old_flux_field": "2.00",
                    "new_flux_field": "2.20",
                    "old_snr": "4.000",
                    "new_snr": "5.500",
                    "unrounded_flux_mjy": "2.199",
                }
            ]
            merge_correction_manifest(manifest, photometry)
            fields, migrated = read_csv(manifest)

        self.assertEqual(tuple(fields), CORRECTION_FIELDS)
        self.assertEqual({row["schema_version"] for row in migrated}, {"2"})
        other = next(row for row in migrated if row["source"] == "OTHER")
        self.assertEqual(
            other["validated_adopted_noise_snr"],
            legacy[0]["validated_background_only_snr"],
        )
        for field in LEGACY_CORRECTION_FIELDS[:-1]:
            self.assertEqual(other[field], legacy[0][field])
        ngc = next(row for row in migrated if row["source"] == "NGC4945")
        self.assertEqual(
            ngc["correction_reasons"],
            "fixed_nuclear_component_source18+single_local_scaled_mad+"
            "source18_10bmaj_local_noise_aperture",
        )
        self.assertEqual(ngc["validated_adopted_noise_snr"], "5.500")
        self.assertNotIn("single_background_rms", ngc["correction_reasons"])
        self.assertEqual(ngc["schema_version"], CORRECTION_SCHEMA_VERSION)

    def run_check(
        self,
        *,
        photometry: Path = PHOTOMETRY,
        components: Path = COMPONENTS,
        corrections: Path = CORRECTIONS,
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--check",
                "--output",
                str(photometry),
                "--component-output",
                str(components),
                "--correction-manifest",
                str(corrections),
            ],
            text=True,
            capture_output=True,
            check=False,
        )

    def test_released_sidecars_pass_strict_validation(self) -> None:
        result = self.run_check()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_deleting_either_free_width_diagnostic_fails(self) -> None:
        fields, released_rows = read_csv(COMPONENTS)
        for removed_id in (
            "source18-free-width-0p18",
            "source18-free-width-0p24",
        ):
            with self.subTest(removed_id=removed_id), tempfile.TemporaryDirectory() as raw:
                tampered = Path(raw) / "components.csv"
                rows = [row for row in released_rows if row["record_id"] != removed_id]
                write_csv(tampered, fields, rows)
                result = self.run_check(components=tampered)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("exactly the seven expected records", result.stderr)

    def test_tampering_with_a_recorded_input_hash_fails(self) -> None:
        fields, rows = read_csv(PHOTOMETRY)
        rows[0]["input_sha256"] = "0" * 64
        with tempfile.TemporaryDirectory() as raw:
            tampered = Path(raw) / "photometry.csv"
            write_csv(tampered, fields, rows)
            result = self.run_check(photometry=tampered)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("SHA-256 mapping is stale or has been tampered", result.stderr)

    def test_tampering_with_flux_adopted_noise_snr_arithmetic_fails(self) -> None:
        fields, rows = read_csv(PHOTOMETRY)
        rows[0]["adopted_noise_mjy_per_beam"] = str(
            2.0 * float(rows[0]["adopted_noise_mjy_per_beam"])
        )
        with tempfile.TemporaryDirectory() as raw:
            tampered = Path(raw) / "photometry.csv"
            write_csv(tampered, fields, rows)
            result = self.run_check(photometry=tampered)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("adopted scaled MAD is stale", result.stderr)

    def test_retained_source18_identity_cannot_become_rejected_source17(self) -> None:
        fields, rows = read_csv(PHOTOMETRY)
        rows[0]["component_identity"] = rows[0]["rejected_seed_identity"]
        with tempfile.TemporaryDirectory() as raw:
            tampered = Path(raw) / "photometry.csv"
            write_csv(tampered, fields, rows)
            result = self.run_check(photometry=tampered)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("unexpected component_identity", result.stderr)

    def test_zero_scaled_mad_fails_closed(self) -> None:
        fields, rows = read_csv(PHOTOMETRY)
        rows[0]["noise_scaled_mad_mjy_per_beam"] = "0"
        with tempfile.TemporaryDirectory() as raw:
            tampered = Path(raw) / "photometry.csv"
            write_csv(tampered, fields, rows)
            result = self.run_check(photometry=tampered)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("noise scaled MAD must be positive", result.stderr)

    def test_source_excluded_pixel_count_cannot_exceed_aperture_count(self) -> None:
        fields, rows = read_csv(PHOTOMETRY)
        rows[0]["noise_npts"] = str(
            int(rows[0]["noise_aperture_pre_exclusion_npts"]) + 1
        )
        with tempfile.TemporaryDirectory() as raw:
            tampered = Path(raw) / "photometry.csv"
            write_csv(tampered, fields, rows)
            result = self.run_check(photometry=tampered)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("source-excluded noise_npts exceeds aperture npts", result.stderr)

    def test_manifest_adopted_snr_must_match_photometry(self) -> None:
        fields, rows = read_csv(CORRECTIONS)
        ngc = next(row for row in rows if row["source"] == "NGC4945")
        ngc["validated_adopted_noise_snr"] = "999.000"
        with tempfile.TemporaryDirectory() as raw:
            tampered = Path(raw) / "corrections.csv"
            write_csv(tampered, fields, rows)
            result = self.run_check(corrections=tampered)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("correction manifest is stale", result.stderr)

    def test_manifest_historical_fields_must_match_photometry(self) -> None:
        fields, rows = read_csv(CORRECTIONS)
        ngc = next(row for row in rows if row["source"] == "NGC4945")
        ngc["old_flux_field"] = "999.00"
        with tempfile.TemporaryDirectory() as raw:
            tampered = Path(raw) / "corrections.csv"
            write_csv(tampered, fields, rows)
            result = self.run_check(corrections=tampered)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("correction manifest is stale", result.stderr)

    def test_previous_photometry_lineage_requires_exact_unique_inventory(self) -> None:
        fields, previous_rows = read_csv(PHOTOMETRY)
        self.assertIn(tuple(fields), {LEGACY_PHOTOMETRY_FIELDS, PHOTOMETRY_FIELDS})
        legacy_rows = previous_rows[:27]
        self.assertEqual(len(legacy_rows), 27)
        current = [
            {
                "product": row["product"],
                "old_flux_field": "placeholder",
                "old_snr": "placeholder",
            }
            for row in previous_rows
        ]
        with tempfile.TemporaryDirectory() as raw:
            previous = Path(raw) / "previous.csv"
            duplicate = [dict(row) for row in legacy_rows]
            duplicate[-1]["product"] = duplicate[-2]["product"]
            write_csv(previous, fields, duplicate)
            with self.assertRaisesRegex(RuntimeError, "duplicate products"):
                preserve_pre_component_fields(current, previous, CORRECTIONS)

    def test_apply_recorded_validates_before_mutating_the_sorted_table(self) -> None:
        fields, rows = read_csv(PHOTOMETRY)
        rows[0]["input_sha256"] = "0" * 64
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            tampered = directory / "photometry.csv"
            sorted_table = directory / "fitsummary.txt"
            write_csv(tampered, fields, rows)
            shutil.copyfile(
                PACKAGE_DIR
                / "fitsumfiles/sorted/fitsummary_fixtochosen_sorted.txt",
                sorted_table,
            )
            before = sorted_table.read_bytes()
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--apply-recorded",
                    "--source-table",
                    str(sorted_table),
                    "--output",
                    str(tampered),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(sorted_table.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
