from __future__ import annotations

import csv
import math
import subprocess
import sys
import unittest
from pathlib import Path


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

import build_s230_astrometry_audit as builder  # noqa: E402


AUDIT = ANALYSIS_DIR / "machinetables/s230_anchor_astrometry_audit.csv"
S230 = ANALYSIS_DIR / "machinetables/S230_extrapolation_comparison.csv"
IDENTITY = ANALYSIS_DIR / "s230_anchor_identity_audit.csv"
WCS = ANALYSIS_DIR / "s230_wcs_evidence.csv"


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


class S230AstrometryAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.fields, cls.rows = read_csv(AUDIT)
        _, cls.s230_rows = read_csv(S230)
        _, cls.identity_rows = read_csv(IDENTITY)
        _, cls.wcs_rows = read_csv(WCS)
        cls.by_source = {row["source"]: row for row in cls.rows}
        cls.s230_by_source = {row["Source"]: row for row in cls.s230_rows}
        cls.identity_by_source = {row["source"]: row for row in cls.identity_rows}
        cls.wcs_by_source = {row["source"]: row for row in cls.wcs_rows}

    def test_exact_current_s230_coverage_order_and_fixed_schema(self) -> None:
        self.assertEqual(self.fields, list(builder.FIELDS))
        self.assertEqual(len(self.rows), 19)
        self.assertEqual(len(self.by_source), 19)
        self.assertEqual(
            [row["source"] for row in self.rows],
            [row["Source"] for row in self.s230_rows],
        )
        self.assertEqual(set(self.by_source), set(self.s230_by_source))
        self.assertEqual(set(self.by_source), set(self.identity_by_source))
        self.assertEqual(set(self.by_source), set(self.wcs_by_source))

    def test_no_stale_selected_status_ids_or_references(self) -> None:
        for source, row in self.by_source.items():
            with self.subTest(source=source):
                selected = self.s230_by_source[source]
                identity = self.identity_by_source[source]
                wcs = self.wcs_by_source[source]
                self.assertEqual(row["selected_flux_status"], selected["flux_status"])
                self.assertEqual(row["selected_observation_id"], selected["nu0_observation_id"])
                self.assertEqual(row["selected_reference"], selected["nu0_reference"])
                self.assertEqual(identity["selected_observation_id"], selected["nu0_observation_id"])
                self.assertEqual(identity["selected_reference"], selected["nu0_reference"])
                self.assertEqual(wcs["selected_observation_id"], selected["nu0_observation_id"])
                self.assertEqual(wcs["selected_reference"], selected["nu0_reference"])

    def test_detection_and_upper_limit_counts_and_coordinate_rules(self) -> None:
        detections = [row for row in self.rows if row["selected_flux_status"] == "detection"]
        upper_limits = [row for row in self.rows if row["selected_flux_status"] == "upper_limit"]
        self.assertEqual(len(detections), 12)
        self.assertEqual(len(upper_limits), 7)

        for row in detections:
            with self.subTest(source=row["source"]):
                if row["source"] == "NGC3079":
                    self.assertEqual(
                        row["component_record_status"],
                        "literature_detection_no_astrometric_component",
                    )
                    self.assertEqual(
                        row["positional_significance_status"],
                        "NOT_COMPUTABLE_NO_LOCAL_IMAGE_ASTROMETRY",
                    )
                    for field in (
                        "component_native_frame",
                        "component_ra_native_deg",
                        "component_dec_native_deg",
                        "component_common_frame",
                        "component_ra_icrs_deg",
                        "component_dec_icrs_deg",
                        "reference_native_frame",
                        "reference_ra_native_deg",
                        "reference_dec_native_deg",
                        "reference_common_frame",
                        "reference_ra_icrs_deg",
                        "reference_dec_icrs_deg",
                        "signed_east_offset_arcsec",
                        "signed_north_offset_arcsec",
                        "total_separation_arcsec",
                        "separation_geometric_beams",
                    ):
                        self.assertEqual(row[field], "NA", field)
                    continue
                self.assertEqual(row["component_record_status"], "detection_component")
                self.assertNotEqual(row["component_coordinate_evidence"], "NA")
                self.assertNotEqual(row["disposition_evidence"], "NA")
                self.assertNotEqual(row["reference_coordinate_citation"], "NA")
                self.assertEqual(row["component_common_frame"], "ICRS")
                self.assertEqual(row["reference_common_frame"], "ICRS")
                for field in (
                    "component_ra_icrs_deg",
                    "component_dec_icrs_deg",
                    "reference_ra_icrs_deg",
                    "reference_dec_icrs_deg",
                    "signed_east_offset_arcsec",
                    "signed_north_offset_arcsec",
                    "total_separation_arcsec",
                    "separation_geometric_beams",
                ):
                    self.assertTrue(math.isfinite(float(row[field])), field)

        for row in upper_limits:
            with self.subTest(source=row["source"]):
                self.assertEqual(row["component_record_status"], "upper_limit_no_component")
                for field in (
                    "component_native_frame",
                    "component_ra_native_deg",
                    "component_dec_native_deg",
                    "component_common_frame",
                    "component_ra_icrs_deg",
                    "component_dec_icrs_deg",
                    "reference_native_frame",
                    "reference_ra_native_deg",
                    "reference_dec_native_deg",
                    "reference_common_frame",
                    "reference_ra_icrs_deg",
                    "reference_dec_icrs_deg",
                    "signed_east_offset_arcsec",
                    "signed_north_offset_arcsec",
                    "total_separation_arcsec",
                    "separation_geometric_beams",
                ):
                    self.assertEqual(row[field], "NA", field)

    def test_material_dispositions_are_exact_and_evidence_backed(self) -> None:
        ambiguous = {
            row["source"] for row in self.rows if row["identity_disposition"] == "ambiguous"
        }
        blended = {
            row["source"] for row in self.rows if row["identity_disposition"] == "blended"
        }
        self.assertEqual(ambiguous, {"IC2560", "NGC1320"})
        self.assertEqual(blended, {"NGC3079", "NGC4258"})
        self.assertIn("published RA uncertainty", self.by_source["IC2560"]["disposition_evidence"])
        self.assertIn("NED is only a search coordinate", self.by_source["NGC1320"]["disposition_evidence"])
        for source in blended:
            self.assertIn("beam", self.by_source[source]["disposition_evidence"])
            self.assertIn("blend", self.by_source[source]["disposition_evidence"])

    def test_authoritative_coordinate_gaps_and_proxies_are_explicit(self) -> None:
        unavailable = {
            source
            for source, row in self.by_source.items()
            if row["selected_flux_status"] == "detection"
            and row["authoritative_nuclear_coordinate_status"].startswith("UNAVAILABLE")
        }
        self.assertEqual(unavailable, {"IC2560", "NGC1320", "NGC3079", "NGC4258"})
        self.assertEqual(
            self.by_source["NGC1194"]["reference_coordinate_role"],
            "maser_position_nuclear_proxy",
        )
        self.assertEqual(
            self.by_source["NGC4388"]["reference_coordinate_role"],
            "maser_position_nuclear_proxy",
        )
        self.assertIn(
            "relative",
            self.by_source["NGC4258"]["reference_position_uncertainty_detail"],
        )

    def test_wcs_is_complete_hash_verified_and_matches_component_frames(self) -> None:
        for row in self.rows:
            with self.subTest(source=row["source"]):
                if row["source"] == "NGC3079":
                    self.assertEqual(
                        row["wcs_evidence_status"],
                        "literature_only_no_local_product",
                    )
                    self.assertEqual(row["product_sha256"], "NA")
                    self.assertEqual(row["product_sha256_algorithm"], "not_applicable")
                    for field in builder.WCS_OUTPUT_FIELDS[4:]:
                        self.assertEqual(
                            row[field], "NOT_APPLICABLE_LITERATURE_ONLY", field
                        )
                    continue
                self.assertEqual(row["wcs_evidence_status"], "local_product_hash_verified")
                self.assertEqual(len(row["product_sha256"]), 64)
                int(row["product_sha256"], 16)
                self.assertEqual(row["wcs_projection"], "SIN")
                self.assertEqual(row["wcs_axis1_name"], "Right Ascension")
                self.assertEqual(row["wcs_axis2_name"], "Declination")
                self.assertEqual(row["wcs_axis1_unit"], "rad")
                self.assertEqual(row["wcs_axis2_unit"], "rad")
                for field in (
                    "wcs_crpix1_zero_based",
                    "wcs_crpix2_zero_based",
                    "wcs_crval1_deg",
                    "wcs_crval2_deg",
                    "wcs_cdelt1_deg_per_pixel",
                    "wcs_cdelt2_deg_per_pixel",
                    "wcs_pc11",
                    "wcs_pc12",
                    "wcs_pc21",
                    "wcs_pc22",
                    "wcs_transform_convention",
                ):
                    self.assertNotEqual(row[field], "NA", field)
                if row["selected_flux_status"] == "detection":
                    expected_frame = (
                        "J2000" if row["component_native_frame"] == "FK5_J2000" else "ICRS"
                    )
                    self.assertEqual(row["wcs_direction_reference_code"], expected_frame)

    def test_transforms_offsets_and_beam_ratios_recompute(self) -> None:
        for row in self.rows:
            if (
                row["selected_flux_status"] != "detection"
                or row["component_record_status"]
                == "literature_detection_no_astrometric_component"
            ):
                continue
            with self.subTest(source=row["source"]):
                component_ra, component_dec, _ = builder._coordinate_to_icrs(
                    row["component_native_frame"],
                    row["component_ra_native_deg"],
                    row["component_dec_native_deg"],
                )
                reference_ra, reference_dec, _ = builder._coordinate_to_icrs(
                    row["reference_native_frame"],
                    row["reference_ra_native_deg"],
                    row["reference_dec_native_deg"],
                )
                east, north, total = builder._offsets(
                    (component_ra, component_dec), (reference_ra, reference_dec)
                )
                self.assertAlmostEqual(float(row["component_ra_icrs_deg"]), component_ra, places=11)
                self.assertAlmostEqual(float(row["component_dec_icrs_deg"]), component_dec, places=11)
                self.assertAlmostEqual(float(row["reference_ra_icrs_deg"]), reference_ra, places=11)
                self.assertAlmostEqual(float(row["reference_dec_icrs_deg"]), reference_dec, places=11)
                self.assertAlmostEqual(float(row["signed_east_offset_arcsec"]), east, places=5)
                self.assertAlmostEqual(float(row["signed_north_offset_arcsec"]), north, places=5)
                self.assertAlmostEqual(float(row["total_separation_arcsec"]), total, places=5)
                ratio = total / float(row["selected_geometric_beam_arcsec"])
                self.assertAlmostEqual(float(row["separation_geometric_beams"]), ratio, places=5)

    def test_fk5_to_icrs_rotation_matches_reference_values(self) -> None:
        ra, dec, method = builder._coordinate_to_icrs(
            "FK5_J2000", "154.077958333333", "-33.563816666667"
        )
        self.assertAlmostEqual(ra, 154.07794940676828, places=12)
        self.assertAlmostEqual(dec, -33.56381197675372, places=12)
        self.assertIn("USNO_Circular_179", method)

    def test_ngc4945_source18_is_accepted_and_source17_is_rejected(self) -> None:
        row = self.by_source["NGC4945"]
        self.assertTrue(row["accepted_component_record"].endswith("#accepted-source18-fixed-canonical"))
        self.assertTrue(row["rejected_component_record"].endswith("#rejected-source17-archived-anchor"))
        self.assertEqual(row["component_identity"], "Emig2020_source18_AGN_maser_nucleus")
        self.assertEqual(row["component_ra_icrs_deg"], row["reference_ra_icrs_deg"])
        self.assertEqual(row["component_dec_icrs_deg"], row["reference_dec_icrs_deg"])
        self.assertNotEqual(row["rejected_component_ra_icrs_deg"], row["component_ra_icrs_deg"])
        self.assertIn("SSC/SNR", row["rejected_component_reason"])

    def test_task28_replacement_components_are_bound_to_their_final_fits(self) -> None:
        expected = {
            "NGC1068": ("40.669621277506", "-0.013318235271"),
            "NGC1194": ("45.954619236446", "-1.103743634277"),
        }
        for source, (ra, dec) in expected.items():
            with self.subTest(source=source):
                row = self.by_source[source]
                self.assertEqual(
                    row["component_coordinate_method"],
                    "fresh_final_fixed_beam_gaussian_centroid",
                )
                self.assertIn(
                    "230316_2_final_publication_noise_audit.csv#",
                    row["component_coordinate_evidence"],
                )
                self.assertEqual(row["component_ra_icrs_deg"], ra)
                self.assertEqual(row["component_dec_icrs_deg"], dec)

    def test_ngc5765b_uses_updated_pesce_center_and_independent_identity_support(self) -> None:
        row = self.by_source["NGC5765B"]
        self.assertEqual(row["reference_native_frame"], "ICRS")
        self.assertEqual(row["reference_ra_native_deg"], "222.714664571062286")
        self.assertEqual(row["reference_dec_native_deg"], "5.114513916805556")
        self.assertIn("10.3847/0004-637X/817/2/128", row["reference_coordinate_citation"])
        self.assertIn("10.1086/339806", row["reference_coordinate_citation"])
        self.assertIn("10.3847/2041-8213/ab75f0", row["reference_coordinate_citation"])
        self.assertEqual(
            row["reference_coordinate_evidence"],
            "Gao2016_Table2_VCS_ICRF_absolute_peak_position_plus_Pesce2020b_Table5_X0_Y0",
        )
        self.assertIn("0.0014 mas RA", row["reference_position_uncertainty_detail"])
        self.assertIn("0.0019 mas Dec", row["reference_position_uncertainty_detail"])
        self.assertIn("10.1093/mnras/stad3675", row["disposition_evidence"])
        self.assertIn("compact nuclear source", row["disposition_evidence"])
        self.assertIn("same ALMA program", row["disposition_evidence"])
        self.assertEqual(row["reference_common_frame"], "ICRS")
        self.assertEqual(row["reference_ra_icrs_deg"], "222.714664571062")
        self.assertEqual(row["reference_dec_icrs_deg"], "5.114513916806")
        self.assertEqual(row["total_separation_arcsec"], "0.055816")
        self.assertEqual(row["separation_geometric_beams"], "0.179838")

    def test_calibrator_and_qa_availability_are_never_implicit(self) -> None:
        sma = [row for row in self.rows if row["selected_telescope"] == "SMA"]
        alma = [row for row in self.rows if row["selected_telescope"] == "ALMA"]
        noema = [row for row in self.rows if row["selected_telescope"] == "NOEMA"]
        self.assertEqual(len(sma), 7)
        self.assertEqual(len(alma), 11)
        self.assertEqual(len(noema), 1)
        self.assertEqual(
            noema[0]["selected_image_phase_calibrator_status"],
            "NOT_APPLICABLE_LITERATURE_MEASUREMENT",
        )
        self.assertEqual(
            noema[0]["selected_image_astrometric_qa_status"],
            "NOT_APPLICABLE_LITERATURE_MEASUREMENT",
        )
        for row in sma:
            self.assertEqual(
                row["selected_image_phase_calibrator_status"],
                "AVAILABLE_IN_REDUCTION_CONFIGURATION",
            )
            self.assertNotEqual(row["selected_image_phase_calibrator"], "NA")
            self.assertEqual(
                row["selected_image_astrometric_qa_status"],
                "UNAVAILABLE_QUANTITATIVE_ASTROMETRIC_QA",
            )
        eso558 = self.by_source["ESO558-G009"]
        self.assertEqual(
            eso558["selected_image_phase_calibrator_status"],
            "AVAILABLE_IN_RELEASED_MEMBER_DOCUMENTATION",
        )
        self.assertEqual(eso558["selected_image_phase_calibrator"], "J0648-1744")
        self.assertIn(
            "member.uid___A001_X122_X600.README.txt#J0648_1744_ph",
            eso558["selected_image_phase_calibrator_evidence"],
        )
        self.assertEqual(
            eso558["selected_image_astrometric_qa_status"],
            "UNAVAILABLE_QUANTITATIVE_ASTROMETRIC_QA",
        )
        for row in alma:
            if row["source"] == "ESO558-G009":
                continue
            self.assertEqual(
                row["selected_image_phase_calibrator_status"],
                "UNAVAILABLE_IN_SELECTED_IMAGE_PRODUCT",
            )
            self.assertEqual(row["selected_image_phase_calibrator"], "NA")
            self.assertEqual(
                row["selected_image_astrometric_qa_status"],
                "UNAVAILABLE_IN_RELEASED_SELECTED_IMAGE",
            )

    def test_release_metadata_has_no_workstation_paths_or_blank_cells(self) -> None:
        for row in self.rows:
            with self.subTest(source=row["source"]):
                self.assertTrue(all(value != "" for value in row.values()))
                self.assertFalse(any("/home/" in value for value in row.values()))

    def test_wcs_extractor_requires_an_explicit_archive_root(self) -> None:
        source = (ANALYSIS_DIR / "extract_s230_wcs_evidence.py").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("/home/", source)
        self.assertNotIn("DEFAULT_ARCHIVE_ROOT", source)
        self.assertRegex(
            source,
            r'"--archive-root",\s*type=Path,\s*required=True,',
        )

    def test_builder_check_mode_accepts_committed_output(self) -> None:
        result = subprocess.run(
            [sys.executable, str(ANALYSIS_DIR / "build_s230_astrometry_audit.py"), "--check"],
            cwd=ANALYSIS_DIR.parents[1],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("19/19", result.stdout)


if __name__ == "__main__":
    unittest.main()
