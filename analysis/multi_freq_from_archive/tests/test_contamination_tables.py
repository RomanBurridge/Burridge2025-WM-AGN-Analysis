from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import build_contamination_tables as contamination  # noqa: E402


class ContaminationTableTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = contamination.build_results()

    def selection(self, source: str, component: str, scenario: str = "standard"):
        matches = [
            row
            for row in self.results["selections"]
            if row["source"] == source
            and row["component"] == component
            and row["scenario"] == scenario
        ]
        self.assertEqual(len(matches), 1)
        return matches[0]

    def test_candidate_boundaries_are_strict(self) -> None:
        self.assertIsNone(contamination.classify_candidate(253.0))
        self.assertIsNone(contamination.classify_candidate(207.0))
        self.assertEqual(contamination.classify_candidate(math.nextafter(253.0, math.inf)), "dust")
        self.assertEqual(contamination.classify_candidate(math.nextafter(207.0, -math.inf)), "jet")

    def test_projection_formulas(self) -> None:
        dust = contamination.projection_scale("dust", 347.0, 15.0, 5.63004)
        expected_dust = (230.0 / 347.0) ** 4.5 * (5.63004 / 15.0) ** 2
        self.assertAlmostEqual(dust, expected_dust, places=15)

        jet_small_beam = contamination.projection_scale("jet", 1.49, 0.01, 5.63004)
        jet_large_beam = contamination.projection_scale("jet", 1.49, 100.0, 5.63004)
        self.assertAlmostEqual(jet_small_beam, (230.0 / 1.49) ** -0.5, places=15)
        self.assertEqual(jet_small_beam, jet_large_beam)

    def test_ngc4258_recalculated_values_and_date(self) -> None:
        dust = self.selection("NGC4258", "dust")["selection"]
        self.assertEqual(dust["record_id"], "archive-L417")
        self.assertAlmostEqual(dust["projected_flux_mjy"], 2.074305214605298, places=12)
        self.assertAlmostEqual(dust["projected_sigma_mjy"], 0.4604451086804213, places=12)
        self.assertAlmostEqual(dust["projected_fraction_s230_percent"], 16.947057358642112, places=11)

        jet = self.selection("NGC4258", "jet")["selection"]
        self.assertAlmostEqual(jet["projected_flux_mjy"], 0.14017613144215016, places=12)
        self.assertAlmostEqual(jet["projected_sigma_mjy"], 0.004305560446053081, places=15)
        self.assertEqual(jet["date"], "1997/01/07_1")
        self.assertEqual(jet["record_id"], "archive-L341")
        self.assertEqual(jet["telescope"], "VLA")

    def test_accepted_exclusion_is_absent_from_candidate_inventory(self) -> None:
        self.assertEqual(len(self.results["candidates"]), 333)
        self.assertEqual(len(self.results["selections"]), 14)
        rejected = (
            "VLA/NGC4258/"
            "22.4I1.92_AH0847_2004MAY22_1_2.38M55.7S.imfits"
        )
        self.assertFalse(
            any(row["reference"] == rejected for row in self.results["candidates"])
        )
        self.assertEqual(self.selection("NGC4258", "jet")["candidate_count"], 82)

    def test_ngc1194_detection_is_a_conservative_attribution_upper_bound(self) -> None:
        jet = self.selection("NGC1194", "jet")["selection"]
        self.assertEqual(jet["record_id"], "archive-L178")
        self.assertEqual(jet["input_status"], "detection")
        self.assertEqual(jet["constraint_status"], "attribution_upper_bound")
        self.assertAlmostEqual(jet["sigma_mjy"], 0.09534186870062653, places=15)
        self.assertAlmostEqual(
            jet["projected_sigma_mjy"], 0.007673842405753396, places=16
        )
        self.assertAlmostEqual(
            jet["projected_flux_mjy"], 0.056341350943041435, places=16
        )

    def test_ngc3079_compact_literature_result_wins(self) -> None:
        jet = self.selection("NGC3079", "jet")["selection"]
        self.assertEqual(jet["record_id"], "lit-ngc3079-sawada2000-1p4")
        self.assertAlmostEqual(jet["projected_flux_mjy"], 0.9050198172223729, places=12)
        self.assertAlmostEqual(jet["projected_sigma_mjy"], 0.18757027668490592, places=12)

    def test_circinus_is_explicitly_atca(self) -> None:
        jet = self.selection("Circinus", "jet")["selection"]
        self.assertEqual(jet["telescope"], "ATCA")
        self.assertEqual(jet["record_id"], "lit-circinus-elmouttie1998-2p37")
        self.assertEqual(jet["components"], "Central peak brightness (0.07 Jy beam^-1)")
        self.assertIn("integrated flux is blank", jet["provenance_note"])
        self.assertEqual(jet["frequency_ghz"], 2.37)
        self.assertEqual(jet["flux_mjy"], 70.0)
        self.assertEqual(jet["sigma_mjy"], 7.0)
        self.assertAlmostEqual(jet["beam_arcsec"], math.sqrt(3.4 * 3.1), places=15)
        self.assertAlmostEqual(jet["projected_flux_mjy"], 7.105723351484077, places=12)
        self.assertAlmostEqual(jet["projected_sigma_mjy"], 0.7105723351484077, places=15)
        self.assertAlmostEqual(
            jet["projected_fraction_s230_percent"], 45.84411343635344, places=14
        )
        elmouttie = [
            row
            for row in self.results["candidates"]
            if row["record_id"].startswith("lit-circinus-elmouttie1998-")
        ]
        self.assertEqual(len(elmouttie), 4)
        for candidate in elmouttie:
            self.assertIn("Central peak brightness", candidate["components"])
            self.assertNotIn("Integrated source", candidate["components"])
            self.assertIn("integrated flux is blank", candidate["provenance_note"])

    def test_target_relative_rule_changes_dust_minima(self) -> None:
        ngc1068 = self.selection("NGC1068", "dust", "central_only")["selection"]
        self.assertEqual(ngc1068["frequency_ghz"], 484.25)
        self.assertAlmostEqual(
            ngc1068["projected_flux_mjy"],
            0.0001376576449375683,
            places=16,
        )

        ngc4945 = self.selection("NGC4945", "dust")["selection"]
        self.assertEqual(ngc4945["frequency_ghz"], 258.8)
        self.assertAlmostEqual(ngc4945["projected_flux_mjy"], 0.07822131957627056, places=14)

    def test_ngc4945_jet_is_calculated_but_not_adopted(self) -> None:
        jet = self.selection("NGC4945", "jet")
        self.assertEqual(jet["adoption_status"], "not_adopted")
        self.assertIsNotNone(jet["selection"])
        self.assertAlmostEqual(
            jet["selection"]["projected_flux_mjy"],
            5.346874910108505,
            places=12,
        )

    def test_ngc1068_uses_source_reported_summed_hsa_flux(self) -> None:
        central = self.selection("NGC1068", "jet", "central_only")
        summed = self.selection("NGC1068", "jet", "summed_components")
        self.assertEqual(central["adoption_status"], "not_adopted")
        self.assertEqual(summed["adoption_status"], "adopted")
        self.assertEqual(central["selection"]["record_id"], "lit-ngc1068-gallimore2024-central")
        self.assertEqual(summed["selection"]["record_id"], "lit-ngc1068-gallimore2024-sum")
        self.assertEqual(summed["selection"]["frequency_ghz"], 21.9)
        self.assertAlmostEqual(
            summed["selection"]["beam_arcsec"],
            math.sqrt(0.00146 * 0.00054),
            places=15,
        )
        self.assertEqual(summed["selection"]["flux_mjy"], 7.0)
        self.assertEqual(summed["selection"]["sigma_mjy"], 1.0)
        self.assertAlmostEqual(central["selection"]["projected_flux_mjy"], 0.13932077898332948, places=13)
        self.assertAlmostEqual(summed["selection"]["projected_flux_mjy"], 2.160012077260922, places=12)
        self.assertAlmostEqual(summed["selection"]["projected_sigma_mjy"], 0.3085731538944174, places=13)
        self.assertIn("six Gaussians plus two point sources", summed["selection"]["components"])
        self.assertIn("Section 2.3", summed["selection"]["reference"])
        self.assertEqual(self.results["method"]["ngc1068_final_adoption"], "summed_components")
        self.assertNotIn(
            "pending_choice",
            {row["adoption_status"] for row in self.results["selections"]},
        )

    def test_modified_blackbody_20_to_2000_k_diagnostic(self) -> None:
        rows = {row["source"]: row for row in self.results["modified_blackbody"]}
        ngc4258 = rows["NGC4258"]
        self.assertEqual(ngc4258["temperature_min_k"], 20)
        self.assertEqual(ngc4258["temperature_max_k"], 2000)
        self.assertAlmostEqual(ngc4258["minimum_flux_mjy"], 3.133896537932096, places=12)
        self.assertAlmostEqual(ngc4258["maximum_flux_mjy"], 3.659452708875212, places=12)
        self.assertAlmostEqual(ngc4258["maximum_fraction_s230_percent"], 29.897699972926674, places=11)

    def test_every_selected_row_has_explicit_status_and_provenance(self) -> None:
        for item in self.results["selections"]:
            selected = item["selection"]
            if selected is None:
                self.assertEqual(item["adoption_status"], "not_available")
                continue
            for field in ("record_id", "input_status", "telescope", "date", "reference"):
                self.assertTrue(str(selected[field]).strip(), (item["source"], field))
            self.assertIn(selected["constraint_status"], {"input_upper_limit", "attribution_upper_bound"})

    def test_generated_records_do_not_leak_author_home_paths(self) -> None:
        self.assertNotIn("/home/", str(self.results))

    def test_generated_files_are_byte_deterministic_and_publish_canonical_ngc1068(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            contamination.write_results(self.results, root / "one")
            contamination.write_results(self.results, root / "two")
            files_one = sorted(path.relative_to(root / "one") for path in (root / "one").rglob("*") if path.is_file())
            files_two = sorted(path.relative_to(root / "two") for path in (root / "two").rglob("*") if path.is_file())
            self.assertEqual(files_one, files_two)
            for relative in files_one:
                self.assertEqual((root / "one" / relative).read_bytes(), (root / "two" / relative).read_bytes())
            self.assertTrue((root / "one/latex/contam_NGC1068.tex").exists())
            self.assertTrue((root / "one/latex/contam_NGC1068_central_only.tex").exists())
            self.assertFalse((root / "one/latex/contam_NGC1068_summed_components.tex").exists())

            expected_note_words = {
                "contam_Circinus.tex": 46,
                "contam_NGC1068.tex": 56,
                "contam_NGC1068_central_only.tex": 48,
                "contam_NGC1194.tex": 28,
                "contam_NGC3079.tex": 48,
                "contam_NGC4258.tex": 37,
                "contam_NGC4945.tex": 27,
            }
            for filename, expected_words in expected_note_words.items():
                table_text = (root / "one/latex" / filename).read_text()
                note = table_text.split(r"\tablecomments{", 1)[1].split(
                    "}\n\\end{deluxetable*}", 1
                )[0]
                with self.subTest(compact_note=filename):
                    self.assertEqual(len(note.split()), expected_words)
                    self.assertIn("maximum modeled contribution", note)
                    self.assertIn(r"Section~\ref{subsec:contam}", note)
                    self.assertNotIn("equality is excluded", note)

            ngc1068_tex = (root / "one/latex/contam_NGC1068.tex").read_text()
            self.assertIn(r"\label{tab:NGC1068_flux}", ngc1068_tex)
            self.assertIn(r"$7\pm1$ mJy", ngc1068_tex)
            self.assertIn(r"\citet{Gallimore2024}", ngc1068_tex)
            self.assertIn("six Gaussian components and two point sources", ngc1068_tex)
            self.assertIn(r"2.16 \pm 0.31", ngc1068_tex)
            self.assertIn(r"0.0001377 \pm 0.0000081", ngc1068_tex)
            self.assertIn(r"($9.47$--$12.1$)", ngc1068_tex)
            self.assertNotIn("e-06", ngc1068_tex)
            self.assertNotIn("Record ID}", ngc1068_tex)
            self.assertIn("Extended jet & HSA &", ngc1068_tex)
            self.assertIn("HSA combines the VLBA, GBT, and phased VLA", ngc1068_tex)
            self.assertNotIn("pending", ngc1068_tex.lower())

            ngc1068_diagnostic = (
                root / "one/latex/contam_NGC1068_central_only.tex"
            ).read_text()
            self.assertIn("central-only diagnostic; not adopted", ngc1068_diagnostic)
            self.assertIn("Extended jet [not adopted]", ngc1068_diagnostic)
            self.assertIn(r"\citet{Gallimore2024}", ngc1068_diagnostic)
            circinus_tex = (root / "one/latex/contam_Circinus.tex").read_text()
            self.assertIn("ATCA", circinus_tex)
            self.assertIn(r"$70\pm7$ mJy beam$^{-1}$ central peak", circinus_tex)
            self.assertIn(r"\citet{Elmouttie1998}", circinus_tex)
            self.assertNotIn("Integrated source", circinus_tex)
            self.assertIn(r"\label{tab:circinus_flux}", circinus_tex)

            ngc4258_tex = (root / "one/latex/contam_NGC4258.tex").read_text()
            self.assertIn(r"\label{tab:NGC4258_flux}", ngc4258_tex)
            self.assertNotIn(r"Thermal dust$^{\ast}$", ngc4258_tex)
            self.assertNotIn("2023/03/16\\_1", ngc4258_tex)
            self.assertIn("$12.2$--$12.3$", ngc4258_tex)
            self.assertIn(r"2.07 \pm 0.46", ngc4258_tex)
            self.assertIn(r"0.1402 \pm 0.0043", ngc4258_tex)

            ngc4945_tex = (root / "one/latex/contam_NGC4945.tex").read_text()
            self.assertNotIn("Extended jet", ngc4945_tex)
            self.assertIn(r"0.0782 \pm 0.0032", ngc4945_tex)
            self.assertNotIn("diagnostic only", ngc4945_tex)

    def test_publication_rounding_does_not_change_machine_precision(self) -> None:
        value, sigma, places = contamination._format_value_and_uncertainty(
            17.527674085288034, 2.485842303969371
        )
        self.assertEqual((value, sigma, places), ("17.5", "2.5", 1))
        self.assertEqual(
            contamination._format_significant_fixed(34.0256), "34.0"
        )
        self.assertEqual(contamination._format_significant_fixed(1.225), "1.23")
        self.assertEqual(
            contamination._format_value_and_uncertainty(1.225, 0.345),
            ("1.23", "0.35", 2),
        )
        self.assertEqual(contamination._format_metadata_number(42.25), "42.3")
        self.assertEqual(contamination._format_metadata_number(5.60759), "5.61")
        self.assertAlmostEqual(
            self.selection("NGC4945", "dust")["selection"]["projected_flux_mjy"],
            0.07822131957627056,
            places=14,
        )

    def test_checked_in_manuscript_tables_use_generated_numerical_rows(self) -> None:
        manuscript_dir = (
            MODULE_DIR.parents[1].parent / "revision_package" / "working_clean"
        )
        if not manuscript_dir.exists():
            self.skipTest("revision-package manuscript is not part of the Zenodo archive")

        manuscript_names = {
            "contam_NGC4258.tex": "contam_NGC4258.tex",
            "contam_NGC1194.tex": "contam_NGC1194.tex",
            "contam_NGC1068.tex": "contam_NGC1068_new.tex",
            "contam_NGC3079.tex": "contam_NGC3079_new.tex",
            "contam_NGC4945.tex": "contam_NGC4945.tex",
            "contam_Circinus.tex": "contam_Circinus_new.tex",
        }
        with tempfile.TemporaryDirectory() as temp_dir:
            generated_root = Path(temp_dir)
            contamination.write_results(self.results, generated_root)
            for generated_name, manuscript_name in manuscript_names.items():
                generated = generated_root / "latex" / generated_name
                manuscript = manuscript_dir / manuscript_name
                with self.subTest(manuscript=manuscript_name):
                    self.assertTrue(manuscript.exists(), manuscript)
                    generated_data = generated.read_text().split(
                        r"\startdata", 1
                    )[1].split(r"\enddata", 1)[0]
                    manuscript_data = manuscript.read_text().split(
                        r"\startdata", 1
                    )[1].split(r"\enddata", 1)[0]
                    self.assertEqual(
                        generated_data,
                        manuscript_data,
                        f"stale manuscript contamination rows: {manuscript}",
                    )

    def test_ngc4258_prose_uses_the_generated_extended_jet_bound(self) -> None:
        audit_root = MODULE_DIR.parents[1].parent
        manuscript = (
            audit_root / "revision_package" / "working_clean" / "contam_NGC4258.tex"
        )
        response = (
            audit_root
            / "revision_package"
            / "response"
            / "CONSOLIDATED_CORRECTION_SUMMARY.tex"
        )
        if not manuscript.exists() or not response.exists():
            self.skipTest("adjacent manuscript and response are not present")

        jet = self.selection("NGC4258", "jet")["selection"]
        value, sigma, _ = contamination._format_value_and_uncertainty(
            jet["projected_flux_mjy"], jet["projected_sigma_mjy"]
        )
        expected = rf"${value}\pm{sigma}$ mJy"
        manuscript_text = manuscript.read_text(encoding="utf-8")
        response_text = response.read_text(encoding="utf-8")
        self.assertIn(
            expected.replace(" mJy", "").replace("$", "").replace(" ", ""),
            manuscript_text.replace(" ", ""),
        )
        self.assertIn(expected, response_text)
        self.assertNotIn(r"$0.137\pm0.016$", manuscript_text)


if __name__ == "__main__":
    unittest.main()
