from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "build_fitsummary_machine.py"
sys.path.insert(0, str(MODULE_PATH.parent))
SPEC = importlib.util.spec_from_file_location("build_fitsummary_machine", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class ReleasedContinuumBuildTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.base_dir = MODULE_PATH.parent
        cls.sorted_rows = MODULE.parse_python_rows(
            cls.base_dir / "fitsumfiles/sorted/fitsummary_fixtochosen_sorted.txt"
        )
        cls.alias_policy = MODULE.load_product_alias_policy(
            cls.base_dir / "continuum_product_alias_policy.csv"
        )
        cls.exclusion_policy = MODULE.load_product_exclusion_policy(
            cls.base_dir / "continuum_product_exclusion_policy.csv"
        )
        cls.bandwidth_smearing_policy = MODULE.load_bandwidth_smearing_policy(
            cls.base_dir / "continuum_bandwidth_smearing_exclusion_policy.csv"
        )
        cls.source_roles = MODULE.load_source_roles(
            cls.base_dir / "continuum_source_roles.csv"
        )
        alias_curated_rows, cls.policy_rows = MODULE.apply_product_alias_policy(
            cls.sorted_rows,
            cls.alias_policy,
        )
        exclusion_curated_rows, cls.exclusion_rows = MODULE.apply_product_exclusion_policy(
            alias_curated_rows,
            cls.exclusion_policy,
        )
        cls.curated_rows, cls.bandwidth_smearing_rows = (
            MODULE.apply_bandwidth_smearing_policy(
                exclusion_curated_rows,
                cls.bandwidth_smearing_policy,
            )
        )
        cls.final_rows = MODULE.build_final_rows(cls.curated_rows)
        MODULE.validate_exclusion_replacements(cls.final_rows, cls.exclusion_policy)

    def test_augmented_release_counts_and_ic2560_alias_policy(self) -> None:
        self.assertEqual(
            (
                MODULE.EXPECTED_SORTED_ROWS,
                MODULE.EXPECTED_PREFILTERED_SORTED_ROWS,
                MODULE.EXPECTED_CURATED_ROWS,
                MODULE.EXPECTED_FINAL_ROWS,
            ),
            (495, 486, 485, 501),
        )
        self.assertEqual(
            (len(self.sorted_rows), len(self.curated_rows), len(self.final_rows)),
            (495, 485, 501),
        )
        MODULE.validate_release_architecture(
            self.sorted_rows,
            self.curated_rows,
            self.final_rows,
            self.exclusion_rows,
            self.bandwidth_smearing_rows,
        )
        products = {row[7] for row in self.final_rows}
        self.assertIn(
            "VLA/IC2560/22.0I4.02_BB0091_1998JAN06_1_168.U56.8S.imfits",
            products,
        )
        self.assertNotIn(
            "VLA/IC2560/22.0I4.02_BB091_1998JAN06_1_168.U56.8S.imfits",
            products,
        )

    def test_freshly_combined_prefiltered_input_has_the_same_final_architecture(self) -> None:
        excluded = {entry.product for entry in self.bandwidth_smearing_policy}
        prefiltered_rows = [
            row for row in self.sorted_rows if row[7] not in excluded
        ]
        alias_rows, _ = MODULE.apply_product_alias_policy(
            prefiltered_rows,
            self.alias_policy,
        )
        exclusion_rows, failed_matches = MODULE.apply_product_exclusion_policy(
            alias_rows,
            self.exclusion_policy,
        )
        curated_rows, smearing_matches = MODULE.apply_bandwidth_smearing_policy(
            exclusion_rows,
            self.bandwidth_smearing_policy,
        )
        final_rows = MODULE.build_final_rows(curated_rows)
        self.assertEqual(
            (len(prefiltered_rows), len(curated_rows), len(final_rows)),
            (486, 485, 501),
        )
        self.assertTrue(all(row is None for row in smearing_matches.values()))
        MODULE.validate_release_architecture(
            prefiltered_rows,
            curated_rows,
            final_rows,
            failed_matches,
            smearing_matches,
        )

    def test_nine_bandwidth_smeared_vla_products_are_removed_from_outputs(self) -> None:
        self.assertEqual(len(self.bandwidth_smearing_policy), 9)
        self.assertEqual(
            {entry.source for entry in self.bandwidth_smearing_policy},
            {"NGC3079", "NGC4258", "UGC6093"},
        )
        by_source = {}
        for entry in self.bandwidth_smearing_policy:
            by_source[entry.source] = by_source.get(entry.source, 0) + 1
            self.assertGreater(
                float(entry.smearing_parameter_x),
                float(entry.maximum_accepted_x),
            )
            self.assertIsNotNone(self.bandwidth_smearing_rows[entry.product])
        self.assertEqual(by_source, {"NGC3079": 4, "NGC4258": 4, "UGC6093": 1})
        products = {row[7] for row in self.final_rows}
        self.assertTrue(
            {entry.product for entry in self.bandwidth_smearing_policy}.isdisjoint(products)
        )

    def test_failed_ngc4258_products_are_absent_with_declared_replacements(self) -> None:
        self.assertEqual(len(self.exclusion_policy), 3)
        excluded_products = {entry.product for entry in self.exclusion_policy}
        self.assertTrue(
            all(self.exclusion_rows[product] is None for product in excluded_products)
        )
        final_products = {row[7] for row in self.final_rows}
        self.assertTrue(excluded_products.isdisjoint(final_products))

        replacements = [
            row
            for row in self.final_rows
            if row[0] == "NGC4258" and row[7] == MODULE.DOI2013_REFERENCE
        ]
        self.assertEqual(len(replacements), 2)
        by_frequency = {row[2]: row for row in replacements}
        self.assertEqual(
            by_frequency["15.365"],
            [
                "NGC4258",
                "VLA",
                "15.365",
                "3.5",
                "1997/01/04_1",
                "0.14967",
                "17.500",
                MODULE.DOI2013_REFERENCE,
            ],
        )
        self.assertEqual(
            by_frequency["22.460"],
            [
                "NGC4258",
                "VLA",
                "22.460",
                "3.1",
                "1996/12/31_1",
                "0.09487",
                "7.750",
                MODULE.DOI2013_REFERENCE,
            ],
        )

    def test_exclusion_policy_preserves_hashes_headers_and_failure_evidence(self) -> None:
        by_frequency = {
            entry.replacement_frequency_ghz: entry
            for entry in self.exclusion_policy
        }
        self.assertEqual(
            by_frequency["15.365"].sha256,
            "216232c0b622ec2c9954f8e2a0284a5b92f01251959160c6d9d03e0469a19b92",
        )
        self.assertEqual(
            by_frequency["22.460"].sha256,
            "73a6fc960a7636a6e3b45fc8ae4386c437cf3b5389326df59c303a38760099d5",
        )
        replaced = [
            entry for entry in self.exclusion_policy if entry.replacement_reference != "NA"
        ]
        self.assertEqual(len(replaced), 2)
        for entry in replaced:
            self.assertGreater(float(entry.centroid_offset_beam_major), 6.0)
            self.assertIn("exclude_failed_nvas_reduction", entry.decision_reason)
            self.assertEqual(entry.replacement_doi, "10.1088/0004-637X/765/1/63")
            self.assertEqual(entry.replacement_arxiv, "1301.4757")
            self.assertEqual(len(entry.sha256), 64)

        ah0847 = next(
            entry
            for entry in self.exclusion_policy
            if entry.exclusion_id == "ngc4258_ah0847_20040522_22ghz"
        )
        self.assertEqual(
            ah0847.sha256,
            "651244e9eb68e58c6fdc86ffb3fe6e7416c6b82080b2502194a5e377345cf27e",
        )
        self.assertEqual(ah0847.header_observer, "AH847")
        self.assertEqual(ah0847.header_date_obs, "2004-05-22")
        self.assertIn("diagonal stripes", ah0847.image_quality_evidence)
        self.assertIn("fresh CASA imfit fails", ah0847.image_quality_evidence)
        self.assertEqual(ah0847.replacement_reference, "NA")
        self.assertTrue(
            all(
                getattr(ah0847, field) == "NA"
                for field in (
                    "replacement_doi",
                    "replacement_arxiv",
                    "replacement_frequency_ghz",
                    "replacement_flux_mjy",
                    "replacement_flux_error_mjy",
                    "replacement_beam_major_arcsec",
                    "replacement_beam_minor_arcsec",
                    "replacement_date_label",
                )
            )
        )

    def test_sorted_and_final_provenance_paths_are_release_portable(self) -> None:
        released_paths = (
            self.base_dir / "fitsumfiles/sorted/fitsummary_fixtochosen_sorted.txt",
            self.base_dir / "fitsumfiles/final/fitsummary_final_withextras.txt",
        )
        for path in released_paths:
            text = path.read_text(encoding="utf-8")
            self.assertNotIn("/home/", text)
        for rows in (self.sorted_rows, self.final_rows):
            self.assertFalse(any(row[7].startswith("/") for row in rows))
            self.assertFalse(any("/home/" in row[7] for row in rows))
        sma_products = [row[7] for row in self.sorted_rows if row[1].upper() == "SMA"]
        self.assertTrue(sma_products)
        self.assertTrue(
            all(product.startswith("analysis/sma_calibration/") for product in sma_products)
        )

    def test_alias_manifest_preserves_both_products_and_header_evidence(self) -> None:
        self.assertEqual(len(self.alias_policy), 2)
        by_role = {entry.role: entry for entry in self.alias_policy}
        self.assertEqual(set(by_role), {"canonical", "excluded_alias"})
        self.assertEqual(by_role["canonical"].header_date_map, "2009-02-12")
        self.assertEqual(by_role["excluded_alias"].header_date_map, "2007-04-25")
        for entry in self.alias_policy:
            self.assertEqual(entry.header_observer, "BB091")
            self.assertEqual(entry.header_date_obs, "1998-01-06")
            self.assertEqual(entry.header_frequency_hz, "22025000000")
            self.assertEqual(entry.header_naxis1, "172")
            self.assertEqual(entry.header_naxis2, "182")
            self.assertEqual(entry.header_clean_niter, "880")
            self.assertEqual(len(entry.sha256), 64)
            self.assertIn("31289_of_31304_pixels_identical", entry.pixel_comparison_evidence)

    def test_source_role_manifest_is_exactly_21_sample_plus_ngc5495(self) -> None:
        MODULE.validate_source_roles(self.final_rows, self.source_roles)
        self.assertEqual(len(self.source_roles), 22)
        self.assertEqual(sum(role.role == "sample" for role in self.source_roles), 21)
        ancillary = [role for role in self.source_roles if role.role == "ancillary"]
        self.assertEqual([role.table_source for role in ancillary], ["NGC5495"])
        self.assertEqual(sum(row[0] == "NGC5495" for row in self.final_rows), 3)
        with (
            self.base_dir.parent / "CountingBHS_and_Spectra/bhs_source_inputs.csv"
        ).open(encoding="utf-8", newline="") as handle:
            bhs_source_keys = {
                row["source_key"] for row in csv.DictReader(handle)
            }
        manifest_sample_keys = {
            role.source_key for role in self.source_roles if role.role == "sample"
        }
        self.assertEqual(manifest_sample_keys, bhs_source_keys)

    def test_generated_continuum_products_match_checked_in_files(self) -> None:
        with tempfile.TemporaryDirectory() as raw_directory:
            directory = Path(raw_directory)
            outputs = {
                "final": directory / "fitsummary_final_withextras.txt",
                "machine": directory / "fitsummary_machine.txt",
                "aas": directory / "fitsummary_aas_mrt.txt",
                "aliases": directory / "continuum_product_aliases.csv",
                "exclusions": directory / "continuum_product_exclusions.csv",
                "bandwidth_smearing": directory / "vla_bandwidth_smearing_exclusions.csv",
            }
            MODULE.write_final_rows(self.final_rows, outputs["final"])
            MODULE.write_fixed_width_machine(self.final_rows, outputs["machine"])
            MODULE.write_aas_mrt(self.final_rows, outputs["aas"])
            MODULE.write_product_alias_report(
                self.alias_policy,
                self.policy_rows,
                outputs["aliases"],
            )
            MODULE.write_product_exclusion_report(
                self.exclusion_policy,
                self.exclusion_rows,
                outputs["exclusions"],
            )
            MODULE.write_bandwidth_smearing_report(
                self.bandwidth_smearing_policy,
                self.bandwidth_smearing_rows,
                outputs["bandwidth_smearing"],
            )

            expected = {
                "final": self.base_dir / "fitsumfiles/final/fitsummary_final_withextras.txt",
                "machine": self.base_dir / "machinetables/fitsummary_machine.txt",
                "aas": self.base_dir / "machinetables/fitsummary_aas_mrt.txt",
                "aliases": self.base_dir / "machinetables/continuum_product_aliases.csv",
                "exclusions": self.base_dir / "machinetables/continuum_product_exclusions.csv",
                "bandwidth_smearing": self.base_dir / "machinetables/vla_bandwidth_smearing_exclusions.csv",
            }
            for name, output in outputs.items():
                self.assertEqual(output.read_bytes(), expected[name].read_bytes())
            self.assertEqual(
                outputs["aas"].read_bytes(),
                (
                    self.base_dir
                    / "machinetables/AAS72796R2_datafile2_corrected.txt"
                ).read_bytes(),
            )

    def test_aas_provenance_column_covers_image_and_literature_rows(self) -> None:
        text = (self.base_dir / "machinetables/fitsummary_aas_mrt.txt").read_text()
        self.assertIn(
            "Continuum Flux-density Measurements and Product or Literature Provenance",
            text,
        )
        self.assertIn("Prov   Archive product or literature citation (5)", text)
        self.assertIn("literature-only rows give their published provenance", text)
        self.assertNotIn("File   File identifier", text)

    def test_aas_high_snr_uses_eight_bytes_without_losing_precision(self) -> None:
        row = [
            "NGC1068",
            "VLA",
            "4.86",
            "2321.97",
            "2003/02/20_1",
            "14.10000",
            "3026.080",
            "VLA/NGC1068/high_snr.imfits",
        ]
        with tempfile.TemporaryDirectory() as raw_directory:
            output = Path(raw_directory) / "high_snr_aas_mrt.txt"
            MODULE.write_aas_mrt([row], output)
            text = output.read_text(encoding="utf-8")

        self.assertIn("  66- 73 A8     ---        SNR", text)
        self.assertIn("  75-220 A146   ---        Prov", text)
        data_line = text.splitlines()[-1]
        self.assertEqual(len(data_line), 220)
        self.assertEqual(data_line[65:73], "3026.080")
        self.assertEqual(data_line[73], " ")
        self.assertTrue(data_line[74:].startswith("high_snr.imfits"))

        with self.assertRaisesRegex(ValueError, "exceeds the A8 field"):
            MODULE.aas_snr("10000.000")

    def test_correction_manuscript_machine_tables_are_exact_copies(self) -> None:
        manuscript_dir = self.base_dir.parents[2] / "revision_package/working_clean"
        if not manuscript_dir.is_dir():
            self.skipTest("adjacent manuscript revision package is not present")
        self.assertEqual(
            (self.base_dir / "machinetables/fitsummary_machine.txt").read_bytes(),
            (manuscript_dir / "tab2.txt").read_bytes(),
        )
        self.assertEqual(
            (self.base_dir / "machinetables/fitsummary_aas_mrt.txt").read_bytes(),
            (manuscript_dir / "datafile2.txt").read_bytes(),
        )


if __name__ == "__main__":
    unittest.main()
