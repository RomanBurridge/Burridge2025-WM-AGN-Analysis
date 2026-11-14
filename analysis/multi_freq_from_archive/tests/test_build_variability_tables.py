from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "build_variability_tables.py"
sys.path.insert(0, str(MODULE_PATH.parent))
SPEC = importlib.util.spec_from_file_location("build_variability_tables", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def input_row(
    source: str,
    frequency: str,
    flux: str,
    date_label: str,
    beam: str,
    snr: str,
    product: str,
) -> str:
    return repr([source, "VLA", frequency, flux, date_label, beam, snr, product + "\n"])


class VariabilityMathTests(unittest.TestCase):
    def test_absolute_symmetric_percentage_difference(self) -> None:
        self.assertAlmostEqual(MODULE.symmetric_percentage_difference(8.0, 12.0), 40.0)
        self.assertAlmostEqual(MODULE.symmetric_percentage_difference(12.0, 8.0), 40.0)

    def test_independent_formal_uncertainty(self) -> None:
        actual = MODULE.detection_contrast_uncertainty(8.0, 12.0, 1.0, 2.0)
        expected = ((12.0**2) + (16.0**2)) ** 0.5
        self.assertAlmostEqual(actual, expected)

    def test_thresholds_are_strict(self) -> None:
        self.assertTrue(MODULE.passes_selected_thresholds(1.3999, 2.5999)[0])
        self.assertFalse(MODULE.passes_selected_thresholds(1.4, 2.5999)[0])
        self.assertFalse(MODULE.passes_selected_thresholds(1.3999, 2.6)[0])
        self.assertFalse(MODULE.passes_selected_thresholds(1.0, None)[0])

    def test_sma_product_paths_become_portable_identifiers(self) -> None:
        raw = (
            "/home/author/repo/analysis/sma_calibration/230316_128/images/clean/"
            "nitercal_10/flag2_cal1/fields/science/NGC1068_test.image"
        )
        self.assertEqual(
            MODULE.portable_product_identifier("SMA", raw),
            "SMA/230316_128/clean/nitercal_10/flag2_cal1/NGC1068_test.image",
        )
        relative = raw.split("/repo/", 1)[1]
        self.assertEqual(
            MODULE.portable_product_identifier("SMA", relative),
            "SMA/230316_128/clean/nitercal_10/flag2_cal1/NGC1068_test.image",
        )


class VariabilityTableTests(unittest.TestCase):
    def _write_input(self, directory: Path) -> Path:
        path = directory / "fitsummary.txt"
        path.write_text(
            "\n".join(
                [
                    input_row("TEST", "100", "12", "2020/01/01_1", "1.0", "12", "product_A"),
                    input_row("TEST", "101", "8", "2021/01/01_1", "1.02", "8", "product_B"),
                    input_row("TEST", "100", "20*", "2022/01/01_1", "1.01", "NA", "product_C"),
                    input_row("TEST", "100", "5*", "2023/01/01_1", "1.01", "NA", "product_D"),
                    input_row("TEST", "200", "4", "2024/01/01_1", "0.19-0.50", "4", "product_E"),
                    input_row("DUP", "22.02", "17.13", "1998/01/06_1", "3.65676", "9.083", "BB0091"),
                    input_row("DUP", "22.02", "17.13", "1998/01/06_2", "3.65676", "9.154", "BB091"),
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        return path

    def _write_source_roles(self, directory: Path) -> Path:
        path = directory / "source_roles.csv"
        rows = [
            ("TEST", "TEST", "sample", "test_sample"),
            ("DUP", "DUP", "sample", "test_sample"),
        ]
        rows.extend(
            (f"SAMPLE{index:02d}", f"SAMPLE{index:02d}", "sample", "test_sample")
            for index in range(1, 20)
        )
        rows.append(("NGC5495", "NGC5495", "ancillary", "test_ancillary"))
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(MODULE.SOURCE_ROLE_COLUMNS)
            writer.writerows(rows)
        return path

    def _append_inventory_singletons(self, input_path: Path) -> None:
        with input_path.open("a", encoding="utf-8") as handle:
            for index in range(1, 20):
                handle.write(
                    input_row(
                        f"SAMPLE{index:02d}",
                        str(300 + index),
                        "1",
                        "2020/01/01_1",
                        "1.0",
                        "10",
                        f"sample_{index:02d}",
                    )
                    + "\n"
                )
            handle.write(
                input_row(
                    "NGC5495",
                    "33",
                    "0.13",
                    "NA",
                    "0.19-0.50",
                    "6.5",
                    "ancillary_product",
                )
                + "\n"
            )

    def test_censoring_ranges_same_day_and_duplicate_flags(self) -> None:
        with tempfile.TemporaryDirectory() as raw_directory:
            directory = Path(raw_directory)
            measurements = MODULE.parse_fit_summary(self._write_input(directory))
            duplicate_groups = MODULE.likely_duplicate_groups(measurements)
            pairs = MODULE.build_pairs(measurements)

            duplicate_pair = next(pair for pair in pairs if pair.first.source == "DUP")
            duplicate_row = MODULE.pair_row(duplicate_pair, duplicate_groups)
            self.assertEqual(duplicate_pair.comparison_class, "same_day_repeat")
            self.assertEqual(duplicate_row["pair_is_likely_duplicate"], "yes")

            test_measurements = [item for item in measurements if item.source == "TEST"]
            range_record = next(item for item in test_measurements if item.product == "product_E")
            self.assertEqual(range_record.beam.range_min, 0.19)
            self.assertEqual(range_record.beam.range_max, 0.50)
            self.assertIsNone(range_record.beam.geometric_mean_fwhm)

            lower_bound_pair = next(
                MODULE.build_pair(first, second)
                for first in test_measurements
                for second in test_measurements
                if first.product == "product_A" and second.product == "product_D"
            )
            self.assertEqual(lower_bound_pair.flux_contrast.relation, "lower_bound")
            self.assertIsNone(lower_bound_pair.flux_contrast.formal_independent_sigma_percent)

            uninformative_pair = next(
                MODULE.build_pair(first, second)
                for first in test_measurements
                for second in test_measurements
                if first.product == "product_A" and second.product == "product_C"
            )
            self.assertEqual(uninformative_pair.flux_contrast.relation, "uninformative")
            self.assertIsNone(uninformative_pair.flux_contrast.percent)
            self.assertFalse(uninformative_pair.selected)
            self.assertIn("censored_flux_not_quantitatively_comparable", uninformative_pair.selection_reason)

    def test_main_writes_deterministic_tables_with_headers(self) -> None:
        with tempfile.TemporaryDirectory() as raw_directory:
            directory = Path(raw_directory)
            input_path = self._write_input(directory)
            self._append_inventory_singletons(input_path)
            source_roles_path = self._write_source_roles(directory)
            outputs = {
                "variability": directory / "variability.txt",
                "correlations": directory / "correlations.txt",
                "duplicates": directory / "duplicates.txt",
                "schema": directory / "schema.txt",
            }
            arguments = [
                "--input",
                str(input_path),
                "--source-roles",
                str(source_roles_path),
                "--variability-output",
                str(outputs["variability"]),
                "--correlations-output",
                str(outputs["correlations"]),
                "--duplicates-output",
                str(outputs["duplicates"]),
                "--schema-output",
                str(outputs["schema"]),
            ]
            self.assertEqual(MODULE.main(arguments), 0)
            first_bytes = {name: path.read_bytes() for name, path in outputs.items()}
            self.assertEqual(MODULE.main(arguments), 0)
            self.assertEqual(first_bytes, {name: path.read_bytes() for name, path in outputs.items()})

            with outputs["variability"].open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, dialect="excel-tab"))
            self.assertTrue(rows)
            self.assertEqual(tuple(rows[0]), MODULE.PAIR_COLUMNS)
            self.assertTrue(any(row["comparison_class"] == "same_day_repeat" for row in rows))
            self.assertTrue(
                all(
                    row["temporal_interpretation"] == "same_day_not_temporal_variability"
                    for row in rows
                    if row["comparison_class"] == "same_day_repeat"
                )
            )
            self.assertTrue(
                all(row["selected_for_matched_pair_screening"] == "yes" for row in rows)
            )
            self.assertTrue(all(row["flux_contrast_relation"] == "exact" for row in rows))
            with outputs["correlations"].open(encoding="utf-8", newline="") as handle:
                correlation_rows = list(csv.DictReader(handle, dialect="excel-tab"))
            self.assertFalse(any(row["source"] == "NGC5495" for row in correlation_rows))


class ReleasedDataRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.input_path = MODULE_PATH.parent / "fitsumfiles" / "final" / "fitsummary_final_withextras.txt"
        cls.all_measurements = MODULE.parse_fit_summary(cls.input_path)
        cls.source_roles = MODULE.load_source_roles(
            MODULE_PATH.parent / "continuum_source_roles.csv"
        )
        cls.measurements = MODULE.select_sample_measurements(
            cls.all_measurements,
            cls.source_roles,
        )
        cls.pairs = MODULE.build_pairs(cls.measurements)
        cls.selected = [pair for pair in cls.pairs if pair.selected]
        cls.duplicates = MODULE.likely_duplicate_groups(cls.measurements)

    def test_released_input_counts_and_censoring_are_explicit(self) -> None:
        self.assertEqual(len(self.all_measurements), 501)
        self.assertEqual(len(self.measurements), 498)
        self.assertEqual(len(self.pairs), 14_193)
        self.assertEqual(len(self.selected), 105)
        self.assertEqual(sum(pair.flux_contrast.relation == "exact" for pair in self.selected), 105)
        self.assertEqual(
            sum(pair.comparison_class == "same_day_repeat" for pair in self.selected),
            2,
        )
        self.assertEqual(
            sum(pair.comparison_class == "different_day" for pair in self.selected),
            103,
        )
        self.assertFalse(any(pair.flux_contrast.relation != "exact" for pair in self.selected))
        geometry_matched_censored = [
            pair
            for pair in self.pairs
            if pair.flux_contrast.relation != "exact"
            and MODULE.passes_selected_thresholds(
                pair.frequency_difference_pct,
                pair.beam_difference_pct,
            )[0]
        ]
        self.assertEqual(len(geometry_matched_censored), 0)
        self.assertFalse(any(pair.first.source == "NGC5495" for pair in self.pairs))

    def test_ngc1320_pair_is_restored_with_distinct_frequencies(self) -> None:
        ngc1320 = [pair for pair in self.selected if pair.first.source == "NGC1320"]
        self.assertEqual(len(ngc1320), 1)
        self.assertEqual((ngc1320[0].first.frequency, ngc1320[0].second.frequency), (8.44, 8.46))

    def test_ic2560_uses_only_the_canonical_later_reduction(self) -> None:
        products = {measurement.product for measurement in self.measurements}
        self.assertIn(
            "VLA/IC2560/22.0I4.02_BB0091_1998JAN06_1_168.U56.8S.imfits",
            products,
        )
        self.assertNotIn(
            "VLA/IC2560/22.0I4.02_BB091_1998JAN06_1_168.U56.8S.imfits",
            products,
        )
        self.assertEqual(self.duplicates, {})

    def test_ngc4258_uses_doi_detections_not_failed_nvas_products(self) -> None:
        ngc4258 = [
            measurement
            for measurement in self.measurements
            if measurement.source == "NGC4258"
        ]
        products = {measurement.product for measurement in ngc4258}
        self.assertNotIn(
            "VLA/NGC4258/15.3I0.14_BR0043_1997JAN04_1_351.U42.7S.imfits",
            products,
        )
        self.assertNotIn(
            "VLA/NGC4258/22.4I0.07_AH0594_1996DEC31_1_294.U42.8S.imfits",
            products,
        )
        self.assertNotIn(
            "VLA/NGC4258/22.4I1.92_AH0847_2004MAY22_1_2.38M55.7S.imfits",
            products,
        )
        doi = [
            measurement
            for measurement in ngc4258
            if measurement.telescope == "VLA"
            and "Doi et al. 2013" in measurement.product
        ]
        self.assertEqual(len(doi), 2)
        by_frequency = {measurement.frequency: measurement for measurement in doi}
        self.assertEqual(by_frequency[15.365].flux, 3.5)
        self.assertEqual(by_frequency[22.46].flux, 3.1)
        self.assertTrue(all(not measurement.upper_limit for measurement in doi))
        self.assertAlmostEqual(by_frequency[15.365].formal_flux_sigma, 0.2)
        self.assertAlmostEqual(by_frequency[22.46].formal_flux_sigma, 0.4)

    def test_no_author_home_path_leaks_into_parsed_products(self) -> None:
        self.assertFalse(any("/home/" in item.product for item in self.measurements))
        ngc1068_sma = next(
            item
            for item in self.measurements
            if item.source == "NGC1068" and item.telescope == "SMA"
        )
        self.assertEqual(
            ngc1068_sma.product,
            "SMA/230316_128/clean/nitercal_10/flag2_cal1/NGC1068_test.image",
        )


if __name__ == "__main__":
    unittest.main()
