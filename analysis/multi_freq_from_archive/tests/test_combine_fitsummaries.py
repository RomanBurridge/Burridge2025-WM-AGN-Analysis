from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "combine_fitsummaries.py"
sys.path.insert(0, str(MODULE_PATH.parent))
SPEC = importlib.util.spec_from_file_location("combine_fitsummaries", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class CombinedFitSummaryPathTests(unittest.TestCase):
    def test_combined_rows_sort_by_calendar_then_beam_then_product(self) -> None:
        def row(date_label: str, beam: str, product: str) -> list[str]:
            return [
                "NGC4258",
                "VLA",
                "4.86",
                "1.00",
                date_label,
                beam,
                "10.000",
                f"VLA/NGC4258/{product}\n",
            ]

        rows = [
            row("1997/01/06_2", "2.00", "same-day-larger-beam.imfits"),
            row("1997/01/05_1", "1.00", "january.imfits"),
            row("1997/01/06_1", "1.00", "same-day-z.imfits"),
            row("1996/12/22_1", "1.00", "december.imfits"),
            row("1997/01/06_3", "1.00", "same-day-a.imfits"),
        ]
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "combined.txt"
            MODULE.write_combined(rows, output)
            products = [item[7].strip() for item in MODULE.parse_sorted_file(output)]

        self.assertEqual(
            products,
            [
                "VLA/NGC4258/december.imfits",
                "VLA/NGC4258/january.imfits",
                "VLA/NGC4258/same-day-a.imfits",
                "VLA/NGC4258/same-day-z.imfits",
                "VLA/NGC4258/same-day-larger-beam.imfits",
            ],
        )

    def test_legacy_sma_path_is_normalized_and_unknown_absolute_fails(self) -> None:
        row = [
            "NGC1068",
            "SMA",
            "225.54",
            "49.45",
            "2023/03/16_1",
            "5.88569",
            "10.583",
            "/home/author/repo/analysis/sma_calibration/230316_128/images/clean/"
            "nitercal_10/flag2_cal1/fields/science/NGC1068_test.image\n",
        ]
        normalized = MODULE.deduplicate_rows([row])
        self.assertEqual(
            normalized[0][7],
            "analysis/sma_calibration/230316_128/images/clean/"
            "nitercal_10/flag2_cal1/fields/science/NGC1068_test.image",
        )
        row[7] = "/tmp/untraceable/NGC1068_test.image"
        with self.assertRaises(ValueError):
            MODULE.deduplicate_rows([row])

    def test_bandwidth_smearing_policy_removes_a_matching_vla_row(self) -> None:
        policy = MODULE.load_policy(
            MODULE_PATH.parent / "continuum_bandwidth_smearing_exclusion_policy.csv"
        )
        entry = next(
            item
            for item in policy
            if item.exclusion_id == "ngc4258_av0096_19850330_lband"
        )
        row = [
            entry.source,
            "VLA",
            entry.catalog_frequency_ghz,
            entry.catalog_flux_mjy,
            entry.catalog_date_label,
            entry.catalog_beam_fwhm_arcsec,
            entry.catalog_snr,
            entry.product,
        ]
        self.assertEqual(MODULE.deduplicate_rows([row], policy), [])


if __name__ == "__main__":
    unittest.main()
