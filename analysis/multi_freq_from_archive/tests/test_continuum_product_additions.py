from __future__ import annotations

import csv
import hashlib
import tempfile
import unittest
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parents[1]
import sys

if str(PACKAGE_DIR) not in sys.path:
    sys.path.insert(0, str(PACKAGE_DIR))

from continuum_product_additions import (  # noqa: E402
    ADDITION_POLICY_FIELDS,
    ADDITION_POLICY_VERSION,
    ADMISSION_DECISION,
    fitter_output_source,
    load_product_addition_policy,
    publication_addition_products,
    publication_addition_products_in_fitter_order,
)
from run_continuum_four_pass import (  # noqa: E402
    HISTORICAL_CANONICAL_RELATIVE,
    REVIEW_ADDITION_FIELDS,
    REVIEW_DELTA_FIELDS,
    _atomic_write_new_csv,
    compare_current_corrected_table,
    compare_reviewed_current_to_generated,
    compare_reviewed_product_additions,
    derive_reviewed_addition_rows,
    derive_reviewed_delta_rows,
    load_reviewed_product_additions,
    read_sorted_table,
)


def write_csv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_rows(path: Path, rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{row!r}\n" for row in rows), encoding="utf-8")


def fits_card(keyword: str, value: str) -> bytes:
    return f"{keyword:<8}= {value}".ljust(80).encode("ascii")


def write_minimal_fits(path: Path, *, date_obs: str, frequency_hz: float) -> None:
    cards = [
        fits_card("SIMPLE", "                    T"),
        fits_card("BITPIX", "                    8"),
        fits_card("NAXIS", "                    0"),
        fits_card("DATE-OBS", f"'{date_obs}'"),
        fits_card("CTYPE3", "'FREQ'"),
        fits_card("CUNIT3", "'Hz'"),
        fits_card("CRVAL3", f"{frequency_hz:20.12E}"),
        b"END".ljust(80),
    ]
    payload = b"".join(cards).ljust(2880, b" ")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)


class ContinuumProductAdditionTests(unittest.TestCase):
    def test_released_additions_are_absent_from_frozen_409_row_baseline(self) -> None:
        historical = PACKAGE_DIR / HISTORICAL_CANONICAL_RELATIVE
        live = PACKAGE_DIR / "fitsumfiles/sorted/fitsummary_fixtochosen_sorted.txt"
        self.assertNotEqual(historical.resolve(), live.resolve())
        self.assertEqual(
            hashlib.sha256(historical.read_bytes()).hexdigest(),
            "c1b3f8e5f3829fab4eecc8bdf8d0e537ee77d68ab89f711189f897f2f2bdf4f8",
        )
        rows = read_sorted_table(historical)
        self.assertEqual(len(rows), 409)
        baseline_products = {str(row[7]).strip() for row in rows}
        additions = load_product_addition_policy(
            PACKAGE_DIR / "continuum_product_addition_policy.csv"
        )
        self.assertEqual(len(additions), 89)
        self.assertTrue(publication_addition_products(additions).isdisjoint(baseline_products))

    def test_released_additions_use_fitter_order_not_policy_path_order(self) -> None:
        historical = PACKAGE_DIR / HISTORICAL_CANONICAL_RELATIVE
        exclusions = PACKAGE_DIR / "continuum_product_exclusion_policy.csv"
        policy = PACKAGE_DIR / "continuum_product_addition_policy.csv"
        additions = load_product_addition_policy(policy)
        expected_order = publication_addition_products_in_fitter_order(additions)
        self.assertNotEqual(expected_order, sorted(expected_order))
        by_product = {row["product"]: row for row in additions}
        ordered_ids = [by_product[product]["admission_id"] for product in expected_order]
        self.assertEqual(len(ordered_ids), 89)
        self.assertEqual(sum(value.startswith("CPA-T28-") for value in ordered_ids), 82)
        self.assertEqual(
            [value for value in ordered_ids if not value.startswith("CPA-T28-")],
            [
                "CPA-002",
                "CPA-004",
                "CPA-003",
                "CPA-001",
                "CPA-006",
                "CPA-005",
                "CPA-007",
            ],
        )

        with exclusions.open(encoding="utf-8", newline="") as handle:
            excluded = {row["product"] for row in csv.DictReader(handle)}
        common = [
            row
            for row in read_sorted_table(historical)
            if str(row[7]) not in excluded
        ]
        generated = list(common)
        for product in expected_order:
            row = by_product[product]
            generated.append(
                [
                    fitter_output_source(row),
                    row["telescope"],
                    row["expected_frequency_ghz"],
                    "1.00",
                    row["expected_date_label"],
                    "1.00000",
                    "5.000",
                    product,
                ]
            )

        with tempfile.TemporaryDirectory() as temporary:
            generated_path = Path(temporary) / "generated.txt"
            write_rows(generated_path, generated)
            comparison = compare_current_corrected_table(
                generated_path,
                historical,
                exclusions,
                policy,
            )
            self.assertTrue(comparison["addition_product_order_matches"])

            write_rows(
                generated_path,
                common + list(reversed(generated[-len(expected_order) :])),
            )
            comparison = compare_current_corrected_table(
                generated_path,
                historical,
                exclusions,
                policy,
            )
            self.assertFalse(comparison["addition_product_order_matches"])
            self.assertFalse(comparison["product_order_matches"])

    def policy_row(
        self,
        analysis: Path,
        *,
        admission_id: str,
        source: str,
        product: str,
        frequency: str,
        date_label: str,
        disposition: str,
    ) -> dict[str, str]:
        path = analysis / product
        return {
            "schema_version": "1",
            "policy_version": ADDITION_POLICY_VERSION,
            "admission_id": admission_id,
            "source": source,
            "telescope": "ALMA",
            "expected_frequency_ghz": frequency,
            "expected_date_label": date_label,
            "product": product,
            "product_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "expected_run_disposition": disposition,
            "admission_decision": ADMISSION_DECISION,
            "reviewer": "test-reviewer",
            "review_date": "2026-08-20",
            "rationale": "Synthetic reviewed public archive product.",
        }

    def make_fixture(self, root: Path) -> dict[str, object]:
        analysis = root / "analysis"
        published = "ALMA/CircinusGalaxy/new.fits"
        input_only = "ALMA/NoRow/no-row.fits"
        write_minimal_fits(
            analysis / published,
            date_obs="2021-02-03T04:05:06.000",
            frequency_hz=500e9,
        )
        write_minimal_fits(
            analysis / input_only,
            date_obs="2022-03-04T05:06:07.000",
            frequency_hz=510e9,
        )
        policy_rows = [
            self.policy_row(
                analysis,
                admission_id="ADD-001",
                source="CircinusGalaxy",
                product=published,
                frequency="500.00",
                date_label="2021/02/03_1",
                disposition="publication_row_required",
            ),
            self.policy_row(
                analysis,
                admission_id="ADD-002",
                source="NoRow",
                product=input_only,
                frequency="510.00",
                date_label="2022/03/04_1",
                disposition="input_only_expected_no_publication_row",
            ),
        ]
        policy = root / "additions.csv"
        write_csv(policy, ADDITION_POLICY_FIELDS, policy_rows)

        canonical_rows = [
            [
                "Test",
                "ALMA",
                f"{index}.00",
                "1.00",
                "2020/01/01_1",
                "0.10000",
                "10.000",
                f"ALMA/Test/product_{index:03d}.fits",
            ]
            for index in range(409)
        ]
        generated_rows = [list(row) for row in canonical_rows[:406]]
        generated_rows[0][6] = "11.111"
        generated_rows.append(
            [
                "Circinus",
                "ALMA",
                "500.00",
                "2.00",
                "2021/02/03_1",
                "0.20000",
                "5.000",
                published,
            ]
        )
        canonical = root / "canonical.txt"
        generated = root / "generated.txt"
        exclusions = root / "exclusions.csv"
        write_rows(canonical, canonical_rows)
        write_rows(generated, generated_rows)
        write_csv(
            exclusions,
            ("product",),
            [
                {"product": f"ALMA/Test/product_{index:03d}.fits"}
                for index in range(406, 409)
            ],
        )
        return {
            "analysis": analysis,
            "policy": policy,
            "policy_rows": policy_rows,
            "canonical": canonical,
            "generated": generated,
            "generated_rows": generated_rows,
            "exclusions": exclusions,
        }

    def test_header_hash_dispositions_and_exact_ledgers_pass(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self.make_fixture(Path(temporary))
            loaded = load_product_addition_policy(
                paths["policy"], analysis_root=paths["analysis"]
            )
            self.assertEqual(len(loaded), 2)
            comparison = compare_current_corrected_table(
                paths["generated"],
                paths["canonical"],
                paths["exclusions"],
                paths["policy"],
            )
            self.assertEqual(comparison["expected_common_rows"], 406)
            self.assertEqual(comparison["expected_corrected_rows"], 407)
            self.assertEqual(comparison["addition_policy_rows"], 2)
            self.assertEqual(
                comparison["publication_addition_products"],
                ["ALMA/CircinusGalaxy/new.fits"],
            )

            delta = Path(temporary) / "delta.csv"
            additions = Path(temporary) / "addition-ledger.csv"
            write_csv(
                delta,
                REVIEW_DELTA_FIELDS,
                derive_reviewed_delta_rows(
                    paths["generated"],
                    paths["canonical"],
                    paths["exclusions"],
                    paths["policy"],
                ),
            )
            write_csv(
                additions,
                REVIEW_ADDITION_FIELDS,
                derive_reviewed_addition_rows(
                    paths["generated"], paths["canonical"], paths["policy"]
                ),
            )
            common = compare_reviewed_current_to_generated(
                paths["generated"],
                paths["canonical"],
                paths["exclusions"],
                delta,
                paths["policy"],
            )
            admitted = compare_reviewed_product_additions(
                paths["generated"],
                paths["canonical"],
                paths["policy"],
                additions,
            )
            self.assertEqual(common["validated_product_identities"], 406)
            self.assertEqual(common["validated_addition_identities"], 1)
            self.assertEqual(admitted["reviewed_input_additions"], 2)
            self.assertEqual(admitted["publication_additions"], 1)
            self.assertEqual(admitted["input_only_additions"], 1)

    def test_policy_and_inventory_mutations_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self.make_fixture(Path(temporary))
            rows = [dict(row) for row in paths["policy_rows"]]
            rows[0]["product_sha256"] = "0" * 64
            write_csv(paths["policy"], ADDITION_POLICY_FIELDS, rows)
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                load_product_addition_policy(
                    paths["policy"], analysis_root=paths["analysis"]
                )

        with tempfile.TemporaryDirectory() as temporary:
            paths = self.make_fixture(Path(temporary))
            rows = [dict(row) for row in paths["policy_rows"]]
            rows[0]["expected_frequency_ghz"] = "499.00"
            write_csv(paths["policy"], ADDITION_POLICY_FIELDS, rows)
            with self.assertRaisesRegex(ValueError, "header frequency"):
                load_product_addition_policy(
                    paths["policy"], analysis_root=paths["analysis"]
                )

        with tempfile.TemporaryDirectory() as temporary:
            paths = self.make_fixture(Path(temporary))
            rows = [dict(paths["policy_rows"][0])]
            rows[0]["source"] = "Test"
            rows[0]["expected_frequency_ghz"] = "230.00"
            rows[0]["expected_date_label"] = "2020/01/01_2"
            # Keep the path/source invariant so the collision guard, rather
            # than schema normalization, owns this mutation.
            rows[0]["product"] = "ALMA/Test/new.fits"
            write_csv(paths["policy"], ADDITION_POLICY_FIELDS, rows)
            with self.assertRaises(ValueError):
                load_reviewed_product_additions(
                    paths["policy"], canonical_rows=[
                        ["Test", "ALMA", "230.00", "1", "2020/01/01_1", "1", "1", "ALMA/Test/base.fits"]
                    ]
                )

        with tempfile.TemporaryDirectory() as temporary:
            paths = self.make_fixture(Path(temporary))
            with self.assertRaisesRegex(ValueError, "shares source/telescope"):
                load_reviewed_product_additions(
                    paths["policy"],
                    canonical_rows=[
                        [
                            "Circinus",
                            "ALMA",
                            "500.00",
                            "1",
                            "2021/02/03_1",
                            "1",
                            "1",
                            "ALMA/Circinus/base.fits",
                        ]
                    ],
                )

        with tempfile.TemporaryDirectory() as temporary:
            paths = self.make_fixture(Path(temporary))
            rows = list(paths["generated_rows"])
            rows.append(
                ["NoRow", "ALMA", "510.00", "1.00", "2022/03/04_1", "1", "4", "ALMA/NoRow/no-row.fits"]
            )
            write_rows(paths["generated"], rows)
            with self.assertRaisesRegex(ValueError, "Input-only addition"):
                derive_reviewed_addition_rows(
                    paths["generated"], paths["canonical"], paths["policy"]
                )

    def test_atomic_csv_helper_is_no_clobber(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "ledger.csv"
            _atomic_write_new_csv(path, ("field",), [{"field": "first"}])
            original = path.read_bytes()
            with self.assertRaises(FileExistsError):
                _atomic_write_new_csv(path, ("field",), [{"field": "second"}])
            self.assertEqual(path.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
