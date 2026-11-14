from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ANALYSIS_DIR.parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

import build_continuum_recovery_manifest as recovery  # noqa: E402


class ContinuumRecoveryManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.structural = recovery.build_structural_rows()
        cls.checksums = recovery._read_csv_by_id(
            recovery.CHECKSUMS, recovery.CHECKSUM_FIELDS
        )
        cls.rows = recovery.merge_manifest(cls.structural, cls.checksums)

    def test_complete_canonical_coverage(self) -> None:
        self.assertEqual(recovery.EXPECTED_ROWS, 501)
        self.assertEqual(len(self.rows), 501)
        self.assertEqual(len({row["record_id"] for row in self.rows}), 501)
        self.assertEqual(
            {row["record_id"] for row in self.rows},
            {recovery.stable_record_id(row) for row in recovery._read_final_table()},
        )

        kinds: dict[str, int] = {}
        for row in self.rows:
            kinds[row["provenance_kind"]] = kinds.get(row["provenance_kind"], 0) + 1
        self.assertEqual(
            kinds,
            {
                "alma_archive_fits": 222,
                "nvas_image_product": 248,
                "sma_casa_image": 15,
                "literature_measurement": 16,
            },
        )

    def test_manifest_exactly_preserves_product_or_citation(self) -> None:
        final_rows = recovery._read_final_table()
        self.assertEqual(
            [row["product_or_citation"] for row in self.rows],
            [row[7] for row in final_rows],
        )

    def test_all_fitted_products_have_locked_sizes_and_hashes(self) -> None:
        local = [row for row in self.rows if row["local_product_kind"] != "none"]
        literature = [row for row in self.rows if row["local_product_kind"] == "none"]
        self.assertEqual(len(local), 485)
        self.assertEqual(len(literature), 16)
        for row in local:
            self.assertGreater(int(row["local_size_bytes"]), 0)
            self.assertRegex(row["sha256"], r"^[0-9a-f]{64}$")
            self.assertEqual(row["recovery_status"], "local_product_hash_verified")
        for row in literature:
            self.assertEqual(row["local_size_bytes"], "")
            self.assertEqual(row["sha256"], "")
            self.assertEqual(row["sha256_algorithm"], "not_applicable")
            self.assertEqual(
                row["recovery_status"], "literature_only_no_local_product"
            )

        summary = recovery.build_summary(self.rows)
        self.assertEqual(summary["local_product_total_bytes"], 6_040_874_990)

    def test_release_artifacts_contain_no_workstation_paths(self) -> None:
        forbidden = ("/home/", "/Users/", "file://", "~/")
        for path in (recovery.CHECKSUMS, recovery.MANIFEST, recovery.SUMMARY):
            text = path.read_text()
            for value in forbidden:
                self.assertNotIn(value, text, f"{value} leaked into {path}")
        for row in self.rows:
            product = row["product_or_citation"]
            if row["local_product_kind"] != "none":
                self.assertFalse(Path(product).is_absolute())
                self.assertNotIn("..", Path(product).parts)

    def test_archive_specific_keys_are_present(self) -> None:
        for row in self.rows:
            if row["provenance_kind"] == "alma_archive_fits":
                self.assertIn("member_ous_uid=uid://", row["archive_key"])
                self.assertEqual(row["retrieval_url"], recovery.ALMA_URL)
            elif row["provenance_kind"] == "nvas_image_product":
                self.assertIn("program=", row["archive_key"])
                self.assertIn(";filename=", row["archive_key"])
                self.assertEqual(row["retrieval_url"], recovery.NVAS_URL)
            elif row["provenance_kind"] == "sma_casa_image":
                self.assertIn("program=2022B-H002", row["archive_key"])
                self.assertIn("date=2023-03-16", row["archive_key"])
                self.assertIn("image_path=analysis/sma_calibration/", row["archive_key"])
                self.assertEqual(row["retrieval_url"], recovery.SMA_URL)
            else:
                self.assertTrue(row["archive_key"].startswith(("doi:", "arxiv:")))
                self.assertTrue(
                    row["retrieval_url"].startswith(
                        ("https://doi.org/", "https://arxiv.org/")
                    )
                )

    def test_cycle_one_alma_product_has_resolved_member_ous(self) -> None:
        matches = [
            row
            for row in self.rows
            if row["product_filename"] == "uid___A002_X5b1929_X57d.continuum.fits"
        ]
        self.assertEqual(len(matches), 1)
        row = matches[0]
        self.assertEqual(row["observation_id"], "uid://A002/X5a9a13/X528")
        self.assertIn("project=2012.1.00912.S", row["archive_key"])
        self.assertIn("execution_block_uid=uid://A002/X5b1929/X57d", row["archive_key"])

    def test_task28_legacy_filenames_use_reviewed_member_ous_mappings(self) -> None:
        expected = {
            "uid___A002_Xa6c1df_X4e58.ms.split.cal.image.continuum.source4.image.pbcor.fits":
                "uid://A002/X7fb989/X60",
            "uid___A002_Xaac494_X32e.ms.split.cal.image.continuum.source4.image.pbcor.fits":
                "uid://A001/X122/X4dc",
            "uid___A002_X828760_X207.ms.split.cal.image.continuum.source5.pbcor.fits":
                "uid://A002/X5a9a13/Xed",
            "uid___A002_X828760_X2c2.ms.split.cal.image.continuum.source9.pbcor.fits":
                "uid://A002/X5a9a13/Xf1",
            "uid___A002_X828760_X398.ms.split.cal.image.continuum.source12.pbcor.fits":
                "uid://A002/X5a9a13/Xe9",
            "uid___A002_X828760_X454.ms.split.cal.image.continuum.source13.pbcor.fits":
                "uid://A002/X5a9a13/Xf5",
            "uid___A002_X828760_X533.ms.split.cal.image.continuum.source6.pbcor.fits":
                "uid://A002/X5a9a13/Xe5",
        }
        actual = {
            row["product_filename"]: row
            for row in self.rows
            if row["product_filename"] in expected
        }
        self.assertEqual(set(actual), set(expected))
        for filename, member_ous_uid in expected.items():
            row = actual[filename]
            self.assertEqual(row["observation_id"], member_ous_uid)
            self.assertIn(f"member_ous_uid={member_ous_uid}", row["archive_key"])
            self.assertIn("execution_block_uid=uid://A002/", row["archive_key"])

    def test_checksum_coverage_must_be_exact(self) -> None:
        shortened = dict(self.checksums)
        shortened.pop(next(iter(shortened)))
        with self.assertRaisesRegex(ValueError, "coverage does not exactly match"):
            recovery.merge_manifest(self.structural, shortened)

    def test_lightweight_check_mode_needs_no_archive(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                str(ANALYSIS_DIR / "build_continuum_recovery_manifest.py"),
                "--check",
            ],
            cwd=REPOSITORY_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("501 continuum records", result.stdout)


class DirectoryTreeHashTests(unittest.TestCase):
    @staticmethod
    def _make_tree(root: Path) -> None:
        (root / "sub/empty").mkdir(parents=True)
        (root / "table.dat").write_bytes(b"CASA table\x00\x01")
        (root / "sub/payload.bin").write_bytes(b"payload")
        os.symlink("../table.dat", root / "sub/table-link")

    def test_tree_hash_depends_on_relative_contents_not_root_or_mtime(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            first = base / "first.image"
            second = base / "elsewhere/second.image"
            self._make_tree(first)
            self._make_tree(second)
            os.utime(first / "table.dat", (1_000_000_000, 1_000_000_000))
            os.utime(second / "table.dat", (1_700_000_000, 1_700_000_000))

            first_result = recovery.hash_product(first, "directory")
            second_result = recovery.hash_product(second, "directory")
            self.assertEqual(first_result, second_result)
            self.assertEqual(first_result[2], "sha256-relpath-tree-v1")

            (second / "sub/payload.bin").write_bytes(b"changed")
            self.assertNotEqual(
                recovery.hash_product(first, "directory"),
                recovery.hash_product(second, "directory"),
            )

    def test_archive_resolution_uses_only_release_relative_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            alma = recovery._resolve_below_archive(root, "ALMA/Test/product.fits")
            sma = recovery._resolve_below_archive(
                root, "analysis/sma_calibration/track/images/source.image"
            )
            self.assertEqual(
                alma,
                (root / "analysis/multi_freq_from_archive/ALMA/Test/product.fits").resolve(),
            )
            self.assertEqual(
                sma,
                (root / "analysis/sma_calibration/track/images/source.image").resolve(),
            )
            with self.assertRaisesRegex(ValueError, "normalized release-relative"):
                recovery._resolve_below_archive(root, "../outside")


if __name__ == "__main__":
    unittest.main()
