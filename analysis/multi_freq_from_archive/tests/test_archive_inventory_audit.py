from __future__ import annotations

import ast
import csv
import hashlib
import json
import os
import re
import unittest
from pathlib import Path


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
AUDIT_DIR = ANALYSIS_DIR / "archive_inventory_audit"
PRE_FIX_SHA256 = "bc898ea9d487d042b01513797d729cb49c96e52bc86a80472b4a67768ad51c22"
EXPECTED_PRODUCTS = {
    "ALMA/ESO558-G009/member.uid___A001_X122_X600.calibrated_final_cont_image-ESO558-G009.pbcor.fits":
        "5256c80efa3e68b6ea1a45d051f1a98100d0f312efe6c43a985d195891ef702d",
    "ALMA/ESO558-G009/member.uid___A001_X133d_X3fe.ESO_558-G009_sci.spw23_25_27_29.cont.I.pbcor.fits":
        "cf2a0a50a354f3a09770ddbef4a748c45fac4994859a8f028d2d96d316efdb8e",
    "ALMA/CircinusGalaxy/member.uid___A001_X1465_X11ac.ESO_97-G13_sci.spw25_27_29_31.cont.I.pbcor.fits":
        "8557e794dc78f367a6fd7540c3235dfe70ae74a722abaf14343bef74b8d382a3",
    "ALMA/CircinusGalaxy/member.uid___A001_X1590_Xe83.ESO097-013_sci.spw17_19_23_25.cont.I.pbcor.fits":
        "f8a919d1a05986df54b63c34dbfcf468fc110d628615d0211f484f274f4ecb63",
    "ALMA/CircinusGalaxy/member.uid___A001_X1590_Xe85.ESO097-013_sci.spw17_19_23_25.cont.I.pbcor.fits":
        "c22605da436bb8059b8060035c908cd5a0ed826df33dfdb101ea6feb9e631d8a",
    "ALMA/CircinusGalaxy/member.uid___A001_X1590_Xe87.ESO097-013_sci.spw16_18_20_22.cont.I.pbcor.fits":
        "09400e1356b3eee2b2be06fb00fcc951aefe964589f084bc84b10a266bdf73ce",
    "ALMA/NGC1068/member.uid___A001_X122_X4d8.NGC_1068_a_03_TE.continuum.image.pbcor.fits":
        "1f853c82efdcda4973a2e7d98d2dc46b6f0b5c2df5c72d992892c980e4fd0cf7",
}
EXPECTED_BUNDLE = {
    "EVIDENCE_SHA256SUMS.txt",
    "QUERY_PROVENANCE.md",
    "admitted_primary_products.csv",
    "archive_query_audit_executed.py",
    "get_urls_pre_boundary_fix.ipynb",
    "historical_getchooser_candidate_results.csv",
    "pre20250127_unmatched_mous_disposition.csv",
}
EXPECTED_TASK28_BUNDLE = {
    "ARCHIVE_FILE_SELECTION_AUDIT.md",
    "approved_alma_input_manifest.csv",
    "approved_alma_products_validated.csv",
    "build_approved_alma_inputs.py",
    "build_task28_identity_contact_sheets.py",
    "build_task28_identity_policy.py",
    "generate_task28_addition_policy.py",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ArchiveInventoryAuditTests(unittest.TestCase):
    def test_fail_closed_census_and_exact_admitted_products(self) -> None:
        with (AUDIT_DIR / "pre20250127_unmatched_mous_disposition.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            census = list(csv.DictReader(handle))
        with (AUDIT_DIR / "admitted_primary_products.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            admitted = list(csv.DictReader(handle))

        eligible = [row for row in census if row["eligible_pbcor_continuum"] == "yes"]
        excluded = [row for row in census if row["eligible_pbcor_continuum"] == "no"]
        self.assertEqual((len(census), len(eligible), len(excluded)), (112, 7, 105))
        self.assertTrue(all(row["datalink_expansion_status"] == "complete" for row in census))
        self.assertFalse(
            any(value.strip().lower() == "unknown" for row in census for value in row.values())
        )

        observed_products = {
            row["release_relative_path"]: row["sha256"] for row in admitted
        }
        self.assertEqual(observed_products, EXPECTED_PRODUCTS)
        self.assertEqual(
            {row["chosen_product_url"] for row in eligible},
            {row["primary_product_url"] for row in admitted},
        )

    def test_historical_and_remediated_notebooks_are_distinct_valid_json(self) -> None:
        pre_fix_path = AUDIT_DIR / "get_urls_pre_boundary_fix.ipynb"
        current_path = ANALYSIS_DIR / "get_urls.ipynb"
        pre_fix = json.loads(pre_fix_path.read_text(encoding="utf-8"))
        current = json.loads(current_path.read_text(encoding="utf-8"))
        pre_fix_source = "".join(
            line for cell in pre_fix["cells"] for line in cell.get("source", [])
        )
        current_source = "".join(
            line for cell in current["cells"] for line in cell.get("source", [])
        )

        self.assertEqual(sha256(pre_fix_path), PRE_FIX_SHA256)
        self.assertNotEqual(sha256(current_path), PRE_FIX_SHA256)
        self.assertIn("if chem in url:", pre_fix_source)
        self.assertNotIn("product_name=os.path.basename(url)", pre_fix_source)
        self.assertIn("product_name=os.path.basename(url)", current_source)
        self.assertIn("re.search(r'(^|[^A-Za-z])SO([^A-Za-z]|$)'", current_source)
        self.assertNotIn("if chem in url:", current_source)

    def test_bundle_is_portable_lightweight_and_checksum_complete(self) -> None:
        files = {path.name for path in AUDIT_DIR.iterdir() if path.is_file()}
        directories = {
            path.name
            for path in AUDIT_DIR.iterdir()
            if path.is_dir() and path.name != "__pycache__"
        }
        self.assertEqual(files, EXPECTED_BUNDLE)
        self.assertEqual(directories, {"task28"})
        self.assertEqual(
            {
                path.name
                for path in (AUDIT_DIR / "task28").iterdir()
                if path.is_file()
            },
            EXPECTED_TASK28_BUNDLE,
        )
        self.assertFalse(any(path.suffix.lower() in {".fit", ".fits"} for path in AUDIT_DIR.rglob("*")))

        forbidden_workstation_path = re.compile(
            r"(?:/home/|/Users/|/Volumes/|/tmp/|[A-Za-z]:[\\\\/]Users[\\\\/])"
        )
        for path in AUDIT_DIR.rglob("*"):
            if path.is_file() and "__pycache__" not in path.parts:
                self.assertIsNone(
                    forbidden_workstation_path.search(path.read_text(encoding="utf-8")),
                    path.relative_to(AUDIT_DIR).as_posix(),
                )

        script = AUDIT_DIR / "archive_query_audit_executed.py"
        ast.parse(script.read_text(encoding="utf-8"), filename=str(script))
        self.assertTrue(os.access(script, os.X_OK))
        self.assertTrue(script.read_text(encoding="utf-8").startswith("#!/usr/bin/env python3"))

        checksum_rows = {}
        for line in (AUDIT_DIR / "EVIDENCE_SHA256SUMS.txt").read_text(
            encoding="utf-8"
        ).splitlines():
            digest, name = line.split("  ", 1)
            checksum_rows[name] = digest
        expected_checked = EXPECTED_BUNDLE - {"EVIDENCE_SHA256SUMS.txt"}
        self.assertEqual(set(checksum_rows), expected_checked)
        self.assertEqual(
            checksum_rows,
            {name: sha256(AUDIT_DIR / name) for name in expected_checked},
        )


if __name__ == "__main__":
    unittest.main()
