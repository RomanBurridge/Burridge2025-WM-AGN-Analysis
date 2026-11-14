from __future__ import annotations

import csv
import hashlib
import shutil
import sys
import tempfile
import unittest
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parents[1]
if str(PACKAGE_DIR) not in sys.path:
    sys.path.insert(0, str(PACKAGE_DIR))

import continuum_boundary_rejections as boundary  # noqa: E402


POLICY_PATH = PACKAGE_DIR / "continuum_boundary_rejection_policy.csv"
EVIDENCE_PATH = PACKAGE_DIR / "continuum_boundary_rejection_evidence.csv"
EXPECTED_EVIDENCE_SHA256 = (
    "097c63b680b0c9409a41df08b4f2e345f9b524cef7b0a1aad4e7bb112c2fd210"
)
EXPECTED = {
    "VLA/NGC4258/1.51I2.67_TEST_1986JUN19_1_59.9U1.46M.imfits": {
        "input_sha256": (
            "e944a4f34b4a7da798fd083a59939b744141c9ece1ea2b01d843dd06d90d9a68"
        ),
        "input_size": "406080",
        "frequency": "1.49",
        "date": "1986/06/19",
        "boundary_percent": "91.0886498131164",
        "model_payload_sha256": (
            "e53161dac9d5b6b3bf459b00ec3f2beeb4177cff32845c33926b3f031c9671b1"
        ),
        "model_payload_size": "396800",
    },
    "VLA/NGC4258/4.89I0.82_AV0096_1985MAR30_1_48.5U1.32M.imfits": {
        "input_sha256": (
            "063524d35cd9725ef629987bca9b108b92d74e83d21ba475443eafd9b187e6d5"
        ),
        "input_size": "4302720",
        "frequency": "4.86",
        "date": "1985/03/30",
        "boundary_percent": "95.00277142484286",
        "model_payload_sha256": (
            "c643fab2ef4876b526dd485fa682788ed5605707d766909e7de685e2b7acb81e"
        ),
        "model_payload_size": "4313904",
    },
    "VLA/NGC4388/14.9I0.40_AW0122_1985MAY08_1_289.U1.43M.imfits": {
        "input_sha256": (
            "f98b15737b3df998889d505f2bb62427161ac9c094a3156ae1535c4ee6045074"
        ),
        "input_size": "11488320",
        "frequency": "14.94",
        "date": "1985/05/08",
        "boundary_percent": "93.54863463328678",
        "model_payload_sha256": (
            "6e5fb82536e88d52fc92083f1790229537a804a2b202433d34cf9bb7d9f6b8d0"
        ),
        "model_payload_size": "11498760",
    },
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class ContinuumBoundaryRejectionTests(unittest.TestCase):
    def test_canonical_policy_and_independent_evidence_pin_exact_three(self) -> None:
        self.assertEqual(
            hashlib.sha256(EVIDENCE_PATH.read_bytes()).hexdigest(),
            EXPECTED_EVIDENCE_SHA256,
        )
        policy = boundary.load_boundary_rejection_policy(POLICY_PATH)
        evidence = boundary.load_boundary_rejection_evidence(EVIDENCE_PATH)
        self.assertEqual(boundary.boundary_rejection_products(policy), set(EXPECTED))
        self.assertEqual([row["product"] for row in policy], sorted(EXPECTED))
        self.assertEqual([row["product"] for row in evidence], sorted(EXPECTED))

        evidence_by_product = {row["product"]: row for row in evidence}
        for row in policy:
            expected = EXPECTED[row["product"]]
            extracted = evidence_by_product[row["product"]]
            self.assertEqual(row["product_sha256"], expected["input_sha256"])
            self.assertEqual(row["size_bytes"], expected["input_size"])
            self.assertEqual(row["expected_frequency_ghz"], expected["frequency"])
            self.assertEqual(row["expected_observation_date"], expected["date"])
            self.assertEqual(
                row["historical_boundary_highest_percent"],
                expected["boundary_percent"],
            )
            self.assertEqual(
                row["historical_model_payload_sha256"],
                expected["model_payload_sha256"],
            )
            self.assertEqual(
                row["independent_evidence_sha256"], EXPECTED_EVIDENCE_SHA256
            )
            self.assertEqual(
                extracted["retry3_model_payload_sha256"],
                expected["model_payload_sha256"],
            )
            self.assertEqual(
                extracted["authoritative7_model_payload_sha256"],
                expected["model_payload_sha256"],
            )
            self.assertEqual(
                extracted["model_payload_size_bytes"],
                expected["model_payload_size"],
            )

    def test_policy_and_evidence_mutations_fail_closed(self) -> None:
        mutations = {
            "missing policy row": lambda policy, evidence: policy.pop(),
            "duplicate policy product": lambda policy, evidence: policy[1].__setitem__(
                "product", policy[0]["product"]
            ),
            "wrong decision": lambda policy, evidence: policy[0].__setitem__(
                "rejection_decision", "publish_anyway"
            ),
            "stale evidence status": lambda policy, evidence: evidence[0].__setitem__(
                "status", "changed"
            ),
            "payload mismatch": lambda policy, evidence: evidence[0].__setitem__(
                "authoritative7_model_payload_sha256", "0" * 64
            ),
            "wrong boundary result": lambda policy, evidence: evidence[0].__setitem__(
                "boundary_result", "false"
            ),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                policy_path = root / POLICY_PATH.name
                evidence_path = root / EVIDENCE_PATH.name
                shutil.copyfile(POLICY_PATH, policy_path)
                shutil.copyfile(EVIDENCE_PATH, evidence_path)
                policy_rows = read_csv(policy_path)
                evidence_rows = read_csv(evidence_path)
                mutate(policy_rows, evidence_rows)
                write_csv(
                    evidence_path,
                    boundary.BOUNDARY_REJECTION_EVIDENCE_FIELDS,
                    evidence_rows,
                )
                # Preserve a valid evidence hash for evidence-semantic mutations;
                # the loader must still reject the changed evidence itself.
                evidence_sha = hashlib.sha256(evidence_path.read_bytes()).hexdigest()
                for row in policy_rows:
                    row["independent_evidence_sha256"] = evidence_sha
                write_csv(
                    policy_path,
                    boundary.BOUNDARY_REJECTION_POLICY_FIELDS,
                    policy_rows,
                )
                with self.assertRaises(ValueError):
                    boundary.load_boundary_rejection_policy(policy_path)

    def test_policy_rejects_a_stale_evidence_hash_binding(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            policy_path = root / POLICY_PATH.name
            evidence_path = root / EVIDENCE_PATH.name
            shutil.copyfile(POLICY_PATH, policy_path)
            shutil.copyfile(EVIDENCE_PATH, evidence_path)
            policy_rows = read_csv(policy_path)
            for row in policy_rows:
                row["independent_evidence_sha256"] = "0" * 64
            write_csv(
                policy_path,
                boundary.BOUNDARY_REJECTION_POLICY_FIELDS,
                policy_rows,
            )
            with self.assertRaisesRegex(ValueError, "does not hash-bind"):
                boundary.load_boundary_rejection_policy(policy_path)

    def test_fresh_audit_identity_and_disposition_are_exact(self) -> None:
        policy = boundary.load_boundary_rejection_policy(POLICY_PATH)[0]
        record = {
            "source": policy["source"],
            "telescope": policy["telescope"],
            "frequency_ghz": policy["expected_frequency_ghz"],
            "product": policy["product"],
            "row_disposition": policy["expected_row_disposition"],
            "classification": policy["expected_classification"],
        }
        boundary.validate_boundary_rejection_record_identity(record, policy)
        for field, value in (
            ("source", "WrongSource"),
            ("telescope", "ALMA"),
            ("frequency_ghz", "1.50"),
            ("row_disposition", "generic_table_retained"),
            ("classification", "ordinary_upper_limit"),
        ):
            with self.subTest(field=field):
                mutated = dict(record)
                mutated[field] = value
                with self.assertRaises(RuntimeError):
                    boundary.validate_boundary_rejection_record_identity(
                        mutated, policy
                    )


if __name__ == "__main__":
    unittest.main()
