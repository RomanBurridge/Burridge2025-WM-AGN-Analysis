from __future__ import annotations

import copy
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from io import StringIO
from pathlib import Path
from unittest.mock import patch


PACKAGE_DIR = Path(__file__).resolve().parents[1]
if str(PACKAGE_DIR) not in sys.path:
    sys.path.insert(0, str(PACKAGE_DIR))

import reviewed_release_evidence as release_evidence  # noqa: E402
from reviewed_release_evidence import (  # noqa: E402
    FINAL_EXACT_COUNTS,
    FROZEN_TO_RELEASE_ARTIFACTS,
    PROMOTED_REVIEW_FILENAMES,
    _expected_frozen_artifact_paths,
    main,
    promoted_review_state,
    validate_promoted_review_evidence,
)
from run_continuum_four_pass import (  # noqa: E402
    REVIEW_FREEZE_ARTIFACT_NAMES,
    sha256_file,
)


class ReviewedReleaseEvidenceTests(unittest.TestCase):
    @staticmethod
    def _write_json(path: Path, value: object) -> None:
        path.write_text(
            json.dumps(value, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    def _make_synthetic_bundle(self, root: Path) -> dict[str, object]:
        run_id = "synthetic-run"
        for _, release_relative in FROZEN_TO_RELEASE_ARTIFACTS.values():
            path = root / release_relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"synthetic {release_relative}\n", encoding="utf-8")
        runner = root / "run_continuum_four_pass.py"
        runner.write_text("# synthetic released runner\n", encoding="utf-8")
        delta_path = root / "reviewed_current_to_generated_deltas.csv"
        addition_path = root / "reviewed_new_product_additions.csv"
        delta_path.write_text("synthetic delta ledger\n", encoding="utf-8")
        addition_path.write_text("synthetic addition ledger\n", encoding="utf-8")

        expected_paths = _expected_frozen_artifact_paths(run_id)
        frozen_artifacts: dict[str, dict[str, str]] = {}
        for name in REVIEW_FREEZE_ARTIFACT_NAMES:
            if name in FROZEN_TO_RELEASE_ARTIFACTS:
                release_relative = FROZEN_TO_RELEASE_ARTIFACTS[name][1]
                digest = sha256_file(root / release_relative)
            elif name == "reviewed_delta_csv":
                digest = sha256_file(delta_path)
            elif name == "reviewed_addition_csv":
                digest = sha256_file(addition_path)
            else:
                digest = "0" * 64
            frozen_artifacts[name] = {
                "path": expected_paths[name],
                "sha256": digest,
            }
        self.assertEqual(set(frozen_artifacts), set(REVIEW_FREEZE_ARTIFACT_NAMES))

        counts = {**FINAL_EXACT_COUNTS, "reviewed_delta_rows": 2}
        freeze = {
            "schema_version": 3,
            "review_scope": "final_local_source_excluded_scaled_mad",
            "review_status": "frozen_after_exact_inventory_delta_and_identity_review",
            "run_id": run_id,
            "expected_counts": counts,
            "artifacts": frozen_artifacts,
        }
        freeze_path = root / "reviewed_local_scaled_mad_freeze.json"
        self._write_json(freeze_path, freeze)

        delta_result = {
            "schema_version": 1,
            "status": "exact_reviewed_match",
            "changed_product_count": 2,
            "reviewed_deltas": [{"product": "a"}, {"product": "b"}],
        }
        addition_result = {
            "schema_version": 1,
            "status": "exact_reviewed_additions_match",
            "reviewed_input_additions": FINAL_EXACT_COUNTS[
                "reviewed_addition_rows"
            ],
        }
        adjudication = {
            "schema_version": 2,
            "status": "exact_reviewed_match",
            "run_id": run_id,
            "review_freeze": (
                "analysis/multi_freq_from_archive/workflow_runs/"
                f"{run_id}/reviewed_local_scaled_mad_freeze.json"
            ),
            "review_freeze_sha256": sha256_file(freeze_path),
            "adjudicator_sha256": sha256_file(runner),
            "frozen_artifacts": frozen_artifacts,
            "validated_counts": counts,
            "delta_comparison": copy.deepcopy(delta_result),
            "addition_comparison": copy.deepcopy(addition_result),
            "policy": (
                "read-only adjudication; neither canonical nor generated table "
                "was modified"
            ),
        }
        self._write_json(root / "reviewed_delta_adjudication.json", adjudication)
        return {
            "freeze": freeze,
            "adjudication": adjudication,
            "delta_result": delta_result,
            "addition_result": addition_result,
        }

    def _validate_synthetic(
        self, root: Path, fixture: dict[str, object]
    ) -> dict[str, object]:
        with (
            patch.object(
                release_evidence,
                "load_reviewed_delta_rows",
                return_value=[{"product": "a"}, {"product": "b"}],
            ),
            patch.object(
                release_evidence,
                "load_reviewed_addition_rows",
                return_value=[
                    {"admission_id": str(index)}
                    for index in range(FINAL_EXACT_COUNTS["reviewed_addition_rows"])
                ],
            ),
            patch.object(
                release_evidence,
                "compare_reviewed_current_to_generated",
                return_value=fixture["delta_result"],
            ),
            patch.object(
                release_evidence,
                "compare_reviewed_product_additions",
                return_value=fixture["addition_result"],
            ),
            patch.object(
                release_evidence,
                "load_boundary_rejection_policy",
                return_value=[
                    {"product": f"VLA/Test/boundary-{index}.fits"}
                    for index in range(3)
                ],
            ),
        ):
            return validate_promoted_review_evidence(root)

    def test_release_cli_rejects_an_absent_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            stderr = StringIO()
            with redirect_stderr(stderr):
                status = main(
                    ["--analysis-dir", temporary, "--check"]
                )
            self.assertEqual(status, 1)
            self.assertIn("has not been promoted", stderr.getvalue())

    def test_partial_canonical_bundle_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / PROMOTED_REVIEW_FILENAMES[0]).write_text(
                "incomplete\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(FileNotFoundError, "partial"):
                promoted_review_state(root)

    def test_complete_but_malformed_bundle_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in PROMOTED_REVIEW_FILENAMES:
                (root / name).write_text("malformed\n", encoding="utf-8")
            self.assertEqual(promoted_review_state(root), "complete")
            with self.assertRaises(ValueError):
                validate_promoted_review_evidence(root)

    def test_exact_synthetic_bundle_passes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = self._make_synthetic_bundle(root)
            result = self._validate_synthetic(root, fixture)
            self.assertEqual(result["status"], "exact_promoted_review_evidence")

    def test_run_only_frozen_artifact_path_mutation_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = self._make_synthetic_bundle(root)
            freeze = fixture["freeze"]
            adjudication = fixture["adjudication"]
            assert isinstance(freeze, dict)
            assert isinstance(adjudication, dict)
            artifacts = freeze["artifacts"]
            assert isinstance(artifacts, dict)
            workflow_plan = artifacts["workflow_plan"]
            assert isinstance(workflow_plan, dict)
            workflow_plan["path"] = (
                "analysis/multi_freq_from_archive/workflow_runs/"
                "other-run/plan.json"
            )
            freeze_path = root / "reviewed_local_scaled_mad_freeze.json"
            self._write_json(freeze_path, freeze)
            adjudication["frozen_artifacts"] = artifacts
            adjudication["review_freeze_sha256"] = sha256_file(freeze_path)
            self._write_json(root / "reviewed_delta_adjudication.json", adjudication)
            with self.assertRaisesRegex(ValueError, "workflow_plan.*wrong run path"):
                self._validate_synthetic(root, fixture)

    def test_json_numeric_type_mutations_fail_closed(self) -> None:
        mutations = {
            "top-level bool schema": lambda value: value.__setitem__(
                "schema_version", True
            ),
            "float count": lambda value: value["validated_counts"].__setitem__(
                "canonical_rows", 409.0
            ),
            "nested bool schema": lambda value: value[
                "delta_comparison"
            ].__setitem__("schema_version", True),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                fixture = self._make_synthetic_bundle(root)
                adjudication = fixture["adjudication"]
                assert isinstance(adjudication, dict)
                mutate(adjudication)
                self._write_json(
                    root / "reviewed_delta_adjudication.json", adjudication
                )
                with self.assertRaises(ValueError):
                    self._validate_synthetic(root, fixture)

    def test_duplicate_json_field_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = self._make_synthetic_bundle(root)
            adjudication_path = root / "reviewed_delta_adjudication.json"
            original = adjudication_path.read_text(encoding="utf-8")
            adjudication_path.write_text(
                original.replace(
                    "{\n", '{\n  "schema_version": 999,\n', 1
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "Duplicate JSON object field"):
                self._validate_synthetic(root, fixture)

    def test_hash_count_and_order_mutations_fail_closed(self) -> None:
        cases = ("hash", "count", "order")
        for case in cases:
            with self.subTest(case=case), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                fixture = self._make_synthetic_bundle(root)
                adjudication = fixture["adjudication"]
                assert isinstance(adjudication, dict)
                if case == "hash":
                    generated = (
                        root
                        / "fitsumfiles/sorted/fitsummary_fixtochosen_sorted.txt"
                    )
                    generated.write_text("post-freeze mutation\n", encoding="utf-8")
                elif case == "count":
                    validated_counts = adjudication["validated_counts"]
                    assert isinstance(validated_counts, dict)
                    validated_counts["generated_rows"] = 412
                    self._write_json(
                        root / "reviewed_delta_adjudication.json", adjudication
                    )
                else:
                    delta_comparison = adjudication["delta_comparison"]
                    assert isinstance(delta_comparison, dict)
                    reviewed = delta_comparison["reviewed_deltas"]
                    assert isinstance(reviewed, list)
                    reviewed.reverse()
                    self._write_json(
                        root / "reviewed_delta_adjudication.json", adjudication
                    )
                with self.assertRaises(ValueError):
                    self._validate_synthetic(root, fixture)

    def test_canonical_bundle_is_schema_exact_and_hash_bound(self) -> None:
        # Before adjudication, the only valid staging state is that none of the
        # four canonical files exists.  As soon as promotion starts, partial or
        # malformed evidence is a hard failure; after promotion, this exercises
        # every CSV/JSON schema, cross-file hash, frozen release-side artifact,
        # and strict table-ledger recomputation.
        state = promoted_review_state(PACKAGE_DIR)
        if state == "absent":
            self.skipTest("authoritative reviewed evidence is not promoted yet")
        result = validate_promoted_review_evidence(PACKAGE_DIR)
        self.assertEqual(result["status"], "exact_promoted_review_evidence")
        self.assertEqual(
            result["generated_rows"], FINAL_EXACT_COUNTS["generated_rows"]
        )
        self.assertEqual(
            result["reviewed_addition_rows"],
            FINAL_EXACT_COUNTS["reviewed_addition_rows"],
        )
        self.assertGreater(result["reviewed_delta_rows"], 0)


if __name__ == "__main__":
    unittest.main()
