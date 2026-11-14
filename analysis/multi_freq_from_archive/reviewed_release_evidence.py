#!/usr/bin/env python3
"""Validate the canonical, promoted four-pass review-evidence bundle.

The authoritative workflow creates review evidence inside an isolated,
run-specific ``workflow_runs`` directory.  A release promotes byte-identical
copies of the four review artifacts to this module's directory.  This checker
binds those copies back to one another and to the promoted science-side
artifacts without requiring the large CASA workspace.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Sequence

from run_continuum_four_pass import (
    HISTORICAL_CANONICAL_RELATIVE,
    REVIEW_FREEZE_ARTIFACT_NAMES,
    compare_reviewed_current_to_generated,
    compare_reviewed_product_additions,
    load_review_freeze,
    load_reviewed_addition_rows,
    load_reviewed_delta_rows,
    sha256_file,
)
from continuum_boundary_rejections import load_boundary_rejection_policy


PROMOTED_REVIEW_FILENAMES = (
    "reviewed_current_to_generated_deltas.csv",
    "reviewed_new_product_additions.csv",
    "reviewed_local_scaled_mad_freeze.json",
    "reviewed_delta_adjudication.json",
)

ADJUDICATION_FIELDS = (
    "schema_version",
    "status",
    "run_id",
    "review_freeze",
    "review_freeze_sha256",
    "adjudicator_sha256",
    "frozen_artifacts",
    "validated_counts",
    "delta_comparison",
    "addition_comparison",
    "policy",
)

FINAL_EXACT_COUNTS = {
    "canonical_rows": 409,
    "reviewed_exclusion_rows": 3,
    "generated_rows": 495,
    "generic_noise_audit_rows": 498,
    "generic_publication_noise_audit_rows": 495,
    "boundary_rejection_policy_rows": 3,
    "boundary_rejected_noise_audit_rows": 3,
    "ngc4945_photometry_rows": 45,
    "ngc4945_component_rows": 7,
    "preapplication_correction_rows": 66,
    "postrun_correction_rows": 84,
    "classification_identity_rows": 100,
    "addition_policy_rows": 89,
    "publication_addition_rows": 89,
    "reviewed_addition_rows": 89,
}

_RUN_ID_PATTERN = re.compile(r"[A-Za-z0-9_.-]+")
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
_WORKSPACE_PREFIX = "analysis/multi_freq_from_archive/"

# Each tuple is (path recorded in the run freeze, canonical release copy).
# The generated table is deliberately renamed on promotion; all other entries
# retain their workspace-relative paths beneath the analysis directory.
FROZEN_TO_RELEASE_ARTIFACTS = {
    "canonical_table": (
        _WORKSPACE_PREFIX + HISTORICAL_CANONICAL_RELATIVE.as_posix(),
        HISTORICAL_CANONICAL_RELATIVE,
    ),
    "generated_table": (
        _WORKSPACE_PREFIX
        + "fitsumfiles/sorted/230316_2_fitsummary_fixedtochosen_sorted.txt",
        Path("fitsumfiles/sorted/fitsummary_fixtochosen_sorted.txt"),
    ),
    "exclusion_policy": (
        _WORKSPACE_PREFIX + "continuum_product_exclusion_policy.csv",
        Path("continuum_product_exclusion_policy.csv"),
    ),
    "addition_policy": (
        _WORKSPACE_PREFIX + "continuum_product_addition_policy.csv",
        Path("continuum_product_addition_policy.csv"),
    ),
    "boundary_rejection_policy": (
        _WORKSPACE_PREFIX + "continuum_boundary_rejection_policy.csv",
        Path("continuum_boundary_rejection_policy.csv"),
    ),
    "generic_final_noise_audit": (
        _WORKSPACE_PREFIX
        + "fitsumfiles/noise_audits/230316_2_final_publication_noise_audit.csv",
        Path("fitsumfiles/noise_audits/230316_2_final_publication_noise_audit.csv"),
    ),
    "ngc4945_nuclear_photometry": (
        _WORKSPACE_PREFIX + "ngc4945_nuclear_photometry.csv",
        Path("ngc4945_nuclear_photometry.csv"),
    ),
    "ngc4945_component_audit": (
        _WORKSPACE_PREFIX + "ngc4945_component_audit.csv",
        Path("ngc4945_component_audit.csv"),
    ),
    "classification_identity_policy": (
        _WORKSPACE_PREFIX + "continuum_classification_identity_policy.csv",
        Path("continuum_classification_identity_policy.csv"),
    ),
    "postrun_cumulative_corrections": (
        _WORKSPACE_PREFIX + "validated_continuum_corrections.csv",
        Path("validated_continuum_corrections.csv"),
    ),
}


def _expected_frozen_artifact_paths(run_id: str) -> dict[str, str]:
    """Return every deterministic workspace path emitted by the runner."""

    evidence = f"{_WORKSPACE_PREFIX}workflow_runs/{run_id}/"
    paths = {
        "workflow_plan": evidence + "plan.json",
        "prepared_archive_inventory": ".continuum-fit-inputs.json",
        "prepared_input_integrity": evidence + "prepared_input_integrity.json",
        "initial_comparison": evidence + "current_corrected_table_comparison.json",
        "final_stage_snapshot": (
            evidence + "ngc4945_fixed_nucleus_check.artifacts.json"
        ),
        "canonical_table": (
            _WORKSPACE_PREFIX + HISTORICAL_CANONICAL_RELATIVE.as_posix()
        ),
        "generated_table": (
            _WORKSPACE_PREFIX
            + "fitsumfiles/sorted/230316_2_fitsummary_fixedtochosen_sorted.txt"
        ),
        "exclusion_policy": (
            _WORKSPACE_PREFIX + "continuum_product_exclusion_policy.csv"
        ),
        "addition_policy": (
            _WORKSPACE_PREFIX + "continuum_product_addition_policy.csv"
        ),
        "boundary_rejection_policy": (
            _WORKSPACE_PREFIX + "continuum_boundary_rejection_policy.csv"
        ),
        "generic_final_noise_audit": (
            _WORKSPACE_PREFIX
            + "fitsumfiles/noise_audits/"
            "230316_2_final_publication_noise_audit.csv"
        ),
        "ngc4945_nuclear_photometry": (
            _WORKSPACE_PREFIX + "ngc4945_nuclear_photometry.csv"
        ),
        "ngc4945_component_audit": (
            _WORKSPACE_PREFIX + "ngc4945_component_audit.csv"
        ),
        "final_noise_inventory_validation": (
            evidence + "final_noise_inventory_validation.json"
        ),
        "correction_manifest_transition_validation": (
            evidence + "correction_manifest_transition_validation.json"
        ),
        "classification_identity_policy": (
            _WORKSPACE_PREFIX + "continuum_classification_identity_policy.csv"
        ),
        "classification_identity_validation": (
            evidence + "classification_identity_validation.json"
        ),
        "analysis_utils_snapshot": evidence + "analysis_utils_snapshot.json",
        "preapplication_cumulative_corrections": (
            evidence + "preapplication_validated_continuum_corrections.csv"
        ),
        "postrun_cumulative_corrections": (
            _WORKSPACE_PREFIX + "validated_continuum_corrections.csv"
        ),
        "reviewed_delta_csv": (
            evidence + "reviewed_current_to_generated_deltas.csv"
        ),
        "reviewed_addition_csv": (
            evidence + "reviewed_new_product_additions.csv"
        ),
    }
    if set(paths) != set(REVIEW_FREEZE_ARTIFACT_NAMES):
        raise AssertionError("Internal frozen-artifact path map is incomplete")
    return paths


def promoted_review_paths(analysis_dir: Path) -> dict[str, Path]:
    """Return the four fixed canonical destinations."""

    return {name: analysis_dir / name for name in PROMOTED_REVIEW_FILENAMES}


def promoted_review_state(analysis_dir: Path) -> str:
    """Return ``absent`` or ``complete`` and reject every partial state."""

    paths = promoted_review_paths(analysis_dir)
    present = {
        name
        for name, path in paths.items()
        if path.exists() or path.is_symlink()
    }
    if not present:
        return "absent"
    expected = set(paths)
    if present != expected:
        raise FileNotFoundError(
            "Promoted review evidence is partial; "
            f"present={sorted(present)}, missing={sorted(expected - present)}"
        )
    for name, path in paths.items():
        if path.is_symlink() or not path.is_file():
            raise ValueError(
                f"Promoted review artifact must be a real regular file: {name}: {path}"
            )
    return "complete"


def _exact_object(value: object, fields: Sequence[str], label: str) -> dict[str, object]:
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object")
    expected = set(fields)
    actual = set(value)
    if actual != expected:
        raise ValueError(
            f"{label} fields differ from the exact schema; "
            f"missing={sorted(expected - actual)}, extra={sorted(actual - expected)}"
        )
    return value


def _require_hash(value: object, label: str) -> str:
    if not isinstance(value, str) or not _SHA256_PATTERN.fullmatch(value):
        raise ValueError(f"{label} must be 64 lowercase hexadecimal characters")
    return value


def _unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"Duplicate JSON object field: {key!r}")
        value[key] = item
    return value


def _reject_nonfinite_json(value: str) -> object:
    raise ValueError(f"Non-finite JSON value is forbidden: {value}")


def _load_exact_json(path: Path, label: str) -> object:
    """Load standards-compliant JSON while rejecting duplicate object fields."""

    try:
        return json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_unique_json_object,
            parse_constant=_reject_nonfinite_json,
        )
    except ValueError as exc:
        raise ValueError(f"{label} is not exact JSON: {exc}") from exc


def _strict_json_equal(actual: object, expected: object) -> bool:
    """Compare JSON values without Python's bool/int/float equivalence."""

    if type(actual) is not type(expected):
        return False
    if isinstance(actual, dict):
        assert isinstance(expected, dict)
        return set(actual) == set(expected) and all(
            _strict_json_equal(actual[key], expected[key]) for key in actual
        )
    if isinstance(actual, list):
        assert isinstance(expected, list)
        return len(actual) == len(expected) and all(
            _strict_json_equal(left, right)
            for left, right in zip(actual, expected)
        )
    return actual == expected


def _require_frozen_release_artifacts(
    analysis_dir: Path,
    frozen_artifacts: dict[str, object],
) -> None:
    for name, artifact_paths in FROZEN_TO_RELEASE_ARTIFACTS.items():
        expected_frozen_path, release_relative = artifact_paths
        record = _exact_object(
            frozen_artifacts[name], ("path", "sha256"), f"frozen artifact {name}"
        )
        if record["path"] != expected_frozen_path:
            raise ValueError(
                f"Frozen artifact {name} has the wrong run path: {record['path']!r}"
            )
        expected_hash = _require_hash(
            record["sha256"], f"frozen artifact {name} SHA-256"
        )
        release_path = analysis_dir / release_relative
        if release_path.is_symlink() or not release_path.is_file():
            raise FileNotFoundError(
                f"Frozen release-side artifact is absent or not a real file: {release_path}"
            )
        actual_hash = sha256_file(release_path)
        if actual_hash != expected_hash:
            raise ValueError(
                f"Promoted {name} differs from the frozen run artifact: "
                f"expected={expected_hash}, found={actual_hash}"
            )


def validate_promoted_review_evidence(analysis_dir: Path) -> dict[str, object]:
    """Fail closed unless the canonical review bundle is exact and coherent."""

    analysis_dir = analysis_dir.resolve()
    if promoted_review_state(analysis_dir) != "complete":
        raise FileNotFoundError(
            "Canonical reviewed evidence has not been promoted to "
            f"{analysis_dir}"
        )
    paths = promoted_review_paths(analysis_dir)

    delta_rows = load_reviewed_delta_rows(
        paths["reviewed_current_to_generated_deltas.csv"]
    )
    addition_rows = load_reviewed_addition_rows(
        paths["reviewed_new_product_additions.csv"]
    )

    raw_freeze = _load_exact_json(
        paths["reviewed_local_scaled_mad_freeze.json"],
        "promoted review freeze",
    )
    if not isinstance(raw_freeze, dict):
        raise ValueError("Promoted review freeze must be a JSON object")
    run_id = raw_freeze.get("run_id")
    if not isinstance(run_id, str) or not _RUN_ID_PATTERN.fullmatch(run_id):
        raise ValueError("Promoted review freeze has an unsafe or missing run_id")
    freeze_hash = sha256_file(paths["reviewed_local_scaled_mad_freeze.json"])
    freeze = load_review_freeze(
        paths["reviewed_local_scaled_mad_freeze.json"], freeze_hash, run_id
    )
    counts = freeze["expected_counts"]
    assert isinstance(counts, dict)
    for field, expected in FINAL_EXACT_COUNTS.items():
        if counts[field] != expected:
            raise ValueError(
                f"Final review-freeze count {field} must equal {expected}, "
                f"found {counts[field]}"
            )
    if counts["reviewed_delta_rows"] != len(delta_rows):
        raise ValueError("Review-freeze delta count differs from promoted delta CSV")
    if counts["reviewed_delta_rows"] < 1:
        raise ValueError("Final reviewed-delta ledger must not be empty")
    if counts["reviewed_addition_rows"] != len(addition_rows):
        raise ValueError(
            "Review-freeze addition count differs from promoted addition CSV"
        )

    frozen_artifacts = freeze["artifacts"]
    assert isinstance(frozen_artifacts, dict)
    if set(frozen_artifacts) != set(REVIEW_FREEZE_ARTIFACT_NAMES):
        raise ValueError("Promoted review freeze has a noncanonical artifact inventory")
    expected_frozen_paths = _expected_frozen_artifact_paths(run_id)
    for frozen_name, expected_path in expected_frozen_paths.items():
        record = _exact_object(
            frozen_artifacts[frozen_name],
            ("path", "sha256"),
            f"frozen artifact {frozen_name}",
        )
        if record["path"] != expected_path:
            raise ValueError(
                f"Frozen artifact {frozen_name} has the wrong run path: "
                f"{record['path']!r}"
            )
        _require_hash(record["sha256"], f"frozen artifact {frozen_name} SHA-256")
    ledger_bindings = {
        "reviewed_delta_csv": "reviewed_current_to_generated_deltas.csv",
        "reviewed_addition_csv": "reviewed_new_product_additions.csv",
    }
    for frozen_name, canonical_name in ledger_bindings.items():
        record = _exact_object(
            frozen_artifacts[frozen_name],
            ("path", "sha256"),
            f"frozen artifact {frozen_name}",
        )
        expected_path = expected_frozen_paths[frozen_name]
        if record["path"] != expected_path:
            raise ValueError(
                f"Frozen {frozen_name} path does not identify the adjudicated run"
            )
        if _require_hash(record["sha256"], f"{frozen_name} SHA-256") != sha256_file(
            paths[canonical_name]
        ):
            raise ValueError(
                f"Promoted {canonical_name} differs from its frozen run copy"
            )
    _require_frozen_release_artifacts(analysis_dir, frozen_artifacts)
    boundary_rows = load_boundary_rejection_policy(
        analysis_dir / "continuum_boundary_rejection_policy.csv"
    )
    if len(boundary_rows) != counts["boundary_rejection_policy_rows"]:
        raise ValueError(
            "Promoted boundary-rejection policy count differs from the freeze"
        )

    adjudication = _exact_object(
        _load_exact_json(
            paths["reviewed_delta_adjudication.json"],
            "promoted adjudication",
        ),
        ADJUDICATION_FIELDS,
        "promoted adjudication",
    )
    if type(adjudication["schema_version"]) is not int or adjudication[
        "schema_version"
    ] != 2:
        raise ValueError("Promoted adjudication schema_version must equal 2")
    if adjudication["status"] != "exact_reviewed_match":
        raise ValueError("Promoted adjudication status is not exact_reviewed_match")
    if adjudication["run_id"] != run_id:
        raise ValueError("Adjudication and review-freeze run IDs differ")
    expected_freeze_path = (
        f"{_WORKSPACE_PREFIX}workflow_runs/{run_id}/"
        "reviewed_local_scaled_mad_freeze.json"
    )
    if adjudication["review_freeze"] != expected_freeze_path:
        raise ValueError("Adjudication review_freeze path is not run-specific and exact")
    if adjudication["review_freeze_sha256"] != freeze_hash:
        raise ValueError("Adjudication does not hash-bind the promoted review freeze")
    if not _strict_json_equal(adjudication["frozen_artifacts"], frozen_artifacts):
        raise ValueError("Adjudication frozen-artifact map differs from the review freeze")
    if not _strict_json_equal(adjudication["validated_counts"], counts):
        raise ValueError("Adjudication counts differ from the review freeze")
    if adjudication["policy"] != (
        "read-only adjudication; neither canonical nor generated table was modified"
    ):
        raise ValueError("Adjudication policy statement differs from the exact contract")
    runner_hash = sha256_file(analysis_dir / "run_continuum_four_pass.py")
    if _require_hash(
        adjudication["adjudicator_sha256"], "adjudicator SHA-256"
    ) != runner_hash:
        raise ValueError("Adjudication was not produced by the released runner bytes")

    generated = analysis_dir / "fitsumfiles/sorted/fitsummary_fixtochosen_sorted.txt"
    canonical = analysis_dir / HISTORICAL_CANONICAL_RELATIVE
    exclusion_policy = analysis_dir / "continuum_product_exclusion_policy.csv"
    addition_policy = analysis_dir / "continuum_product_addition_policy.csv"
    fresh_delta = compare_reviewed_current_to_generated(
        generated,
        canonical,
        exclusion_policy,
        paths["reviewed_current_to_generated_deltas.csv"],
        addition_policy,
    )
    fresh_additions = compare_reviewed_product_additions(
        generated,
        canonical,
        addition_policy,
        paths["reviewed_new_product_additions.csv"],
    )
    if not _strict_json_equal(adjudication["delta_comparison"], fresh_delta):
        raise ValueError("Promoted delta adjudication differs from strict recomputation")
    if not _strict_json_equal(adjudication["addition_comparison"], fresh_additions):
        raise ValueError(
            "Promoted addition adjudication differs from strict recomputation"
        )

    return {
        "schema_version": 1,
        "status": "exact_promoted_review_evidence",
        "run_id": run_id,
        "review_freeze_sha256": freeze_hash,
        "reviewed_delta_rows": len(delta_rows),
        "reviewed_addition_rows": len(addition_rows),
        "generated_rows": counts["generated_rows"],
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--analysis-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="canonical analysis/multi_freq_from_archive directory",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="validate only (the checker never writes, with or without this flag)",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        result = validate_promoted_review_evidence(args.analysis_dir)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
