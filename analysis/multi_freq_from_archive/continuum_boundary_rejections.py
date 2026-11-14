#!/usr/bin/env python3
"""Fail-closed policy for converged fits rejected by the plot-boundary gate.

These products are discovery inputs, but they have never been publication
rows.  The fitter nevertheless records their final local-noise diagnostics
before the historical half-maximum model-contour boundary check rejects them.
This module keeps that diagnostic evidence without allowing arbitrary audit
rows to drift outside the publication table.
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path, PurePosixPath
from typing import Mapping

from continuum_product_additions import fits_header_output_identity, sha256_file


BOUNDARY_REJECTION_POLICY_FIELDS = (
    "schema_version",
    "policy_version",
    "rejection_id",
    "source",
    "telescope",
    "expected_frequency_ghz",
    "expected_observation_date",
    "product",
    "product_sha256",
    "size_bytes",
    "independent_evidence_record",
    "independent_evidence_sha256",
    "historical_boundary_highest_percent",
    "historical_model_payload_sha256",
    "expected_row_disposition",
    "expected_classification",
    "rejection_gate",
    "rejection_decision",
    "reviewer",
    "review_date",
    "historical_evidence",
    "rationale",
)
BOUNDARY_REJECTION_POLICY_VERSION = "continuum_boundary_rejection_v1"
EXPECTED_ROW_DISPOSITION = "rejected_on_boundary"
EXPECTED_CLASSIFICATION = "detection"
REJECTION_GATE = "half_maximum_model_contour_intersects_outer_ten_percent"
REJECTION_DECISION = "retain_noise_diagnostic_exclude_publication_row"
EXPECTED_BOUNDARY_REJECTIONS = 3

BOUNDARY_REJECTION_EVIDENCE_FIELDS = (
    "schema_version",
    "evidence_version",
    "rejection_id",
    "product",
    "independent_source_file",
    "independent_source_sha256",
    "independent_source_record_index",
    "fit_converged",
    "boundary_result",
    "rerun_outcome",
    "historical_boundary_diagnostic",
    "boundary_highest_percent",
    "archived_final_model_exists",
    "archived_final_model_flux_mjy",
    "rerun_model_flux_mjy",
    "model_flux_delta_mjy",
    "retry3_model_payload_sha256",
    "authoritative7_model_payload_sha256",
    "model_payload_size_bytes",
    "input_metadata_unchanged",
    "status",
)
BOUNDARY_REJECTION_EVIDENCE_VERSION = "continuum_boundary_rejection_evidence_v1"
BOUNDARY_REJECTION_EVIDENCE_NAME = "continuum_boundary_rejection_evidence.csv"

DATE_PATTERN = re.compile(r"\d{4}-\d{2}-\d{2}")
OBSERVATION_DATE_PATTERN = re.compile(r"\d{4}/\d{2}/\d{2}")
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")


def _portable_vla_product(value: object) -> str:
    text = str(value).strip().replace("\\", "/")
    path = PurePosixPath(text)
    windows_absolute = len(text) >= 3 and text[1:3] == ":/"
    if (
        not text
        or path.is_absolute()
        or windows_absolute
        or ".." in path.parts
        or path.as_posix() != text
        or len(path.parts) < 3
        or path.parts[0] != "VLA"
        or not text.casefold().endswith((".fits", ".imfits"))
    ):
        raise ValueError(
            f"boundary-rejection product must be a normalized VLA FITS path: {value!r}"
        )
    return text


def load_boundary_rejection_evidence(path: Path) -> list[dict[str, str]]:
    """Validate the shipped extraction of the independent full-fit rerun."""

    if path.is_symlink() or not path.is_file():
        raise ValueError(f"boundary-rejection evidence must be a real file: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        actual_fields = tuple(reader.fieldnames or ())
        if actual_fields != BOUNDARY_REJECTION_EVIDENCE_FIELDS:
            raise ValueError(
                f"{path}: boundary-evidence fields differ from exact schema"
            )
        rows = [dict(row) for row in reader]
    if len(rows) != EXPECTED_BOUNDARY_REJECTIONS:
        raise ValueError(
            f"{path}: expected {EXPECTED_BOUNDARY_REJECTIONS} boundary-evidence rows"
        )
    normalized: list[dict[str, str]] = []
    products: set[str] = set()
    rejection_ids: set[str] = set()
    source_hashes: set[str] = set()
    for line_number, raw in enumerate(rows, start=2):
        if any(value is None or not str(value).strip() for value in raw.values()):
            raise ValueError(f"{path}:{line_number}: evidence fields must be nonblank")
        row = {
            field: str(raw[field]).strip()
            for field in BOUNDARY_REJECTION_EVIDENCE_FIELDS
        }
        product = _portable_vla_product(row["product"])
        rejection_id = row["rejection_id"]
        if (
            row["schema_version"] != "1"
            or row["evidence_version"] != BOUNDARY_REJECTION_EVIDENCE_VERSION
        ):
            raise ValueError(f"{path}:{line_number}: stale evidence version")
        if rejection_id in rejection_ids or product in products:
            raise ValueError(f"{path}:{line_number}: duplicate evidence ID or product")
        rejection_ids.add(rejection_id)
        products.add(product)
        row["product"] = product
        exact_values = {
            "fit_converged": "true",
            "boundary_result": "true",
            "rerun_outcome": "on_boundary",
            "historical_boundary_diagnostic": "true",
            "archived_final_model_exists": "true",
            "model_flux_delta_mjy": "0.0",
            "input_metadata_unchanged": "true",
            "status": "ok",
        }
        for field, expected in exact_values.items():
            if row[field] != expected:
                raise ValueError(
                    f"{path}:{line_number}: evidence {field} must equal {expected!r}"
                )
        if not SHA256_PATTERN.fullmatch(row["independent_source_sha256"]):
            raise ValueError(f"{path}:{line_number}: invalid source evidence SHA-256")
        source_hashes.add(row["independent_source_sha256"])
        for field in (
            "retry3_model_payload_sha256",
            "authoritative7_model_payload_sha256",
        ):
            if not SHA256_PATTERN.fullmatch(row[field]):
                raise ValueError(f"{path}:{line_number}: invalid {field}")
        if (
            row["retry3_model_payload_sha256"]
            != row["authoritative7_model_payload_sha256"]
        ):
            raise ValueError(
                f"{path}:{line_number}: retry3/A7 model payloads differ"
            )
        try:
            record_index = int(row["independent_source_record_index"])
            model_size = int(row["model_payload_size_bytes"])
            highest_percent = float(row["boundary_highest_percent"])
            archived_flux = float(row["archived_final_model_flux_mjy"])
            rerun_flux = float(row["rerun_model_flux_mjy"])
        except ValueError as exc:
            raise ValueError(f"{path}:{line_number}: invalid numeric evidence") from exc
        if (
            record_index <= 0
            or row["independent_source_record_index"] != str(record_index)
            or model_size <= 0
            or row["model_payload_size_bytes"] != str(model_size)
            or not 90.0 <= highest_percent <= 100.0
            or not math.isfinite(archived_flux)
            or not math.isfinite(rerun_flux)
            or archived_flux != rerun_flux
        ):
            raise ValueError(f"{path}:{line_number}: stale boundary evidence values")
        normalized.append(row)
    if len(source_hashes) != 1:
        raise ValueError(f"{path}: evidence rows do not share one source-file SHA-256")
    if [row["product"] for row in normalized] != sorted(products):
        raise ValueError(f"{path}: boundary-evidence rows must be product-sorted")
    return normalized


def load_boundary_rejection_policy(
    path: Path,
    *,
    analysis_root: Path | None = None,
) -> list[dict[str, str]]:
    """Read and validate the exact reviewed boundary-rejection inventory."""

    if path.is_symlink() or not path.is_file():
        raise ValueError(f"boundary-rejection policy must be a real file: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        actual_fields = tuple(reader.fieldnames or ())
        if actual_fields != BOUNDARY_REJECTION_POLICY_FIELDS:
            raise ValueError(
                f"{path}: boundary-rejection fields differ from exact schema; "
                f"expected={list(BOUNDARY_REJECTION_POLICY_FIELDS)}, "
                f"actual={list(actual_fields)}"
            )
        raw_rows = [dict(row) for row in reader]

    if len(raw_rows) != EXPECTED_BOUNDARY_REJECTIONS:
        raise ValueError(
            f"{path}: expected {EXPECTED_BOUNDARY_REJECTIONS} reviewed boundary "
            f"rejections, found {len(raw_rows)}"
        )

    resolved_root = analysis_root.resolve() if analysis_root is not None else None
    products: set[str] = set()
    rejection_ids: set[str] = set()
    normalized: list[dict[str, str]] = []
    for line_number, raw in enumerate(raw_rows, start=2):
        if any(value is None or not str(value).strip() for value in raw.values()):
            raise ValueError(
                f"{path}:{line_number}: boundary-rejection fields must be nonblank"
            )
        row = {
            field: str(raw[field]).strip()
            for field in BOUNDARY_REJECTION_POLICY_FIELDS
        }
        product = _portable_vla_product(row["product"])
        product_path = PurePosixPath(product)
        if (
            row["schema_version"] != "1"
            or row["policy_version"] != BOUNDARY_REJECTION_POLICY_VERSION
        ):
            raise ValueError(
                f"{path}:{line_number}: stale boundary-rejection policy version"
            )
        if row["telescope"] != "VLA" or product_path.parts[1] != row["source"]:
            raise ValueError(
                f"{path}:{line_number}: source/telescope differs from product path"
            )
        try:
            frequency = float(row["expected_frequency_ghz"])
        except ValueError as exc:
            raise ValueError(
                f"{path}:{line_number}: invalid expected_frequency_ghz"
            ) from exc
        if (
            not math.isfinite(frequency)
            or frequency <= 0.0
            or row["expected_frequency_ghz"] != f"{frequency:.2f}"
        ):
            raise ValueError(
                f"{path}:{line_number}: expected_frequency_ghz must use exact "
                "two-decimal output format"
            )
        if not OBSERVATION_DATE_PATTERN.fullmatch(
            row["expected_observation_date"]
        ):
            raise ValueError(f"{path}:{line_number}: invalid observation date")
        if not SHA256_PATTERN.fullmatch(row["product_sha256"]):
            raise ValueError(f"{path}:{line_number}: invalid product SHA-256")
        if not SHA256_PATTERN.fullmatch(row["independent_evidence_sha256"]):
            raise ValueError(f"{path}:{line_number}: invalid evidence SHA-256")
        if not SHA256_PATTERN.fullmatch(row["historical_model_payload_sha256"]):
            raise ValueError(f"{path}:{line_number}: invalid model-payload SHA-256")
        expected_evidence_record = (
            f"{BOUNDARY_REJECTION_EVIDENCE_NAME}#{row['rejection_id']}"
        )
        if row["independent_evidence_record"] != expected_evidence_record:
            raise ValueError(
                f"{path}:{line_number}: independent evidence locator is stale"
            )
        try:
            boundary_highest_percent = float(
                row["historical_boundary_highest_percent"]
            )
        except ValueError as exc:
            raise ValueError(
                f"{path}:{line_number}: invalid historical boundary percentage"
            ) from exc
        if not 90.0 <= boundary_highest_percent <= 100.0:
            raise ValueError(
                f"{path}:{line_number}: historical boundary percentage is outside "
                "the rejected outer-ten-percent region"
            )
        try:
            size = int(row["size_bytes"])
        except ValueError as exc:
            raise ValueError(f"{path}:{line_number}: invalid size_bytes") from exc
        if size <= 0 or row["size_bytes"] != str(size):
            raise ValueError(
                f"{path}:{line_number}: size_bytes must be a canonical positive integer"
            )
        exact_values = {
            "expected_row_disposition": EXPECTED_ROW_DISPOSITION,
            "expected_classification": EXPECTED_CLASSIFICATION,
            "rejection_gate": REJECTION_GATE,
            "rejection_decision": REJECTION_DECISION,
        }
        for field, expected in exact_values.items():
            if row[field] != expected:
                raise ValueError(
                    f"{path}:{line_number}: {field} must equal {expected!r}"
                )
        if not DATE_PATTERN.fullmatch(row["review_date"]):
            raise ValueError(f"{path}:{line_number}: invalid review_date")
        rejection_id = row["rejection_id"]
        if rejection_id in rejection_ids or product in products:
            raise ValueError(
                f"{path}:{line_number}: duplicate rejection ID or product"
            )
        rejection_ids.add(rejection_id)
        products.add(product)
        row["product"] = product

        if resolved_root is not None:
            candidate = resolved_root / product
            resolved = candidate.resolve()
            try:
                resolved.relative_to(resolved_root)
            except ValueError as exc:
                raise ValueError(
                    f"boundary-rejection product escapes analysis root: {product}"
                ) from exc
            if candidate.is_symlink() or not resolved.is_file():
                raise ValueError(
                    f"boundary-rejection product is not a real file: {candidate}"
                )
            if resolved.stat().st_size != size:
                raise ValueError(
                    f"boundary-rejection product size differs for {product}"
                )
            observed_sha256 = sha256_file(resolved)
            if observed_sha256 != row["product_sha256"]:
                raise ValueError(
                    f"boundary-rejection product SHA-256 differs for {product}: "
                    f"expected {row['product_sha256']}, found {observed_sha256}"
                )
            header_frequency, header_date = fits_header_output_identity(resolved)
            if (
                header_frequency != row["expected_frequency_ghz"]
                or header_date != row["expected_observation_date"]
            ):
                raise ValueError(
                    f"boundary-rejection FITS header identity differs for {product}"
                )
        normalized.append(row)

    if [row["product"] for row in normalized] != sorted(products):
        raise ValueError(f"{path}: boundary-rejection rows must be product-sorted")
    # The immutable FITS bytes may live in a separate prepared archive, while
    # the compact evidence extraction always travels beside policy.
    evidence_path = path.parent / BOUNDARY_REJECTION_EVIDENCE_NAME
    evidence_sha256 = sha256_file(evidence_path)
    if any(
        row["independent_evidence_sha256"] != evidence_sha256
        for row in normalized
    ):
        raise ValueError(
            "boundary-rejection policy does not hash-bind the shipped evidence"
        )
    evidence_rows = load_boundary_rejection_evidence(evidence_path)
    evidence_by_id = {row["rejection_id"]: row for row in evidence_rows}
    if set(evidence_by_id) != rejection_ids:
        raise ValueError(
            "boundary-rejection policy/evidence rejection IDs differ"
        )
    for row in normalized:
        evidence = evidence_by_id[row["rejection_id"]]
        if (
            evidence["product"] != row["product"]
            or evidence["boundary_highest_percent"]
            != row["historical_boundary_highest_percent"]
            or evidence["retry3_model_payload_sha256"]
            != row["historical_model_payload_sha256"]
        ):
            raise ValueError(
                f"boundary policy/evidence values differ for {row['product']}"
            )
    return normalized


def boundary_rejection_products(rows: list[Mapping[str, str]]) -> set[str]:
    """Return the exact normalized product inventory from validated rows."""

    return {_portable_vla_product(row["product"]) for row in rows}


def validate_boundary_rejection_record_identity(
    record: Mapping[str, str],
    policy: Mapping[str, str],
) -> None:
    """Require one fresh audit record to match its exact reviewed policy row."""

    product = _portable_vla_product(record["product"])
    if product != _portable_vla_product(policy["product"]):
        raise RuntimeError("boundary-rejection record/product policy mismatch")
    try:
        frequency = float(record["frequency_ghz"])
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeError(
            f"boundary-rejection audit frequency is invalid for {product}"
        ) from exc
    if (
        not math.isfinite(frequency)
        or f"{frequency:.2f}" != policy["expected_frequency_ghz"]
        or record["source"] != policy["source"]
        or record["telescope"] != policy["telescope"]
        or record["row_disposition"] != policy["expected_row_disposition"]
        or record["classification"] != policy["expected_classification"]
    ):
        raise RuntimeError(
            f"boundary-rejection audit identity/disposition is stale for {product}"
        )
