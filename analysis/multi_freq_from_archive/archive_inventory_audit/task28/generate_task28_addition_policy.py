#!/usr/bin/env python3
"""Merge the 82 validated Task 28 ALMA products into the admission policy."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = (
    "schema_version",
    "policy_version",
    "admission_id",
    "source",
    "telescope",
    "expected_frequency_ghz",
    "expected_date_label",
    "product",
    "product_sha256",
    "expected_run_disposition",
    "admission_decision",
    "reviewer",
    "review_date",
    "rationale",
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--existing-policy", type=Path, required=True)
    parser.add_argument("--validated-products", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    existing = read_csv(arguments.existing_policy)
    validated = read_csv(arguments.validated_products)
    if len(existing) != 7 or len(validated) != 82:
        raise RuntimeError(
            f"unexpected policy inventory: existing={len(existing)}, task28={len(validated)}"
        )
    if any(row["validation_status"] != "validated" for row in validated):
        raise RuntimeError("Task 28 policy cannot include an unvalidated product")
    products = {row["product"] for row in existing}
    additions = []
    for index, row in enumerate(sorted(validated, key=lambda item: item["output_relative_path"]), start=1):
        product = row["output_relative_path"]
        if product in products:
            raise RuntimeError(f"Task 28 product already exists in policy: {product}")
        products.add(product)
        rationale_parts = [
            "task28_filename_neutral_source_mous_continuum_admission",
            row["screening_disposition"],
        ]
        if row["transformation"].startswith("extract_"):
            rationale_parts.append("target_stokes_i_extracted_from_hash_recorded_iquv_parent")
        if ".split.cal." in product:
            rationale_parts.append("anonymous_source_index_confirmed_by_fits_object_and_wcs")
        additions.append(
            {
                "schema_version": "1",
                "policy_version": "continuum_product_addition_v1",
                "admission_id": f"CPA-T28-{index:03d}",
                "source": row["source"],
                "telescope": "ALMA",
                "expected_frequency_ghz": f"{float(row['header_frequency_ghz']):.2f}",
                "expected_date_label": row["header_date_obs"][:10].replace("-", "/") + "_1",
                "product": product,
                "product_sha256": row["output_sha256"],
                "expected_run_disposition": "publication_row_required",
                "admission_decision": "admit_reviewed_public_archive_continuum_product",
                "reviewer": "author_approved_task28_archive_audit",
                "review_date": "2026-08-25",
                "rationale": "+".join(rationale_parts),
            }
        )
    merged = sorted(existing + additions, key=lambda row: row["product"])
    if len(merged) != 89 or len({row["product"] for row in merged}) != 89:
        raise RuntimeError("combined addition policy is not exactly 89 unique products")
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    with arguments.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(merged)
    print(f"wrote {len(merged)} reviewed additions to {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
