"""Fail-closed identity review for scaled-MAD continuum classifications."""

from __future__ import annotations

import ast
import csv
import hashlib
import math
import os
import re
from pathlib import Path
from typing import Any

from continuum_measurements import portable_release_product
from continuum_noise_inventory import GENERIC_FINAL_NOISE_AUDIT_FIELDS
from continuum_product_additions import (
    addition_products,
    fitter_output_source,
    load_product_addition_policy,
    publication_addition_products,
    publication_addition_products_in_fitter_order,
)
from continuum_boundary_rejections import (
    boundary_rejection_products,
    load_boundary_rejection_policy,
    validate_boundary_rejection_record_identity,
)
from rerun_ngc4945_nuclear import EXPECTED_ROWS as EXPECTED_NGC4945_ROWS
from rerun_ngc4945_nuclear import PHOTOMETRY_FIELDS as NGC4945_PHOTOMETRY_FIELDS


POLICY_FIELDS = (
    "schema_version",
    "policy_version",
    "review_id",
    "trigger",
    "source",
    "telescope",
    "frequency_ghz",
    "product",
    "product_sha256",
    "reference_role",
    "reference_frame",
    "reference_ra_deg",
    "reference_dec_deg",
    "reference_evidence",
    "identity_evidence_class",
    "component_identity_status",
    "component_label",
    "measurement_semantics",
    "background_morphology",
    "expected_scaled_mad_classification",
    "review_decision",
    "reviewer",
    "review_date",
    "rationale",
)
POLICY_VERSION = "scaled_mad_identity_v1"
EXPECTED_CANONICAL_ROWS = 409
EXPECTED_GENERATED_ROWS = 406
EXPECTED_EXCLUSIONS = 3
DATE_PATTERN = re.compile(r"\d{4}-\d{2}-\d{2}")

IDENTITY_CLASSES = {
    "a_priori_position_match",
    "published_same_observation_component",
    "published_same_observation_context",
    "target_associated_blend",
}
IDENTITY_STATUSES = {
    "confirmed_target_component",
    "confirmed_target_associated_blend",
}
MEASUREMENT_SEMANTICS = {
    "compact_component_flux",
    "target_associated_blended_flux",
}
BACKGROUND_MORPHOLOGIES = {"noise_like", "structured_source_emission"}


def _sexagesimal_ra_deg(value: str) -> float:
    """Convert an hours:minutes:seconds right ascension to degrees."""

    fields = str(value).strip().split(":")
    if len(fields) != 3:
        raise RuntimeError(f"invalid sexagesimal right ascension: {value!r}")
    hours, minutes, seconds = map(float, fields)
    if not (0.0 <= hours < 24.0 and 0.0 <= minutes < 60.0 and 0.0 <= seconds < 60.0):
        raise RuntimeError(f"invalid sexagesimal right ascension: {value!r}")
    return 15.0 * (hours + minutes / 60.0 + seconds / 3600.0)


def _sexagesimal_dec_deg(value: str) -> float:
    """Convert a signed degrees:minutes:seconds declination to degrees."""

    text = str(value).strip()
    sign = -1.0 if text.startswith("-") else 1.0
    fields = text.lstrip("+-").split(":")
    if len(fields) != 3:
        raise RuntimeError(f"invalid sexagesimal declination: {value!r}")
    degrees, minutes, seconds = map(float, fields)
    if not (0.0 <= degrees <= 90.0 and 0.0 <= minutes < 60.0 and 0.0 <= seconds < 60.0):
        raise RuntimeError(f"invalid sexagesimal declination: {value!r}")
    result = sign * (degrees + minutes / 60.0 + seconds / 3600.0)
    if not -90.0 <= result <= 90.0:
        raise RuntimeError(f"invalid sexagesimal declination: {value!r}")
    return result


def load_ngc4945_source18_identity_records(path: Path) -> dict[str, dict[str, str]]:
    """Return source-18 photometry in the fields used by identity review.

    The generic final-noise audit deliberately records the rejected source-17
    fits for NGC 4945.  Publication rows are replaced later by the fixed
    source-18 refit, so component review must use the latter coordinates,
    classifications, and aperture geometry rather than the generic seed.
    """

    rows = _read_exact_csv(
        path,
        NGC4945_PHOTOMETRY_FIELDS,
        "NGC 4945 source-18 photometry",
    )
    if len(rows) != EXPECTED_NGC4945_ROWS:
        raise RuntimeError(
            f"expected {EXPECTED_NGC4945_ROWS} NGC 4945 source-18 rows; "
            f"found {len(rows)}"
        )
    records: dict[str, dict[str, str]] = {}
    for row in rows:
        product = portable_release_product(row["product"])
        if product in records:
            raise RuntimeError(f"duplicate NGC 4945 source-18 product: {product}")
        if (
            row["source"] != "NGC4945"
            or row["telescope"] != "ALMA"
            or row["coordinate_frame"] != "ICRS"
            or row["component_identity"] != "Emig2020_source18_AGN_maser_nucleus"
            or row["centroid_policy"] != "fixed_xy"
            or row["row_disposition"] != "retained_source18_replacement"
            or row["noise_guard_status"] != "passed"
            or row["status"] != "detection"
        ):
            raise RuntimeError(f"invalid NGC 4945 source-18 identity row: {product}")
        records[product] = {
            "source": row["source"],
            "telescope": row["telescope"],
            "frequency_ghz": row["frequency_ghz"],
            "product": product,
            "fit_center_reference_frame": row["coordinate_frame"],
            "fit_center_ra_deg": f"{_sexagesimal_ra_deg(row['fixed_ra']):.15f}",
            "fit_center_dec_deg": f"{_sexagesimal_dec_deg(row['fixed_dec']):.15f}",
            "beam_major_arcsec": row["beam_major_arcsec"],
            "beam_minor_arcsec": row["beam_minor_arcsec"],
            "beam_pa_deg": row["beam_pa_deg"],
            "aperture_half_width_x_pixels": row[
                "noise_aperture_half_width_x_pixels"
            ],
            "aperture_half_width_y_pixels": row[
                "noise_aperture_half_width_y_pixels"
            ],
            "pixel_width_arcsec": row["pixel_width_arcsec"],
            "pixel_height_arcsec": row["pixel_height_arcsec"],
            "flux_mjy": row["unrounded_flux_mjy"],
            "snr": row["new_snr"],
            "classification": "detection",
            "measurement_scope": "publication_fixed_source18_source_excluded",
            "row_disposition": row["row_disposition"],
            "guard_status": row["noise_guard_status"],
        }
    return records


def _sha256_path(path: Path) -> str:
    """Hash an immutable FITS file or CASA table using workflow semantics."""

    digest = hashlib.sha256()
    if path.is_file() and not path.is_symlink():
        with path.open("rb") as handle:
            while chunk := handle.read(16 * 1024 * 1024):
                digest.update(chunk)
        return digest.hexdigest()
    if not path.is_dir() or path.is_symlink():
        raise RuntimeError(f"identity-policy product is not an isolated path: {path}")
    for item in sorted(path.rglob("*"), key=lambda value: value.relative_to(path).as_posix()):
        relative = item.relative_to(path).as_posix().encode("utf-8")
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        if item.is_symlink():
            target = os.readlink(item).encode("utf-8")
            digest.update(b"L" + len(target).to_bytes(8, "big") + target)
        elif item.is_file():
            digest.update(b"F")
            with item.open("rb") as handle:
                while chunk := handle.read(16 * 1024 * 1024):
                    digest.update(chunk)
        elif item.is_dir():
            digest.update(b"D")
    return digest.hexdigest()


def _read_sorted(path: Path, expected: int, label: str) -> list[list[Any]]:
    rows = []
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        if not line.strip():
            continue
        try:
            row = ast.literal_eval(line)
        except (SyntaxError, ValueError) as error:
            raise RuntimeError(f"invalid {label} row at {path}:{line_number}") from error
        if not isinstance(row, list) or len(row) != 8:
            raise RuntimeError(f"{label} row is not an exact eight-field list")
        row[7] = portable_release_product(row[7])
        rows.append(row)
    if len(rows) != expected:
        raise RuntimeError(f"expected {expected} {label} rows; found {len(rows)}")
    products = [str(row[7]) for row in rows]
    if "" in products or len(products) != len(set(products)):
        raise RuntimeError(f"{label} products must be nonblank and unique")
    return rows


def _read_exact_csv(
    path: Path,
    fields: tuple[str, ...],
    label: str,
) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        actual = tuple(reader.fieldnames or ())
        if actual != fields:
            raise RuntimeError(f"{label} fields differ from exact schema: {actual!r}")
        rows = [dict(row) for row in reader]
    return rows


def _finite(record: dict[str, str], field: str, *, positive: bool = False) -> float:
    try:
        value = float(record[field])
    except (KeyError, TypeError, ValueError) as error:
        raise RuntimeError(f"invalid {field} for {record.get('product')}") from error
    if not math.isfinite(value) or (positive and value <= 0.0):
        raise RuntimeError(f"invalid {field} for {record.get('product')}")
    return value


def _is_detection(row: list[Any]) -> bool:
    return not str(row[3]).endswith("*") and str(row[6]) != "NA"


def _coordinate(record: dict[str, str], prefix: str) -> tuple[str, float, float]:
    frame_field = (
        "fit_center_reference_frame" if prefix == "fit_center" else f"{prefix}_frame"
    )
    frame = record[frame_field]
    ra = _finite(record, f"{prefix}_ra_deg")
    dec = _finite(record, f"{prefix}_dec_deg")
    if frame not in {"ICRS", "J2000"}:
        raise RuntimeError(f"unsupported {prefix} frame for {record['product']}")
    if not 0.0 <= ra < 360.0 or not -90.0 <= dec <= 90.0:
        raise RuntimeError(f"invalid {prefix} coordinate for {record['product']}")
    return frame, ra, dec


def _separation_arcsec(
    ra1_deg: float,
    dec1_deg: float,
    ra2_deg: float,
    dec2_deg: float,
) -> float:
    """Return stable great-circle separation without an optional dependency."""

    ra1, dec1, ra2, dec2 = map(
        math.radians, (ra1_deg, dec1_deg, ra2_deg, dec2_deg)
    )
    delta_ra = ra2 - ra1
    delta_dec = dec2 - dec1
    haversine = (
        math.sin(delta_dec / 2.0) ** 2
        + math.cos(dec1) * math.cos(dec2) * math.sin(delta_ra / 2.0) ** 2
    )
    angle = 2.0 * math.asin(min(1.0, math.sqrt(max(0.0, haversine))))
    return math.degrees(angle) * 3600.0


def validate_classification_identity_policy(
    *,
    workspace: Path,
    canonical_table: Path,
    generated_table: Path,
    exclusion_policy: Path,
    generic_noise_audit: Path,
    identity_policy: Path,
    addition_policy: Path | None = None,
    boundary_rejection_policy: Path | None = None,
) -> dict[str, object]:
    """Validate the exact fresh-result identity-review inventory."""

    workspace = workspace.resolve()
    analysis_root = workspace / "analysis/multi_freq_from_archive"
    canonical = _read_sorted(
        canonical_table, EXPECTED_CANONICAL_ROWS, "canonical continuum"
    )
    additions = (
        []
        if addition_policy is None
        else load_product_addition_policy(
            addition_policy,
            analysis_root=analysis_root,
        )
    )
    admitted = addition_products(additions)
    publication_additions = publication_addition_products(additions)
    boundary_rows = (
        []
        if boundary_rejection_policy is None
        else load_boundary_rejection_policy(
            boundary_rejection_policy,
            analysis_root=analysis_root,
        )
    )
    boundary_by_product = {row["product"]: row for row in boundary_rows}
    boundary_products = boundary_rejection_products(boundary_rows)
    generated = _read_sorted(
        generated_table,
        EXPECTED_GENERATED_ROWS + len(publication_additions),
        "generated continuum",
    )
    with exclusion_policy.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if "product" not in (reader.fieldnames or ()):
            raise RuntimeError("exclusion policy lacks product")
        excluded = [portable_release_product(row["product"]) for row in reader]
    if len(excluded) != EXPECTED_EXCLUSIONS or len(excluded) != len(set(excluded)):
        raise RuntimeError("classification identity review requires three exclusions")
    excluded_set = set(excluded)
    corrected = [row for row in canonical if str(row[7]) not in excluded_set]
    generated_common = [row for row in generated if str(row[7]) not in admitted]
    generated_additions = [str(row[7]) for row in generated if str(row[7]) in admitted]
    if [str(row[7]) for row in corrected] != [str(row[7]) for row in generated_common]:
        raise RuntimeError("canonical-minus-exclusions/generated identity order differs")
    if generated_additions != publication_addition_products_in_fitter_order(
        additions
    ):
        raise RuntimeError("generated product additions differ from exact policy order")

    audit = _read_exact_csv(
        generic_noise_audit,
        GENERIC_FINAL_NOISE_AUDIT_FIELDS,
        "generic final-noise audit",
    )
    audit_by_product = {
        portable_release_product(row["product"]): row for row in audit
    }
    if len(audit_by_product) != len(audit):
        raise RuntimeError("generic final-noise audit contains duplicate products")
    generated_products = [str(row[7]) for row in generated]
    if set(generated_products) & boundary_products:
        raise RuntimeError(
            "boundary-rejection products must be absent from the generated table"
        )
    expected_audit_products = set(generated_products) | boundary_products
    if len(audit) != len(expected_audit_products) or set(audit_by_product) != expected_audit_products:
        raise RuntimeError("generic final-noise audit inventory differs from generated table")
    for product, policy_row in boundary_by_product.items():
        validate_boundary_rejection_record_identity(
            audit_by_product[product], policy_row
        )

    required: dict[str, dict[str, object]] = {}
    for old_row, new_row in zip(corrected, generated_common):
        product = str(new_row[7])
        record = audit_by_product.get(product)
        if record is None:
            raise RuntimeError(f"missing generic final-noise audit for {product}")
        if (
            record["source"] != str(new_row[0])
            or record["telescope"] != str(new_row[1])
            or not math.isclose(
                _finite(record, "frequency_ghz", positive=True),
                float(new_row[2]),
                rel_tol=0.0,
                abs_tol=5e-3,
            )
        ):
            raise RuntimeError(f"generic final-noise audit identity is stale for {product}")
        if record["source"] == "NGC4945" or record["forced_policy"] != "no":
            continue
        flux = _finite(record, "flux_mjy")
        rms = _finite(record, "rms_mjy_per_beam", positive=True)
        scaled_mad = _finite(record, "scaled_mad_mjy_per_beam", positive=True)
        old_detection = _is_detection(old_row)
        mad_detection = flux / scaled_mad >= 3.0
        rms_detection = flux / rms >= 3.0
        if old_detection and not mad_detection:
            raise RuntimeError(
                "scaled-MAD lost-detection transition has no explicit reviewed "
                f"policy schema: {product}"
            )
        new_detection = (not old_detection) and mad_detection
        crossing = mad_detection != rms_detection
        if not (new_detection or crossing):
            continue
        trigger = "+".join(
            label
            for label, active in (
                ("new_detection", new_detection),
                ("rms_mad_threshold_crossing", crossing),
            )
            if active
        )
        expected_class = "detection" if mad_detection else "ordinary_upper_limit"
        if record["classification"] != expected_class:
            raise RuntimeError(f"scaled-MAD classification is stale for {product}")
        required[product] = {
            "trigger": trigger,
            "expected_class": expected_class,
            "record": record,
        }

    generated_by_product = {str(row[7]): row for row in generated}
    additions_by_product = {row["product"]: row for row in additions}
    identity_record_by_product = dict(audit_by_product)
    ngc4945_additions = {
        product
        for product in publication_additions
        if fitter_output_source(additions_by_product[product]) == "NGC4945"
    }
    if ngc4945_additions:
        source18_records = load_ngc4945_source18_identity_records(
            analysis_root / "ngc4945_nuclear_photometry.csv"
        )
        generated_ngc4945_products = {
            product
            for product, row in generated_by_product.items()
            if str(row[0]) == "NGC4945"
        }
        if set(source18_records) != generated_ngc4945_products:
            raise RuntimeError(
                "NGC 4945 source-18 photometry inventory differs from the "
                "publication table"
            )
        identity_record_by_product.update(source18_records)
    for product in sorted(publication_additions):
        new_row = generated_by_product[product]
        record = identity_record_by_product[product]
        addition = additions_by_product[product]
        if (
            record["source"] != str(new_row[0])
            or record["source"] != fitter_output_source(addition)
            or record["telescope"] != "ALMA"
            or str(new_row[1]) != "ALMA"
            or str(new_row[2]) != addition["expected_frequency_ghz"]
            or str(new_row[4]) != addition["expected_date_label"]
        ):
            raise RuntimeError(f"new-product audit identity is stale for {product}")
        detection = _is_detection(new_row)
        expected_class = "detection" if detection else "ordinary_upper_limit"
        if record["classification"] != expected_class:
            raise RuntimeError(f"new-product classification is stale for {product}")
        # Upper limits have admission provenance in the separate product ledger
        # but no adopted component whose identity could be reviewed.
        if detection:
            required[product] = {
                "trigger": "new_product_admission",
                "expected_class": expected_class,
                "record": record,
            }

    policy = _read_exact_csv(identity_policy, POLICY_FIELDS, "identity policy")
    products = [portable_release_product(row["product"]) for row in policy]
    if products != sorted(products) or len(products) != len(set(products)):
        raise RuntimeError("identity-policy products must be unique and sorted")
    if set(products) != set(required):
        raise RuntimeError(
            "identity-policy inventory differs from fresh required products: "
            f"missing={sorted(set(required) - set(products))}, "
            f"extra={sorted(set(products) - set(required))}"
        )

    validation_rows: list[dict[str, object]] = []
    review_ids: set[str] = set()
    for row, product in zip(policy, products):
        row["product"] = product
        if row["schema_version"] != "1" or row["policy_version"] != POLICY_VERSION:
            raise RuntimeError(f"stale identity-policy version for {product}")
        if not row["review_id"] or row["review_id"] in review_ids:
            raise RuntimeError(f"identity-policy review IDs must be unique for {product}")
        review_ids.add(row["review_id"])
        expected = required[product]
        record = expected["record"]
        assert isinstance(record, dict)
        generated_row = generated_by_product[product]
        if (
            row["trigger"] != expected["trigger"]
            or row["source"] != str(generated_row[0])
            or row["telescope"] != str(generated_row[1])
            or not math.isclose(
                _finite(row, "frequency_ghz", positive=True),
                float(generated_row[2]),
                rel_tol=0.0,
                abs_tol=5e-3,
            )
            or row["expected_scaled_mad_classification"]
            != expected["expected_class"]
            or row["review_decision"] != "accept_scaled_mad_classification"
        ):
            raise RuntimeError(f"identity-policy classification fields are stale for {product}")

        product_path = analysis_root / product
        if row["product_sha256"] != _sha256_path(product_path):
            raise RuntimeError(f"identity-policy product SHA-256 is stale for {product}")
        if row["identity_evidence_class"] not in IDENTITY_CLASSES:
            raise RuntimeError(f"unknown identity evidence class for {product}")
        if row["component_identity_status"] not in IDENTITY_STATUSES:
            raise RuntimeError(f"unresolved component identity for {product}")
        if row["measurement_semantics"] not in MEASUREMENT_SEMANTICS:
            raise RuntimeError(f"unknown measurement semantics for {product}")
        if row["background_morphology"] not in BACKGROUND_MORPHOLOGIES:
            raise RuntimeError(f"unknown background morphology for {product}")
        for field in (
            "reference_role",
            "reference_evidence",
            "component_label",
            "reviewer",
            "rationale",
        ):
            if not row[field].strip():
                raise RuntimeError(f"blank identity-policy {field} for {product}")
        if not DATE_PATTERN.fullmatch(row["review_date"]):
            raise RuntimeError(f"invalid review date for {product}")

        evidence_class = row["identity_evidence_class"]
        if evidence_class in {
            "published_same_observation_component",
            "published_same_observation_context",
        }:
            if "doi:" not in row["reference_evidence"].casefold():
                raise RuntimeError(f"published identity evidence lacks DOI for {product}")
            max_separation_beams = 1.0
        else:
            max_separation_beams = 0.5
        if evidence_class == "target_associated_blend":
            if (
                row["component_identity_status"]
                != "confirmed_target_associated_blend"
                or row["measurement_semantics"] != "target_associated_blended_flux"
            ):
                raise RuntimeError(f"blend identity semantics are inconsistent for {product}")
        elif (
            row["component_identity_status"] != "confirmed_target_component"
            or row["measurement_semantics"] != "compact_component_flux"
        ):
            raise RuntimeError(f"compact identity semantics are inconsistent for {product}")

        fit_frame, fit_ra, fit_dec = _coordinate(record, "fit_center")
        reference_frame, reference_ra, reference_dec = _coordinate(row, "reference")
        if fit_frame != reference_frame:
            raise RuntimeError(
                f"identity-policy reference frame differs from fitted native frame "
                f"for {product}: {reference_frame} vs {fit_frame}"
            )
        separation_arcsec = _separation_arcsec(
            fit_ra,
            fit_dec,
            reference_ra,
            reference_dec,
        )
        beam_major = _finite(record, "beam_major_arcsec", positive=True)
        separation_beams = separation_arcsec / beam_major
        if separation_beams > max_separation_beams:
            raise RuntimeError(
                f"identity-policy centroid separation exceeds {max_separation_beams} "
                f"BMAJ for {product}: {separation_beams}"
            )
        validation_rows.append(
            {
                "review_id": row["review_id"],
                "product": product,
                "trigger": row["trigger"],
                "scaled_mad_classification": expected["expected_class"],
                "reference_frame": fit_frame,
                "fit_ra_deg": fit_ra,
                "fit_dec_deg": fit_dec,
                "reference_ra_deg": reference_ra,
                "reference_dec_deg": reference_dec,
                "separation_arcsec": separation_arcsec,
                "separation_bmaj": separation_beams,
                "maximum_separation_bmaj": max_separation_beams,
                "identity_evidence_class": evidence_class,
                "component_identity_status": row["component_identity_status"],
            }
        )

    return {
        "schema_version": 1,
        "status": "exact_required_identity_reviews_resolved",
        "required_review_rows": len(required),
        "new_detection_rows": sum(
            "new_detection" in str(item["trigger"]) for item in required.values()
        ),
        "rms_mad_threshold_crossing_rows": sum(
            "rms_mad_threshold_crossing" in str(item["trigger"])
            for item in required.values()
        ),
        "new_product_admission_rows": sum(
            item["trigger"] == "new_product_admission"
            for item in required.values()
        ),
        "admitted_product_upper_limit_rows": sum(
            not _is_detection(generated_by_product[product])
            for product in publication_additions
        ),
        "validated_reviews": validation_rows,
    }
