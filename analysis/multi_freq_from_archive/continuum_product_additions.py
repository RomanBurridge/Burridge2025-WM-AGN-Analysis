#!/usr/bin/env python3
"""Fail-closed admission policy for continuum products absent from the baseline.

The 409-row historical continuum table remains immutable.  A newly recovered
public ALMA continuum product can enter a fresh fit inventory only through the
exact policy schema below, with a release-relative path and byte-level SHA-256.
An empty policy (header only) is valid and represents zero additions.
"""

from __future__ import annotations

import csv
import hashlib
import math
import re
from pathlib import Path, PurePosixPath
from typing import Mapping

from continuum_measurements import calendar_date_sort_key


ADDITION_POLICY_FIELDS = (
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
ADDITION_POLICY_VERSION = "continuum_product_addition_v1"
ADMISSION_DECISION = "admit_reviewed_public_archive_continuum_product"
EXPECTED_RUN_DISPOSITIONS = {
    "publication_row_required",
    "input_only_expected_no_publication_row",
}
DATE_PATTERN = re.compile(r"\d{4}-\d{2}-\d{2}")
DATE_LABEL_PATTERN = re.compile(r"\d{4}/\d{2}/\d{2}_[1-9]\d*")
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _fits_card_value(card: bytes) -> object:
    text = card[10:80].decode("ascii", errors="strict").rstrip()
    if text.startswith("'"):
        value = text[1:]
        end = value.find("'")
        if end < 0:
            raise ValueError("unterminated FITS string value")
        return value[:end].replace("''", "'").strip()
    token = text.split("/", maxsplit=1)[0].strip()
    if not token:
        raise ValueError("blank FITS value")
    try:
        return float(token.replace("D", "E"))
    except ValueError:
        return token


def fits_header_output_identity(path: Path) -> tuple[str, str]:
    """Read DATE-OBS and the frequency-axis CRVAL from a primary FITS header."""

    header: dict[str, object] = {}
    with path.open("rb") as handle:
        for _block_index in range(128):
            block = handle.read(2880)
            if len(block) != 2880:
                raise ValueError(f"truncated FITS primary header: {path}")
            for offset in range(0, 2880, 80):
                card = block[offset : offset + 80]
                keyword = card[:8].decode("ascii", errors="strict").strip()
                if keyword == "END":
                    block = b""
                    break
                if keyword and card[8:10] == b"= ":
                    header[keyword] = _fits_card_value(card)
            if not block:
                break
        else:
            raise ValueError(f"FITS primary header has no END card: {path}")

    date_obs = str(header.get("DATE-OBS", ""))
    if not re.match(r"^\d{4}-\d{2}-\d{2}(?:T|$)", date_obs):
        raise ValueError(f"FITS header lacks a parseable DATE-OBS: {path}")
    frequency_axis = None
    for key, value in header.items():
        match = re.fullmatch(r"CTYPE(\d+)", key)
        if match and str(value).strip().upper().startswith("FREQ"):
            frequency_axis = match.group(1)
            break
    if frequency_axis is None:
        raise ValueError(f"FITS header has no frequency axis: {path}")
    try:
        frequency = float(header[f"CRVAL{frequency_axis}"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"FITS frequency axis lacks CRVAL: {path}") from exc
    unit = str(header.get(f"CUNIT{frequency_axis}", "Hz")).strip().casefold()
    factors = {"hz": 1e-9, "khz": 1e-6, "mhz": 1e-3, "ghz": 1.0}
    if unit not in factors or not math.isfinite(frequency) or frequency <= 0.0:
        raise ValueError(f"unsupported or invalid FITS frequency unit/value: {path}")
    return f"{frequency * factors[unit]:.2f}", date_obs[:10].replace("-", "/")


def portable_product(value: object) -> str:
    """Return one normalized release-relative ALMA FITS product path."""

    text = str(value).strip().replace("\\", "/")
    path = PurePosixPath(text)
    windows_absolute = len(text) >= 3 and text[1:3] == ":/"
    if (
        not text
        or path.is_absolute()
        or windows_absolute
        or ".." in path.parts
        or path.as_posix() != text
    ):
        raise ValueError(f"addition product must be a normalized relative path: {value!r}")
    if len(path.parts) < 3 or path.parts[0] != "ALMA":
        raise ValueError(f"addition product must be beneath ALMA/<source>/: {text!r}")
    if not text.casefold().endswith((".fits", ".imfits")):
        raise ValueError(f"addition product must be a FITS file: {text!r}")
    return text


def load_product_addition_policy(
    path: Path,
    *,
    analysis_root: Path | None = None,
) -> list[dict[str, str]]:
    """Read and strictly validate the reviewed product-admission policy.

    When ``analysis_root`` is supplied, every admitted product must be a real
    non-symlink file below that root and must match its recorded SHA-256.
    """

    if path.is_symlink() or not path.is_file():
        raise ValueError(f"addition policy must be a real file: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        actual_fields = tuple(reader.fieldnames or ())
        if actual_fields != ADDITION_POLICY_FIELDS:
            raise ValueError(
                f"{path}: addition-policy fields differ from exact schema; "
                f"expected={list(ADDITION_POLICY_FIELDS)}, actual={list(actual_fields)}"
            )
        rows = [dict(row) for row in reader]

    products: set[str] = set()
    admission_ids: set[str] = set()
    normalized: list[dict[str, str]] = []
    resolved_root = analysis_root.resolve() if analysis_root is not None else None
    for line_number, raw in enumerate(rows, start=2):
        if any(value is None or not str(value).strip() for value in raw.values()):
            raise ValueError(f"{path}:{line_number}: addition-policy fields must be nonblank")
        row = {field: str(raw[field]).strip() for field in ADDITION_POLICY_FIELDS}
        product = portable_product(row["product"])
        source = row["source"]
        product_path = PurePosixPath(product)
        if row["schema_version"] != "1" or row["policy_version"] != ADDITION_POLICY_VERSION:
            raise ValueError(f"{path}:{line_number}: stale addition-policy version")
        if row["telescope"] != "ALMA":
            raise ValueError(f"{path}:{line_number}: only reviewed ALMA additions are allowed")
        try:
            frequency = float(row["expected_frequency_ghz"])
        except ValueError as exc:
            raise ValueError(f"{path}:{line_number}: invalid expected_frequency_ghz") from exc
        if (
            not math.isfinite(frequency)
            or frequency <= 0.0
            or row["expected_frequency_ghz"] != f"{frequency:.2f}"
        ):
            raise ValueError(
                f"{path}:{line_number}: expected_frequency_ghz must use exact two-decimal output format"
            )
        if not DATE_LABEL_PATTERN.fullmatch(row["expected_date_label"]):
            raise ValueError(f"{path}:{line_number}: invalid expected_date_label")
        if product_path.parts[1] != source:
            raise ValueError(
                f"{path}:{line_number}: source does not match the product directory"
            )
        if row["admission_decision"] != ADMISSION_DECISION:
            raise ValueError(f"{path}:{line_number}: product is not explicitly admitted")
        if row["expected_run_disposition"] not in EXPECTED_RUN_DISPOSITIONS:
            raise ValueError(f"{path}:{line_number}: invalid expected_run_disposition")
        if not DATE_PATTERN.fullmatch(row["review_date"]):
            raise ValueError(f"{path}:{line_number}: invalid review_date")
        if not SHA256_PATTERN.fullmatch(row["product_sha256"]):
            raise ValueError(f"{path}:{line_number}: invalid product_sha256")
        admission_id = row["admission_id"]
        if admission_id in admission_ids or product in products:
            raise ValueError(f"{path}:{line_number}: duplicate admission ID or product")
        admission_ids.add(admission_id)
        products.add(product)
        row["product"] = product

        if resolved_root is not None:
            candidate = resolved_root / product
            resolved = candidate.resolve()
            try:
                resolved.relative_to(resolved_root)
            except ValueError as exc:
                raise ValueError(f"admitted product escapes analysis root: {product}") from exc
            if candidate.is_symlink() or not resolved.is_file():
                raise ValueError(f"admitted product is not a real file: {candidate}")
            observed = sha256_file(resolved)
            if observed != row["product_sha256"]:
                raise ValueError(
                    f"admitted product SHA-256 mismatch for {product}: "
                    f"expected {row['product_sha256']}, found {observed}"
                )
            header_frequency, header_date = fits_header_output_identity(resolved)
            if header_frequency != row["expected_frequency_ghz"]:
                raise ValueError(
                    f"admitted FITS header frequency differs from policy for {product}: "
                    f"expected {row['expected_frequency_ghz']}, found {header_frequency}"
                )
            expected_date = row["expected_date_label"].split("_", maxsplit=1)[0]
            if header_date != expected_date:
                raise ValueError(
                    f"admitted FITS DATE-OBS differs from policy for {product}: "
                    f"expected {expected_date}, found {header_date}"
                )
        normalized.append(row)

    if [row["product"] for row in normalized] != sorted(products):
        raise ValueError(f"{path}: addition-policy rows must be product-sorted")
    return normalized


def addition_products(rows: list[Mapping[str, str]]) -> set[str]:
    return {portable_product(row["product"]) for row in rows}


def publication_addition_products(rows: list[Mapping[str, str]]) -> set[str]:
    return {
        portable_product(row["product"])
        for row in rows
        if row["expected_run_disposition"] == "publication_row_required"
    }


def fitter_output_source(row: Mapping[str, str]) -> str:
    """Return the source label emitted by the historical continuum fitter.

    The on-disk ALMA directory is part of the product identity and therefore
    remains ``CircinusGalaxy`` in the admission policy.  The fitter removes
    the literal ``Galaxy`` token while constructing its publication source
    list, so its corresponding table/audit label is ``Circinus``.  Keep that
    legacy normalization explicit here instead of weakening either identity.
    """

    source = str(row["source"]).replace("Galaxy", "")
    if not source:
        raise ValueError("addition source normalizes to a blank fitter label")
    return source


def publication_addition_products_in_fitter_order(
    rows: list[Mapping[str, str]],
) -> list[str]:
    """Return publication additions in the historical fitter's row order.

    Admission-policy rows are intentionally stored in product-path order, but
    the fitter writes its table by normalized source, numeric frequency, and
    calendar date before using the product as a final deterministic tie-break.
    The admission gate forbids same-source/frequency/date collisions, so a
    beam-size key is neither available in policy nor needed to determine the
    relative order of admitted rows.
    """

    selected = [
        row
        for row in rows
        if row["expected_run_disposition"] == "publication_row_required"
    ]
    ordered = sorted(
        selected,
        key=lambda row: (
            fitter_output_source(row),
            float(row["expected_frequency_ghz"]),
            calendar_date_sort_key(row["expected_date_label"]),
            portable_product(row["product"]),
        ),
    )
    return [portable_product(row["product"]) for row in ordered]
