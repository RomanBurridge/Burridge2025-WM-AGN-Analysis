#!/usr/bin/env python3

"""Build the continuum summary tables from the combined fit summary.

This script replaces the fragile notebook-only table assembly logic.  It keeps
the literature extras in one place, applies the audited archive-product alias
and exclusion policies, verifies that rejected products are absent from a
fresh fitter summary, verifies the 21-source analysis sample and the separate
NGC 5495 catalog entry, and generates both the repo-local fixed-width table
and an AAS/CDS MRT version.
"""

from __future__ import annotations

import argparse
import ast
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from continuum_measurements import portable_release_product
from vla_bandwidth_smearing import (
    apply_policy as apply_bandwidth_smearing_policy,
    load_policy as load_bandwidth_smearing_policy,
    write_report as write_bandwidth_smearing_report,
)


@dataclass(frozen=True)
class ExtraPoint:
    source: str
    telescope: str
    freq_ghz: str
    flux_mjy: str
    flux_err_mjy: str
    beam_arcsec: str
    date: str
    reference: str


@dataclass(frozen=True)
class ProductAliasPolicy:
    alias_group_id: str
    source: str
    role: str
    product: str
    sha256: str
    size_bytes: str
    header_observer: str
    header_date_obs: str
    header_date_map: str
    header_frequency_hz: str
    header_naxis1: str
    header_naxis2: str
    header_bmaj_deg: str
    header_bmin_deg: str
    header_clean_niter: str
    pixel_comparison_evidence: str
    decision_reason: str


@dataclass(frozen=True)
class ProductExclusionPolicy:
    exclusion_id: str
    source: str
    product: str
    sha256: str
    size_bytes: str
    header_object: str
    header_observer: str
    header_date_obs: str
    header_date_map: str
    header_frequency_hz: str
    header_naxis1: str
    header_naxis2: str
    header_bmaj_deg: str
    header_bmin_deg: str
    header_bpa_deg: str
    header_clean_niter: str
    header_actnoise_jy_per_beam: str
    free_fit_flux_density_mjy: str
    free_fit_snr: str
    free_fit_ra_j2000: str
    free_fit_dec_j2000: str
    nuclear_ra_j2000: str
    nuclear_dec_j2000: str
    centroid_offset_arcsec: str
    centroid_offset_beam_major: str
    nuclear_photometry_evidence: str
    image_quality_evidence: str
    decision_reason: str
    replacement_reference: str
    replacement_doi: str
    replacement_arxiv: str
    replacement_frequency_ghz: str
    replacement_flux_mjy: str
    replacement_flux_error_mjy: str
    replacement_beam_major_arcsec: str
    replacement_beam_minor_arcsec: str
    replacement_date_label: str


@dataclass(frozen=True)
class SourceRole:
    table_source: str
    source_key: str
    role: str
    note: str


ALIAS_POLICY_COLUMNS = tuple(ProductAliasPolicy.__dataclass_fields__)
EXCLUSION_POLICY_COLUMNS = tuple(ProductExclusionPolicy.__dataclass_fields__)
SOURCE_ROLE_COLUMNS = tuple(SourceRole.__dataclass_fields__)
ALIAS_REPORT_COLUMNS = ALIAS_POLICY_COLUMNS + (
    "fit_telescope",
    "fit_frequency_ghz",
    "fit_flux_density_mjy",
    "fit_date_label",
    "fit_beam_fwhm_arcsec",
    "fit_snr",
)
EXCLUSION_REPORT_COLUMNS = EXCLUSION_POLICY_COLUMNS + (
    "legacy_input_row_present",
    "legacy_fit_telescope",
    "legacy_fit_frequency_ghz",
    "legacy_fit_flux_density_mjy",
    "legacy_fit_date_label",
    "legacy_fit_beam_fwhm_arcsec",
    "legacy_fit_snr",
)


DOI2013_REFERENCE = (
    "Doi et al. 2013, ApJ, 765:63, Tables 1-2; reprocessed VLA data"
)

EXPECTED_SORTED_ROWS = 495
EXPECTED_PREFILTERED_SORTED_ROWS = 486
EXPECTED_CURATED_ROWS = 485
EXPECTED_FINAL_ROWS = 501
PROVENANCE_WIDTH = 146


EXTRAS = [
    ExtraPoint("NGC4258", "VLA", "15.365", "3.5", "0.2", "0.14967", "1997/01/04_1", DOI2013_REFERENCE),
    ExtraPoint("NGC4258", "VLA", "22.460", "3.1", "0.4", "0.09487", "1996/12/31_1", DOI2013_REFERENCE),
    ExtraPoint(
        "NGC4258",
        "SCUBA-JCMT",
        "347",
        "93.7",
        "20.8",
        "15",
        "NA",
        "Doi et al. 2013, ApJ, 765:63, Table 2; nuclear-region flux-density upper bound",
    ),
    ExtraPoint("NGC4258", "NMA", "96", "11", "2.7", "7", "NA", "Doi et al, 2005, MNRAS, 363:2"),
    ExtraPoint("ESO558-G009", "VLA", "33", "0.8", "0.04", "0.19-0.50", "NA", "Kamali, 2017, A&A, 605:A84"),
    ExtraPoint("IC2560", "VLA", "33", "2.0", "0.10", "0.19-0.50", "NA", "Kamali, 2017, A&A, 605:A84"),
    ExtraPoint("MRK1029", "VLA", "33", "1.18", "0.14", "0.19-0.50", "NA", "Kamali, 2017, A&A, 605:A84"),
    ExtraPoint("NGC1194", "VLA", "33", "1.08", "0.04", "0.19-0.50", "NA", "Kamali, 2017, A&A, 605:A84"),
    ExtraPoint("NGC2273", "VLA", "33", "2.69", "0.29", "0.19-0.50", "NA", "Kamali, 2017, A&A, 605:A84"),
    ExtraPoint("NGC3393", "VLA", "33", "5.30", "0.36", "0.19-0.50", "NA", "Kamali, 2017, A&A, 605:A84"),
    ExtraPoint("NGC4388", "VLA", "33", "8.57", "0.43", "0.19-0.50", "NA", "Kamali, 2017, A&A, 605:A84"),
    ExtraPoint("NGC5495", "VLA", "33", "0.13", "0.02", "0.19-0.50", "NA", "Kamali, 2017, A&A, 605:A84"),
    ExtraPoint("UGC3789", "VLA", "33", "0.21", "0.02", "0.19-0.50", "NA", "Kamali, 2017, A&A, 605:A84"),
    ExtraPoint("NGC2960", "VLA", "33", "0.36", "0.06", "0.19-0.50", "NA", "Kamali, 2017, A&A, 605:A84"),
    ExtraPoint(
        "NGC3079",
        "SMA",
        "193.310",
        "53",
        "10.6",
        "4.48",
        "2005/03/01_1",
        "Humphreys et al. 2005, ApJL, 634:L133",
    ),
    ExtraPoint(
        "NGC3079",
        "NOEMA",
        "229.66",
        "14.92",
        "NA",
        "0.503",
        "2018/02/07_1",
        "Lin et al. 2026, ApJ accepted, arXiv:2604.27159",
    ),
]


def parse_python_rows(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    for lineno, raw_line in enumerate(path.read_text().splitlines(), 1):
        line = raw_line.strip()
        if not line:
            continue
        if not line.startswith("["):
            raise ValueError(f"{path}:{lineno}: expected a Python list row")
        row = ast.literal_eval(line)
        if len(row) < 8:
            raise ValueError(f"{path}:{lineno}: expected at least 8 fields, got {len(row)}")
        clean = [str(item).strip().replace("\\n", "") for item in row[:8]]
        clean[7] = portable_release_product(clean[7])
        rows.append(clean)
    return rows


def _read_csv_rows(path: Path, expected_columns: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != expected_columns:
            raise ValueError(
                f"{path}: expected columns {expected_columns!r}, found {reader.fieldnames!r}"
            )
        rows = []
        for lineno, row in enumerate(reader, 2):
            if any(row.get(name) is None for name in expected_columns):
                raise ValueError(f"{path}:{lineno}: missing policy field")
            clean = {name: str(row[name]).strip() for name in expected_columns}
            if any(not value for value in clean.values()):
                raise ValueError(f"{path}:{lineno}: blank policy field")
            rows.append(clean)
    return rows


def load_product_alias_policy(path: Path) -> list[ProductAliasPolicy]:
    entries = [
        ProductAliasPolicy(**row)
        for row in _read_csv_rows(path, ALIAS_POLICY_COLUMNS)
    ]
    if not entries:
        raise ValueError(f"{path}: product-alias policy is empty")

    grouped: dict[str, list[ProductAliasPolicy]] = defaultdict(list)
    seen_products: set[str] = set()
    for entry in entries:
        if entry.role not in {"canonical", "excluded_alias"}:
            raise ValueError(f"{path}: invalid alias role {entry.role!r}")
        if entry.product in seen_products:
            raise ValueError(f"{path}: duplicate product policy for {entry.product}")
        seen_products.add(entry.product)
        if len(entry.sha256) != 64 or any(character not in "0123456789abcdef" for character in entry.sha256):
            raise ValueError(f"{path}: invalid SHA-256 for {entry.product}")
        if int(entry.size_bytes) <= 0:
            raise ValueError(f"{path}: nonpositive file size for {entry.product}")
        grouped[entry.alias_group_id].append(entry)

    matching_evidence_fields = (
        "source",
        "size_bytes",
        "header_observer",
        "header_date_obs",
        "header_frequency_hz",
        "header_naxis1",
        "header_naxis2",
        "header_bmaj_deg",
        "header_bmin_deg",
        "header_clean_niter",
        "pixel_comparison_evidence",
    )
    for group_id, group in grouped.items():
        if sum(entry.role == "canonical" for entry in group) != 1:
            raise ValueError(f"{path}: alias group {group_id!r} must have one canonical product")
        if not any(entry.role == "excluded_alias" for entry in group):
            raise ValueError(f"{path}: alias group {group_id!r} has no excluded alias")
        for field in matching_evidence_fields:
            if len({getattr(entry, field) for entry in group}) != 1:
                raise ValueError(
                    f"{path}: alias group {group_id!r} disagrees in matching evidence {field}"
                )
    return entries


def load_product_exclusion_policy(path: Path) -> list[ProductExclusionPolicy]:
    entries = [
        ProductExclusionPolicy(**row)
        for row in _read_csv_rows(path, EXCLUSION_POLICY_COLUMNS)
    ]
    if not entries:
        raise ValueError(f"{path}: product-exclusion policy is empty")

    seen_ids: set[str] = set()
    seen_products: set[str] = set()
    for entry in entries:
        if entry.exclusion_id in seen_ids:
            raise ValueError(f"{path}: duplicate exclusion_id {entry.exclusion_id!r}")
        if entry.product in seen_products:
            raise ValueError(f"{path}: duplicate excluded product {entry.product!r}")
        seen_ids.add(entry.exclusion_id)
        seen_products.add(entry.product)
        if len(entry.sha256) != 64 or any(
            character not in "0123456789abcdef" for character in entry.sha256
        ):
            raise ValueError(f"{path}: invalid SHA-256 for {entry.product}")
        if int(entry.size_bytes) <= 0:
            raise ValueError(f"{path}: nonpositive file size for {entry.product}")
        positive_fields = (
            "header_frequency_hz",
            "header_naxis1",
            "header_naxis2",
            "header_bmaj_deg",
            "header_bmin_deg",
            "header_clean_niter",
            "header_actnoise_jy_per_beam",
            "free_fit_flux_density_mjy",
            "free_fit_snr",
            "centroid_offset_arcsec",
            "centroid_offset_beam_major",
        )
        for field in positive_fields:
            if float(getattr(entry, field)) <= 0:
                raise ValueError(
                    f"{path}: {field} must be positive for {entry.product}"
                )
        replacement_fields = (
            "replacement_reference",
            "replacement_doi",
            "replacement_arxiv",
            "replacement_frequency_ghz",
            "replacement_flux_mjy",
            "replacement_flux_error_mjy",
            "replacement_beam_major_arcsec",
            "replacement_beam_minor_arcsec",
            "replacement_date_label",
        )
        replacement_values = tuple(getattr(entry, field) for field in replacement_fields)
        if all(value == "NA" for value in replacement_values):
            continue
        if any(value == "NA" for value in replacement_values):
            raise ValueError(
                f"{path}: replacement provenance must be complete or all NA for "
                f"{entry.product}"
            )
        for field in (
            "replacement_frequency_ghz",
            "replacement_flux_mjy",
            "replacement_flux_error_mjy",
            "replacement_beam_major_arcsec",
            "replacement_beam_minor_arcsec",
        ):
            if float(getattr(entry, field)) <= 0:
                raise ValueError(
                    f"{path}: {field} must be positive for {entry.product}"
                )
        if entry.replacement_reference != DOI2013_REFERENCE:
            raise ValueError(
                f"{path}: unexpected replacement reference for {entry.product}"
            )
        if entry.replacement_doi != "10.1088/0004-637X/765/1/63":
            raise ValueError(f"{path}: unexpected replacement DOI for {entry.product}")
        if entry.replacement_arxiv != "1301.4757":
            raise ValueError(f"{path}: unexpected replacement arXiv ID for {entry.product}")
    return entries


def load_source_roles(path: Path) -> list[SourceRole]:
    roles = [SourceRole(**row) for row in _read_csv_rows(path, SOURCE_ROLE_COLUMNS)]
    if len({entry.table_source for entry in roles}) != len(roles):
        raise ValueError(f"{path}: duplicate table_source")
    if len({entry.source_key for entry in roles}) != len(roles):
        raise ValueError(f"{path}: duplicate source_key")
    if any(entry.role not in {"sample", "ancillary"} for entry in roles):
        raise ValueError(f"{path}: roles must be sample or ancillary")
    if sum(entry.role == "sample" for entry in roles) != 21:
        raise ValueError(f"{path}: expected exactly 21 sample sources")
    ancillary = [entry for entry in roles if entry.role == "ancillary"]
    if len(ancillary) != 1 or ancillary[0].table_source != "NGC5495":
        raise ValueError(f"{path}: NGC5495 must be the sole ancillary source")
    return roles


def _observation_signature(row: list[str]) -> tuple[str, ...]:
    return (
        row[0].casefold(),
        row[1].casefold(),
        row[2],
        row[3],
        row[4].partition("_")[0],
        row[5],
    )


def apply_product_alias_policy(
    rows: list[list[str]],
    policy: list[ProductAliasPolicy],
) -> tuple[list[list[str]], dict[str, list[str]]]:
    by_product: dict[str, list[list[str]]] = defaultdict(list)
    for row in rows:
        by_product[row[7]].append(row)

    policy_rows: dict[str, list[str]] = {}
    grouped: dict[str, list[ProductAliasPolicy]] = defaultdict(list)
    for entry in policy:
        matches = by_product.get(entry.product, [])
        if len(matches) != 1:
            raise ValueError(
                f"product-alias policy expected one row for {entry.product!r}, found {len(matches)}"
            )
        row = matches[0]
        if row[0] != entry.source:
            raise ValueError(
                f"product-alias policy source {entry.source!r} does not match row {row[0]!r}"
            )
        policy_rows[entry.product] = row
        grouped[entry.alias_group_id].append(entry)

    for group_id, entries in grouped.items():
        signatures = {
            _observation_signature(policy_rows[entry.product])
            for entry in entries
        }
        if len(signatures) != 1:
            raise ValueError(
                f"product-alias group {group_id!r} does not describe one fitted observation"
            )

    excluded = {entry.product for entry in policy if entry.role == "excluded_alias"}
    canonical = {entry.product for entry in policy if entry.role == "canonical"}
    curated = [row for row in rows if row[7] not in excluded]
    curated_products = {row[7] for row in curated}
    if not canonical.issubset(curated_products):
        raise AssertionError("canonical product missing after archive-alias filtering")
    if excluded.intersection(curated_products):
        raise AssertionError("excluded archive alias survived filtering")
    return curated, policy_rows


def apply_product_exclusion_policy(
    rows: list[list[str]],
    policy: list[ProductExclusionPolicy],
) -> tuple[list[list[str]], dict[str, list[str] | None]]:
    """Remove failed products while preserving any legacy input rows for audit.

    A fresh fitter run skips these products before fitting, so a policy product
    is absent from the new sorted summary.  The frozen historical 409-row
    summary still contains exactly one row for each product; when that legacy
    input is inspected, those rows are captured in the generated audit sidecar
    and excluded from scientific tables.
    """

    by_product: dict[str, list[list[str]]] = defaultdict(list)
    for row in rows:
        by_product[row[7]].append(row)

    policy_rows: dict[str, list[str] | None] = {}
    for entry in policy:
        matches = by_product.get(entry.product, [])
        if len(matches) > 1:
            raise ValueError(
                f"product-exclusion policy expected at most one row for "
                f"{entry.product!r}, found {len(matches)}"
            )
        row = matches[0] if matches else None
        if row is not None and row[0] != entry.source:
            raise ValueError(
                f"product-exclusion policy source {entry.source!r} does not "
                f"match row {row[0]!r}"
            )
        policy_rows[entry.product] = row

    excluded = {entry.product for entry in policy}
    curated = [row for row in rows if row[7] not in excluded]
    if excluded.intersection(row[7] for row in curated):
        raise AssertionError("failed archive product survived exclusion policy")
    return curated, policy_rows


def validate_source_roles(rows: list[list[str]], roles: list[SourceRole]) -> None:
    actual_sources = {row[0] for row in rows}
    expected_sources = {entry.table_source for entry in roles}
    if actual_sources != expected_sources:
        missing = sorted(expected_sources - actual_sources, key=str.casefold)
        extra = sorted(actual_sources - expected_sources, key=str.casefold)
        raise ValueError(f"continuum source-role mismatch: missing={missing}, extra={extra}")


def validate_release_architecture(
    sorted_rows: list[list[str]],
    curated_rows: list[list[str]],
    final_rows: list[list[str]],
    exclusion_policy_rows: dict[str, list[str] | None],
    bandwidth_smearing_policy_rows: dict[str, list[str] | None],
) -> None:
    """Require the reviewed augmented-release row architecture."""

    actual = (len(sorted_rows), len(curated_rows), len(final_rows))
    present_failed_products = sorted(
        product
        for product, row in exclusion_policy_rows.items()
        if row is not None
    )
    if present_failed_products:
        raise ValueError(
            "fresh sorted continuum summary contains upstream-rejected products: "
            + ", ".join(present_failed_products)
        )
    present_smearing_products = sorted(
        product
        for product, row in bandwidth_smearing_policy_rows.items()
        if row is not None
    )
    expected_smearing_products = sorted(bandwidth_smearing_policy_rows)
    reviewed_input = (
        actual
        == (EXPECTED_SORTED_ROWS, EXPECTED_CURATED_ROWS, EXPECTED_FINAL_ROWS)
        and present_smearing_products == expected_smearing_products
    )
    freshly_prefiltered_input = (
        actual
        == (
            EXPECTED_PREFILTERED_SORTED_ROWS,
            EXPECTED_CURATED_ROWS,
            EXPECTED_FINAL_ROWS,
        )
        and not present_smearing_products
    )
    if not (reviewed_input or freshly_prefiltered_input):
        raise ValueError(
            "continuum release architecture must use either the frozen reviewed "
            "495-row input with all nine bandwidth-smearing rows present or a "
            "freshly combined 486-row input with all nine already removed; "
            f"found sorted/curated/final={actual} and "
            f"{len(present_smearing_products)} matching bandwidth-smearing rows"
        )


def write_product_alias_report(
    policy: list[ProductAliasPolicy],
    policy_rows: dict[str, list[str]],
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ALIAS_REPORT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for entry in policy:
            row = policy_rows[entry.product]
            output = {name: getattr(entry, name) for name in ALIAS_POLICY_COLUMNS}
            output.update(
                {
                    "fit_telescope": row[1],
                    "fit_frequency_ghz": row[2],
                    "fit_flux_density_mjy": row[3],
                    "fit_date_label": row[4],
                    "fit_beam_fwhm_arcsec": row[5],
                    "fit_snr": row[6],
                }
            )
            writer.writerow(output)


def write_product_exclusion_report(
    policy: list[ProductExclusionPolicy],
    policy_rows: dict[str, list[str] | None],
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=EXCLUSION_REPORT_COLUMNS,
            lineterminator="\n",
        )
        writer.writeheader()
        for entry in policy:
            row = policy_rows[entry.product]
            output = {name: getattr(entry, name) for name in EXCLUSION_POLICY_COLUMNS}
            output.update(
                {
                    "legacy_input_row_present": "yes" if row is not None else "no",
                    "legacy_fit_telescope": row[1] if row is not None else "NA",
                    "legacy_fit_frequency_ghz": row[2] if row is not None else "NA",
                    "legacy_fit_flux_density_mjy": row[3] if row is not None else "NA",
                    "legacy_fit_date_label": row[4] if row is not None else "NA",
                    "legacy_fit_beam_fwhm_arcsec": row[5] if row is not None else "NA",
                    "legacy_fit_snr": row[6] if row is not None else "NA",
                }
            )
            writer.writerow(output)


def date_key(value: str) -> tuple[int, int, int, int]:
    value = str(value).strip()
    if value.upper() == "NA" or not value:
        return (0, 0, 0, 0)
    ymd, _, suffix = value.partition("_")
    parts = ymd.split("/")
    if len(parts) != 3:
        return (0, 0, 0, 0)
    try:
        year, month, day = (int(part) for part in parts)
        index = int(suffix) if suffix else 0
    except ValueError:
        return (0, 0, 0, 0)
    return (year, month, day, index)


def float_key(value: str) -> float:
    value = str(value).strip().replace("*", "")
    if value.upper() == "NA" or not value:
        return -1.0
    if "-" in value:
        low, high = value.split("-", 1)
        return (float(low) + float(high)) / 2.0
    return float(value)


def source_sort_key(source: str) -> str:
    return source.casefold()


def row_sort_key(row: list[str]) -> tuple[float, tuple[int, int, int, int], float, str]:
    return (float_key(row[2]), date_key(row[4]), float_key(row[5]), row[7])


def canonical_source_lookup(rows: list[list[str]]) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for row in rows:
        lookup.setdefault(row[0].casefold(), row[0])
    return lookup


def extra_to_row(extra: ExtraPoint, source_name: str) -> list[str]:
    if "*" in extra.flux_mjy or extra.flux_err_mjy in {"", "NA"}:
        snr = "NA"
    else:
        snr = f"{float(extra.flux_mjy) / float(extra.flux_err_mjy):.3f}"
    return [
        source_name,
        extra.telescope,
        extra.freq_ghz,
        extra.flux_mjy,
        extra.date,
        extra.beam_arcsec,
        snr,
        extra.reference,
    ]


def build_final_rows(sorted_rows: list[list[str]]) -> list[list[str]]:
    source_lookup = canonical_source_lookup(sorted_rows)
    grouped: dict[str, list[list[str]]] = defaultdict(list)
    for row in sorted_rows:
        grouped[row[0]].append(row)

    inserted: list[ExtraPoint] = []
    for extra in EXTRAS:
        source_name = source_lookup.get(extra.source.casefold(), extra.source)
        grouped[source_name].append(extra_to_row(extra, source_name))
        inserted.append(extra)

    if len(inserted) != len(EXTRAS):
        raise AssertionError(f"Inserted {len(inserted)} extras but expected {len(EXTRAS)}")

    final_rows: list[list[str]] = []
    for source in sorted(grouped, key=source_sort_key):
        final_rows.extend(sorted(grouped[source], key=row_sort_key))

    duplicate_keys = duplicate_measurement_keys(final_rows)
    if duplicate_keys:
        formatted = "\n".join(f"{count} x {key}" for key, count in duplicate_keys)
        raise ValueError(f"Duplicate measurement rows after table build:\n{formatted}")

    return final_rows


def validate_exclusion_replacements(
    rows: list[list[str]],
    policy: list[ProductExclusionPolicy],
) -> None:
    """Require each declared replacement and keep every failed product excluded."""

    for entry in policy:
        if any(row[7] == entry.product for row in rows):
            raise AssertionError(f"excluded archive product survived: {entry.product}")
        if entry.replacement_reference == "NA":
            continue
        beam = (
            float(entry.replacement_beam_major_arcsec)
            * float(entry.replacement_beam_minor_arcsec)
        ) ** 0.5
        expected = [
            entry.source,
            "VLA",
            entry.replacement_frequency_ghz,
            entry.replacement_flux_mjy,
            entry.replacement_date_label,
            f"{beam:.5f}",
            f"{float(entry.replacement_flux_mjy) / float(entry.replacement_flux_error_mjy):.3f}",
            entry.replacement_reference,
        ]
        matches = [row for row in rows if row == expected]
        if len(matches) != 1:
            raise ValueError(
                f"expected one exact literature replacement for {entry.product!r}, "
                f"found {len(matches)}"
            )


def duplicate_measurement_keys(rows: list[list[str]]) -> list[tuple[tuple[str, ...], int]]:
    counts = Counter((row[0], row[1], row[2], row[3], row[4], row[5]) for row in rows)
    return [(key, count) for key, count in counts.items() if count > 1]


def write_final_rows(rows: list[list[str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    grouped: dict[str, list[list[str]]] = defaultdict(list)
    for row in rows:
        grouped[row[0]].append(row)
    with path.open("w") as handle:
        for source in sorted(grouped, key=source_sort_key):
            for row in grouped[source]:
                out = list(row)
                out[7] = out[7] + "\n"
                handle.write(f"{out!r}\n")
            handle.write("\n\n")


def display_file_label(row: list[str]) -> str:
    telescope = row[1]
    raw = row[7].strip()
    if telescope.upper() == "SMA":
        return "SMA_this_work"
    if raw.startswith(("ALMA/", "VLA/")):
        return raw.split("/")[-1]
    if raw.startswith("/"):
        return Path(raw).name
    return raw


def write_fixed_width_machine(rows: list[list[str]], path: Path) -> None:
    widths = [15, 15, 19, 16, 15, 19, 15, PROVENANCE_WIDTH]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for row in rows:
            fields = [
                row[0],
                row[1],
                row[2],
                row[3],
                row[4],
                row[5],
                row[6],
                display_file_label(row),
            ]
            for field, width in zip(fields, widths, strict=True):
                if len(field) > width:
                    raise ValueError(f"Field exceeds width {width}: {field!r}")
            handle.write("".join(field.ljust(width) for field, width in zip(fields, widths, strict=True)).rstrip() + "\n")


def aas_beam_and_flag(beam: str) -> tuple[str, str]:
    if "-" in beam:
        low, high = (float(part) for part in beam.split("-", 1))
        return (f"{(low + high) / 2.0:8.5f}", "*")
    return (f"{float(beam):8.5f}", "")


def aas_flux_and_flag(flux: str) -> tuple[str, str]:
    if "*" in flux:
        return (f"{float(flux.replace('*', '')):7.2f}", "*")
    return (f"{float(flux):7.2f}", "")


def aas_snr(snr: str) -> str:
    if snr == "NA":
        return "NA"
    rendered = f"{float(snr):8.3f}"
    if len(rendered) > 8:
        raise ValueError(f"AAS SNR exceeds the A8 field: {snr!r}")
    return rendered


def write_aas_mrt(rows: list[list[str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = """Title: When the Shadow Meets Its Measure: Assessing the Feasibility of
       Submillimeter Black Hole Shadow Imaging in Megamaser Disk AGN
Authors: Burridge R.N., Bower G.C.
Table: Continuum Flux-density Measurements and Product or Literature Provenance
================================================================================
Byte-by-byte Description of file: datafile2.txt
--------------------------------------------------------------------------------
   Bytes Format Units      Label  Explanations
--------------------------------------------------------------------------------
   1- 12 A12    ---        Name   Source name
  14- 23 A10    ---        Tel    Telescope identifier
  25- 30 F6.2   GHz        Freq   Frequency
  32- 38 F7.2   mJy        Flux   Fitted component or literature flux density
  40- 40 A1     ---      f_Flux   [*] Flag on Flux (1)
  42- 53 A12    ---        Date   Observation date; YYYY/MM/DD_N (2)
  55- 62 F8.5   arcsec     FWHM   Beam Full-Width at Half-Maximum size
  64- 64 A1     ---      f_FWHM   [*] Flag on FWHM (3)
  66- 73 A8     ---        SNR    Signal-to-Noise
  75-220 A146   ---        Prov   Archive product or literature citation (5)
--------------------------------------------------------------------------------
Note (1):
    * = denotes a non-detection.
Note (2): If the number following the underscore in the date (i.e., the
          \"N\" value) is not 1, it distinguishes an additional continuum image
          product for the same source, frequency, and calendar date; it does
          not by itself identify an independent observing epoch.
Note (3):
    * = The beam size is 0.19-0.50 arcseconds.
Note (4): The continuum catalog includes NGC5495 in addition to the 21-source
          analysis sample. Sample-level calculations use the 21-source sample.
Note (5): Image-based rows give the archive product filename or this-work
          identifier; literature-only rows give their published provenance.
--------------------------------------------------------------------------------
"""
    with path.open("w") as handle:
        handle.write(header)
        for row in rows:
            flux, flux_flag = aas_flux_and_flag(row[3])
            beam, beam_flag = aas_beam_and_flag(row[5])
            fields = [
                row[0].ljust(12),
                row[1].ljust(10),
                f"{float(row[2]):6.2f}",
                flux,
                flux_flag.ljust(1),
                row[4].ljust(12),
                beam,
                beam_flag.ljust(1),
                aas_snr(row[6]).ljust(8),
                display_file_label(row).ljust(PROVENANCE_WIDTH),
            ]
            line = " ".join(fields)
            if len(line) != 220:
                raise AssertionError(f"AAS MRT row has length {len(line)}, expected 220: {line!r}")
            handle.write(line + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--sorted-input", type=Path, default=None)
    parser.add_argument("--alias-policy", type=Path, default=None)
    parser.add_argument("--exclusion-policy", type=Path, default=None)
    parser.add_argument("--bandwidth-smearing-policy", type=Path, default=None)
    parser.add_argument("--source-roles", type=Path, default=None)
    parser.add_argument("--final-output", type=Path, default=None)
    parser.add_argument("--machine-output", type=Path, default=None)
    parser.add_argument("--aas-output", type=Path, default=None)
    parser.add_argument("--aliases-output", type=Path, default=None)
    parser.add_argument("--exclusions-output", type=Path, default=None)
    parser.add_argument("--bandwidth-smearing-output", type=Path, default=None)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    base_dir = args.base_dir.resolve()
    sorted_input = args.sorted_input or base_dir / "fitsumfiles" / "sorted" / "fitsummary_fixtochosen_sorted.txt"
    alias_policy_path = args.alias_policy or base_dir / "continuum_product_alias_policy.csv"
    exclusion_policy_path = (
        args.exclusion_policy or base_dir / "continuum_product_exclusion_policy.csv"
    )
    bandwidth_smearing_policy_path = (
        args.bandwidth_smearing_policy
        or base_dir / "continuum_bandwidth_smearing_exclusion_policy.csv"
    )
    source_roles_path = args.source_roles or base_dir / "continuum_source_roles.csv"
    final_output = args.final_output or base_dir / "fitsumfiles" / "final" / "fitsummary_final_withextras.txt"
    machine_output = args.machine_output or base_dir / "machinetables" / "fitsummary_machine.txt"
    aas_output = args.aas_output or base_dir / "machinetables" / "fitsummary_aas_mrt.txt"
    aliases_output = args.aliases_output or base_dir / "machinetables" / "continuum_product_aliases.csv"
    exclusions_output = (
        args.exclusions_output
        or base_dir / "machinetables" / "continuum_product_exclusions.csv"
    )
    bandwidth_smearing_output = (
        args.bandwidth_smearing_output
        or base_dir / "machinetables" / "vla_bandwidth_smearing_exclusions.csv"
    )

    sorted_rows = parse_python_rows(sorted_input)
    alias_policy = load_product_alias_policy(alias_policy_path)
    exclusion_policy = load_product_exclusion_policy(exclusion_policy_path)
    bandwidth_smearing_policy = load_bandwidth_smearing_policy(
        bandwidth_smearing_policy_path
    )
    source_roles = load_source_roles(source_roles_path)
    curated_rows, alias_policy_rows = apply_product_alias_policy(sorted_rows, alias_policy)
    curated_rows, exclusion_policy_rows = apply_product_exclusion_policy(
        curated_rows,
        exclusion_policy,
    )
    curated_rows, bandwidth_smearing_policy_rows = apply_bandwidth_smearing_policy(
        curated_rows,
        bandwidth_smearing_policy,
    )
    final_rows = build_final_rows(curated_rows)
    validate_exclusion_replacements(final_rows, exclusion_policy)
    validate_source_roles(final_rows, source_roles)
    validate_release_architecture(
        sorted_rows,
        curated_rows,
        final_rows,
        exclusion_policy_rows,
        bandwidth_smearing_policy_rows,
    )

    write_final_rows(final_rows, final_output)
    write_fixed_width_machine(final_rows, machine_output)
    write_aas_mrt(final_rows, aas_output)
    write_product_alias_report(alias_policy, alias_policy_rows, aliases_output)
    write_product_exclusion_report(
        exclusion_policy,
        exclusion_policy_rows,
        exclusions_output,
    )
    write_bandwidth_smearing_report(
        bandwidth_smearing_policy,
        bandwidth_smearing_policy_rows,
        bandwidth_smearing_output,
    )

    print(f"Read {len(sorted_rows)} sorted rows")
    print(
        f"Excluded {sum(entry.role == 'excluded_alias' for entry in alias_policy)} "
        f"older archive-product alias(es) under {alias_policy_path}"
    )
    print(
        f"Verified all {len(exclusion_policy)} failed archive products were skipped "
        f"upstream under {exclusion_policy_path}"
    )
    print(
        f"Excluded all {len(bandwidth_smearing_policy)} VLA products that exceed "
        f"x=0.5 under {bandwidth_smearing_policy_path}"
    )
    print(f"Inserted {len(EXTRAS)} literature extras")
    print(
        f"Validated {sum(role.role == 'sample' for role in source_roles)} sample sources and "
        f"{sum(role.role == 'ancillary' for role in source_roles)} additional catalog "
        "source outside the analysis sample"
    )
    print(f"Wrote {len(final_rows)} final rows to {final_output}")
    print(f"Wrote fixed-width table to {machine_output}")
    print(f"Wrote AAS MRT table to {aas_output}")
    print(f"Wrote archive-product alias audit to {aliases_output}")
    print(f"Wrote failed-product exclusion audit to {exclusions_output}")
    print(f"Wrote VLA bandwidth-smearing audit to {bandwidth_smearing_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
