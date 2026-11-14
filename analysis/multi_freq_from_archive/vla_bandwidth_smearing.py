#!/usr/bin/env python3

"""Apply the audited VLA bandwidth-smearing exclusion policy."""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class BandwidthSmearingExclusion:
    exclusion_id: str
    source: str
    product: str
    catalog_frequency_ghz: str
    catalog_flux_mjy: str
    catalog_date_label: str
    catalog_beam_fwhm_arcsec: str
    catalog_snr: str
    center_frequency_mhz: str
    channel_width_mhz: str
    channel_width_evidence: str
    source_offset_arcsec: str
    beam_major_fwhm_arcsec: str
    beam_minor_fwhm_arcsec: str
    beam_pa_deg: str
    radial_position_angle_deg: str
    radial_beam_fwhm_arcsec: str
    smearing_parameter_x: str
    maximum_accepted_x: str
    decision_reason: str


POLICY_COLUMNS = tuple(BandwidthSmearingExclusion.__dataclass_fields__)
REPORT_COLUMNS = POLICY_COLUMNS + (
    "input_row_present",
    "input_telescope",
    "input_frequency_ghz",
    "input_flux_mjy",
    "input_date_label",
    "input_beam_fwhm_arcsec",
    "input_snr",
)


def load_policy(path: Path) -> list[BandwidthSmearingExclusion]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != POLICY_COLUMNS:
            raise ValueError(
                f"{path}: expected columns {POLICY_COLUMNS!r}, found {reader.fieldnames!r}"
            )
        raw_rows = []
        for lineno, row in enumerate(reader, 2):
            clean = {name: str(row.get(name, "")).strip() for name in POLICY_COLUMNS}
            if any(not value for value in clean.values()):
                raise ValueError(f"{path}:{lineno}: blank bandwidth-smearing field")
            raw_rows.append(clean)

    entries = [BandwidthSmearingExclusion(**row) for row in raw_rows]
    if not entries:
        raise ValueError(f"{path}: bandwidth-smearing policy is empty")
    if len({entry.exclusion_id for entry in entries}) != len(entries):
        raise ValueError(f"{path}: duplicate exclusion_id")
    if len({entry.product for entry in entries}) != len(entries):
        raise ValueError(f"{path}: duplicate product")

    for entry in entries:
        if not entry.product.startswith(f"VLA/{entry.source}/"):
            raise ValueError(
                f"{path}: product/source mismatch for {entry.exclusion_id}"
            )
        center_frequency = float(entry.center_frequency_mhz)
        channel_width = float(entry.channel_width_mhz)
        source_offset = float(entry.source_offset_arcsec)
        radial_beam = float(entry.radial_beam_fwhm_arcsec)
        reported_x = float(entry.smearing_parameter_x)
        maximum_x = float(entry.maximum_accepted_x)
        if min(center_frequency, channel_width, source_offset, radial_beam, maximum_x) <= 0:
            raise ValueError(f"{path}: nonpositive smearing input for {entry.exclusion_id}")
        calculated_x = (channel_width / center_frequency) * (
            source_offset / radial_beam
        )
        if not math.isclose(calculated_x, reported_x, rel_tol=0.0, abs_tol=1e-6):
            raise ValueError(
                f"{path}: x mismatch for {entry.exclusion_id}: "
                f"calculated {calculated_x:.9f}, reported {reported_x:.9f}"
            )
        if reported_x <= maximum_x:
            raise ValueError(
                f"{path}: exclusion {entry.exclusion_id} does not exceed x={maximum_x}"
            )
    return entries


def _validate_catalog_row(
    entry: BandwidthSmearingExclusion,
    row: list[str],
) -> None:
    expected = (
        entry.source,
        "VLA",
        entry.catalog_frequency_ghz,
        entry.catalog_flux_mjy,
        entry.catalog_date_label,
        entry.catalog_beam_fwhm_arcsec,
        entry.catalog_snr,
        entry.product,
    )
    actual = tuple(str(value).strip() for value in row[:8])
    if actual != expected:
        raise ValueError(
            f"bandwidth-smearing policy row mismatch for {entry.product}: "
            f"expected {expected!r}, found {actual!r}"
        )


def apply_policy(
    rows: list[list[str]],
    policy: list[BandwidthSmearingExclusion],
) -> tuple[list[list[str]], dict[str, list[str] | None]]:
    """Exclude policy products and return their matched input rows for audit."""

    by_product: dict[str, list[list[str]]] = defaultdict(list)
    for row in rows:
        by_product[str(row[7]).strip()].append(row)

    matched: dict[str, list[str] | None] = {}
    for entry in policy:
        matches = by_product.get(entry.product, [])
        if len(matches) > 1:
            raise ValueError(
                f"bandwidth-smearing policy expected at most one row for "
                f"{entry.product!r}, found {len(matches)}"
            )
        row = matches[0] if matches else None
        if row is not None:
            _validate_catalog_row(entry, row)
        matched[entry.product] = row

    excluded = {entry.product for entry in policy}
    curated = [row for row in rows if str(row[7]).strip() not in excluded]
    if excluded.intersection(str(row[7]).strip() for row in curated):
        raise AssertionError("bandwidth-smeared VLA product survived exclusion policy")
    return curated, matched


def write_report(
    policy: list[BandwidthSmearingExclusion],
    matched: dict[str, list[str] | None],
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_COLUMNS)
        writer.writeheader()
        for entry in policy:
            row = matched[entry.product]
            output = {name: getattr(entry, name) for name in POLICY_COLUMNS}
            output.update(
                {
                    "input_row_present": "yes" if row is not None else "no",
                    "input_telescope": row[1] if row is not None else "NA",
                    "input_frequency_ghz": row[2] if row is not None else "NA",
                    "input_flux_mjy": row[3] if row is not None else "NA",
                    "input_date_label": row[4] if row is not None else "NA",
                    "input_beam_fwhm_arcsec": row[5] if row is not None else "NA",
                    "input_snr": row[6] if row is not None else "NA",
                }
            )
            writer.writerow(output)
