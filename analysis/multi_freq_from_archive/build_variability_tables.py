#!/usr/bin/env python3

"""Build deterministic pairwise continuum-comparison machine tables.

This script replaces the table-writing sections of ``Variability_Search.ipynb``.
It deliberately has no NED, distance, or black-hole-parameter dependency because
none of those quantities enters a comparison of frequency, beam, flux, or date.

The full correlations table contains every unordered, within-source pair for
the 21-source analysis sample.  NGC 5495 remains in the continuum catalog but
is outside these sample-level calculations under
``continuum_source_roles.csv``.
The historically named variability table is a matched-pair screening subset with a
frequency separation strictly below
1.4 percent and a geometric-mean beam-FWHM separation strictly below 2.6
percent.  Beam ranges are preserved but are not treated as point estimates for
selection.  Same-day pairs are explicitly non-temporal.  Likely duplicate
records that survive the upstream continuum-product policy are reported and
flagged.  The audited older IC 2560 reduction is excluded upstream while both
product identities remain in the continuum alias sidecar. Three failed NGC
4258 NVAS reductions are likewise excluded and preserved in the failed-product
sidecar. Two corresponding Doi et al. (2013) reprocessed measurements enter as
literature detections; the artifact-dominated 2004-05-22 AH0847 image has no
replacement row. Nine additional off-axis VLA measurements are excluded under
the audited bandwidth-smearing policy before these pairs are constructed. This
screen is not a complete quantitative variability analysis.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import math
from collections import defaultdict
from dataclasses import dataclass
from datetime import date
from pathlib import Path, PurePosixPath
from typing import Iterable, Sequence

from continuum_measurements import portable_release_product


FREQUENCY_THRESHOLD_PCT = 1.4
BEAM_THRESHOLD_PCT = 2.6
DAYS_PER_YEAR = 365.25
SOURCE_ROLE_COLUMNS = ("table_source", "source_key", "role", "note")


PAIR_COLUMNS = (
    "pair_id",
    "source",
    "comparison_class",
    "pair_is_likely_duplicate",
    "selected_for_matched_pair_screening",
    "selection_reason",
    "temporal_interpretation",
    "record1_id",
    "telescope1",
    "product1",
    "frequency1_ghz",
    "flux1_density_mjy",
    "flux1_status",
    "flux1_measurement_sigma_mjy",
    "date1_label",
    "beam1_fwhm_input_arcsec",
    "beam1_geometric_mean_fwhm_arcsec",
    "beam1_range_min_arcsec",
    "beam1_range_max_arcsec",
    "record1_likely_duplicate_group",
    "record2_id",
    "telescope2",
    "product2",
    "frequency2_ghz",
    "flux2_density_mjy",
    "flux2_status",
    "flux2_measurement_sigma_mjy",
    "date2_label",
    "beam2_fwhm_input_arcsec",
    "beam2_geometric_mean_fwhm_arcsec",
    "beam2_range_min_arcsec",
    "beam2_range_max_arcsec",
    "record2_likely_duplicate_group",
    "frequency_symmetric_difference_pct",
    "beam_geometric_mean_symmetric_difference_pct",
    "flux_contrast_relation",
    "flux_symmetric_contrast_pct",
    "flux_symmetric_contrast_formal_independent_sigma_pct",
    "date_separation_days",
    "date_separation_years",
)


DUPLICATE_COLUMNS = (
    "likely_duplicate_group",
    "source",
    "reason",
    "record_id",
    "telescope",
    "product",
    "frequency_ghz",
    "flux_density_mjy",
    "flux_status",
    "date_label",
    "beam_fwhm_input_arcsec",
    "fit_snr",
)


SCHEMA_TEXT = """Variability and correlations machine-table schema

Files
-----
variability_machine.txt
    Detection--detection pairs retained by the historical matched-pair screen.
    Screening is strict:
    frequency_symmetric_difference_pct < 1.4 and
    beam_geometric_mean_symmetric_difference_pct < 2.6, and both fluxes must
    be detections with a calculable symmetric contrast.

correlations_machine.txt
    Every unordered, within-source measurement pair for the 21-source analysis
    sample, including upper-limit pairs, pairs containing beam ranges,
    same-day pairs, and likely duplicates. NGC 5495 is outside these
    sample-level calculations under continuum_source_roles.csv.

variability_likely_duplicates.txt
    Records with the same source, telescope, frequency, flux/censoring status,
    calendar date, and beam metadata. These are candidates for manual review;
    they are flagged in, and never silently removed from, either pair table.
    The audited older IC 2560 NVAS reduction is excluded before this stage and
    is preserved with the retained reduction in continuum_product_aliases.csv;
    the three failed NGC 4258 NVAS reductions are preserved in
    continuum_product_exclusions.csv. Two are replaced by explicitly identified
    Doi et al. (2013) literature measurements, while the 2004-05-22 AH0847
    artifact has no replacement. This file is a safety-net report for any
    remaining candidate duplicates.

Conventions
-----------
All files are UTF-8 tab-separated tables with a header row.

The symmetric percentage difference of positive x and y is
    200 * abs(x - y) / (x + y).

For two detections, flux_contrast_relation is ``exact``. Its uncertainty is
the first-order propagation of the two flux-density measurement uncertainties
while treating them as independent. For image-fit inputs, the measurement
uncertainty is the fitted component flux density divided by the final SNR from
the local scaled-MAD classifier; literature-input errors are carried as
reported.

For one detection and one 3-sigma upper limit U below that detection x,
flux_contrast_relation is ``lower_bound`` and the reported contrast is
    200 * (x - U) / (x + U).
No uncertainty is attached to this censored bound. If U is not below x, or if
both values are upper limits, the relation is ``uninformative``. Unavailable
fields are serialized as the explicit string ``NA``.

comparison_class is ``same_day_repeat``, ``different_day``, or
``date_unknown``. Date suffixes such as _1 and _2 identify input records but
do not turn the same calendar day into a variability time baseline.
``temporal_interpretation`` therefore labels same-day rows as non-temporal and
does not call different-day rows confirmed variability.

The scalar beam values in the fit summary are the geometric-mean synthesized
beam FWHM values used by the selection. Literature ranges such as 0.19-0.50
arcsec are retained in the input/min/max columns, have an unavailable scalar
value serialized as ``NA``, and are not eligible for beam-threshold selection.

Stable record and pair IDs are SHA-256-derived identifiers made from the input
content, rather than row numbers, so reordering the input does not change them.
The product1/product2 fields preserve the archive product or literature
identifier that generated each record. Author-local SMA paths are exported as
portable identifiers of the form SMA/track/[clean/]niter/flag/source-image.

This screening table does not scale measurements to one exact frequency and
does not match beam shape, position angle, uv coverage, aperture, or spatial
filtering. Its propagated measurement-scale terms omit absolute
flux-calibration systematics and within-dataset covariance. It is a diagnostic
matched-pair inventory, not a complete quantitative variability measurement or
a demonstrated null test.
"""


@dataclass(frozen=True)
class Beam:
    raw: str
    geometric_mean_fwhm: float | None
    range_min: float | None
    range_max: float | None


@dataclass(frozen=True)
class Measurement:
    record_id: str
    source: str
    telescope: str
    frequency: float
    flux: float
    upper_limit: bool
    date_label: str
    calendar_date: date | None
    beam: Beam
    fit_snr: float | None
    product: str

    @property
    def formal_flux_sigma(self) -> float | None:
        if self.upper_limit or self.fit_snr is None or self.fit_snr <= 0:
            return None
        return self.flux / self.fit_snr


@dataclass(frozen=True)
class FluxContrast:
    relation: str
    percent: float | None
    formal_independent_sigma_percent: float | None


@dataclass(frozen=True)
class Pair:
    first: Measurement
    second: Measurement
    pair_id: str
    frequency_difference_pct: float
    beam_difference_pct: float | None
    flux_contrast: FluxContrast
    comparison_class: str
    date_separation_days: int | None
    selected: bool
    selection_reason: str


def _short_hash(prefix: str, parts: Iterable[str]) -> str:
    canonical = "\x1f".join(parts).encode("utf-8")
    return f"{prefix}-{hashlib.sha256(canonical).hexdigest()[:16]}"


def _normal_float(value: float) -> str:
    return format(value, ".15g")


def _output_float(value: float | None) -> str:
    if value is None:
        return "NA"
    return format(value, ".10g")


def parse_beam(raw_value: object) -> Beam:
    raw = str(raw_value).strip()
    if "-" in raw:
        pieces = raw.split("-", 1)
        if len(pieces) != 2:
            raise ValueError(f"Invalid beam range: {raw!r}")
        low, high = (float(piece) for piece in pieces)
        if low <= 0 or high <= 0 or low > high:
            raise ValueError(f"Invalid beam range: {raw!r}")
        return Beam(raw, None, low, high)
    value = float(raw)
    if value <= 0:
        raise ValueError(f"Beam FWHM must be positive: {raw!r}")
    return Beam(raw, value, None, None)


def parse_calendar_date(raw_value: str) -> date | None:
    clean = raw_value.split("_", 1)[0]
    if clean.upper() == "NA" or not clean:
        return None
    try:
        year, month, day = (int(piece) for piece in clean.split("/"))
        return date(year, month, day)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid observation date {raw_value!r}") from error


def portable_product_identifier(telescope: str, raw_product: str) -> str:
    """Remove author-local path prefixes while retaining reproducible identity."""
    product = portable_release_product(raw_product)
    if telescope.upper() != "SMA":
        return product
    if product.startswith(("SMA/", "SMA_this_work/")):
        return product

    parts = PurePosixPath(product).parts
    try:
        calibration_index = parts.index("sma_calibration")
        track = parts[calibration_index + 1]
        niter = next(part for part in parts[calibration_index + 2 :] if part.startswith("nitercal_"))
        flag = next(part for part in parts[calibration_index + 2 :] if part.startswith("flag"))
    except (StopIteration, ValueError, IndexError):
        return f"SMA_this_work/{PurePosixPath(product).name}"

    portable_parts = ["SMA", track]
    if "clean" in parts[calibration_index + 2 :]:
        portable_parts.append("clean")
    portable_parts.extend((niter, flag, PurePosixPath(product).name))
    return "/".join(portable_parts)


def _record_id(fields: Sequence[str]) -> str:
    return _short_hash("measurement", fields)


def parse_fit_summary(path: Path) -> list[Measurement]:
    measurements: list[Measurement] = []
    seen_ids: set[str] = set()

    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw_line.strip()
        if not line:
            continue
        try:
            row = ast.literal_eval(line)
        except (SyntaxError, ValueError) as error:
            raise ValueError(f"{path}:{line_number}: invalid Python-list row") from error
        if not isinstance(row, list) or len(row) < 8:
            raise ValueError(f"{path}:{line_number}: expected at least eight fields")

        source, telescope, frequency_text, flux_text, date_label, beam_text, snr_text, product = (
            str(value).strip() for value in row[:8]
        )
        product = portable_product_identifier(telescope, product)
        flux_upper_limit = flux_text.endswith("*")
        flux_numeric = flux_text[:-1] if flux_upper_limit else flux_text
        frequency = float(frequency_text)
        flux = float(flux_numeric)
        if frequency <= 0 or flux < 0:
            raise ValueError(f"{path}:{line_number}: frequency must be positive and flux non-negative")
        fit_snr = None if snr_text.upper() == "NA" else float(snr_text)
        if fit_snr is not None and fit_snr <= 0:
            raise ValueError(f"{path}:{line_number}: fit SNR must be positive")
        if flux_upper_limit and fit_snr is not None:
            raise ValueError(f"{path}:{line_number}: an upper limit cannot have a fit SNR")

        identity_fields = (
            source,
            telescope,
            _normal_float(frequency),
            _normal_float(flux),
            "upper_limit_3sigma" if flux_upper_limit else "detection",
            date_label,
            str(beam_text),
            "NA" if fit_snr is None else _normal_float(fit_snr),
            product,
        )
        record_id = _record_id(identity_fields)
        if record_id in seen_ids:
            raise ValueError(f"{path}:{line_number}: exactly duplicated input record {record_id}")
        seen_ids.add(record_id)
        measurements.append(
            Measurement(
                record_id=record_id,
                source=source,
                telescope=telescope,
                frequency=frequency,
                flux=flux,
                upper_limit=flux_upper_limit,
                date_label=date_label,
                calendar_date=parse_calendar_date(date_label),
                beam=parse_beam(beam_text),
                fit_snr=fit_snr,
                product=product,
            )
        )

    return sorted(measurements, key=measurement_sort_key)


def load_source_roles(path: Path) -> dict[str, str]:
    """Load and validate the 21-source sample plus the NGC 5495 catalog row."""
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != SOURCE_ROLE_COLUMNS:
            raise ValueError(
                f"{path}: expected columns {SOURCE_ROLE_COLUMNS!r}, found {reader.fieldnames!r}"
            )
        rows = list(reader)

    roles: dict[str, str] = {}
    source_keys: set[str] = set()
    for line_number, row in enumerate(rows, 2):
        if any(row.get(name) is None for name in SOURCE_ROLE_COLUMNS):
            raise ValueError(f"{path}:{line_number}: missing source-role field")
        clean = {name: str(row[name]).strip() for name in SOURCE_ROLE_COLUMNS}
        if any(not value for value in clean.values()):
            raise ValueError(f"{path}:{line_number}: blank source-role field")
        table_source = clean["table_source"]
        source_key = clean["source_key"]
        role = clean["role"]
        if table_source in roles:
            raise ValueError(f"{path}:{line_number}: duplicate table_source {table_source!r}")
        if source_key in source_keys:
            raise ValueError(f"{path}:{line_number}: duplicate source_key {source_key!r}")
        if role not in {"sample", "ancillary"}:
            raise ValueError(f"{path}:{line_number}: invalid source role {role!r}")
        roles[table_source] = role
        source_keys.add(source_key)

    if sum(role == "sample" for role in roles.values()) != 21:
        raise ValueError(f"{path}: expected exactly 21 sample sources")
    ancillary = [source for source, role in roles.items() if role == "ancillary"]
    if ancillary != ["NGC5495"]:
        raise ValueError(f"{path}: NGC5495 must be the sole ancillary source")
    return roles


def select_sample_measurements(
    measurements: Sequence[Measurement],
    source_roles: dict[str, str],
) -> list[Measurement]:
    """Validate the continuum inventory and retain analysis-sample rows only."""
    actual_sources = {measurement.source for measurement in measurements}
    expected_sources = set(source_roles)
    if actual_sources != expected_sources:
        missing = sorted(expected_sources - actual_sources, key=str.casefold)
        extra = sorted(actual_sources - expected_sources, key=str.casefold)
        raise ValueError(f"continuum source-role mismatch: missing={missing}, extra={extra}")
    return [
        measurement
        for measurement in measurements
        if source_roles[measurement.source] == "sample"
    ]


def measurement_sort_key(item: Measurement) -> tuple[object, ...]:
    return (
        item.source.casefold(),
        item.frequency,
        item.calendar_date or date.min,
        item.date_label,
        item.beam.raw,
        item.telescope.casefold(),
        item.product.casefold(),
        item.record_id,
    )


def symmetric_percentage_difference(first: float, second: float) -> float:
    denominator = first + second
    if denominator == 0:
        return 0.0
    return 200.0 * abs(first - second) / denominator


def detection_contrast_uncertainty(
    first: float,
    second: float,
    first_sigma: float,
    second_sigma: float,
) -> float:
    """First-order sigma for the absolute symmetric percentage contrast."""
    denominator_squared = (first + second) ** 2
    first_term = 400.0 * second * first_sigma / denominator_squared
    second_term = 400.0 * first * second_sigma / denominator_squared
    return math.hypot(first_term, second_term)


def flux_contrast(first: Measurement, second: Measurement) -> FluxContrast:
    if not first.upper_limit and not second.upper_limit:
        percent = symmetric_percentage_difference(first.flux, second.flux)
        first_sigma = first.formal_flux_sigma
        second_sigma = second.formal_flux_sigma
        sigma = None
        if first_sigma is not None and second_sigma is not None:
            sigma = detection_contrast_uncertainty(
                first.flux,
                second.flux,
                first_sigma,
                second_sigma,
            )
        return FluxContrast("exact", percent, sigma)

    if first.upper_limit and second.upper_limit:
        return FluxContrast("uninformative", None, None)

    detection = second if first.upper_limit else first
    upper_limit = first if first.upper_limit else second
    if upper_limit.flux >= detection.flux:
        return FluxContrast("uninformative", None, None)

    lower_bound = symmetric_percentage_difference(detection.flux, upper_limit.flux)
    return FluxContrast("lower_bound", lower_bound, None)


def passes_selected_thresholds(
    frequency_difference_pct: float,
    beam_difference_pct: float | None,
) -> tuple[bool, str]:
    reasons: list[str] = []
    if not frequency_difference_pct < FREQUENCY_THRESHOLD_PCT:
        reasons.append("frequency_not_below_1.4_pct")
    if beam_difference_pct is None:
        reasons.append("beam_range_not_point_comparable")
    elif not beam_difference_pct < BEAM_THRESHOLD_PCT:
        reasons.append("beam_not_below_2.6_pct")
    return (not reasons, "selected" if not reasons else ";".join(reasons))


def build_pair(first: Measurement, second: Measurement) -> Pair:
    if measurement_sort_key(second) < measurement_sort_key(first):
        first, second = second, first
    frequency_difference = symmetric_percentage_difference(first.frequency, second.frequency)
    if first.beam.geometric_mean_fwhm is None or second.beam.geometric_mean_fwhm is None:
        beam_difference = None
    else:
        beam_difference = symmetric_percentage_difference(
            first.beam.geometric_mean_fwhm,
            second.beam.geometric_mean_fwhm,
        )
    selected, selection_reason = passes_selected_thresholds(frequency_difference, beam_difference)
    contrast = flux_contrast(first, second)
    if contrast.relation != "exact":
        selected = False
        censoring_reason = (
            "censored_flux_not_quantitatively_comparable"
            if contrast.relation == "uninformative"
            else "censored_flux_contrast_is_bound_only"
        )
        selection_reason = (
            censoring_reason
            if selection_reason == "selected"
            else f"{selection_reason};{censoring_reason}"
        )

    if first.calendar_date is None or second.calendar_date is None:
        comparison_class = "date_unknown"
        separation_days = None
    else:
        separation_days = abs((second.calendar_date - first.calendar_date).days)
        comparison_class = "same_day_repeat" if separation_days == 0 else "different_day"

    return Pair(
        first=first,
        second=second,
        pair_id=_short_hash("pair", sorted((first.record_id, second.record_id))),
        frequency_difference_pct=frequency_difference,
        beam_difference_pct=beam_difference,
        flux_contrast=contrast,
        comparison_class=comparison_class,
        date_separation_days=separation_days,
        selected=selected,
        selection_reason=selection_reason,
    )


def build_pairs(measurements: Sequence[Measurement]) -> list[Pair]:
    by_source: dict[str, list[Measurement]] = defaultdict(list)
    for measurement in measurements:
        by_source[measurement.source].append(measurement)

    pairs: list[Pair] = []
    for source in sorted(by_source, key=str.casefold):
        source_measurements = sorted(by_source[source], key=measurement_sort_key)
        for first_index, first in enumerate(source_measurements):
            for second in source_measurements[first_index + 1 :]:
                pairs.append(build_pair(first, second))
    return pairs


def likely_duplicate_groups(measurements: Sequence[Measurement]) -> dict[str, str]:
    """Map candidate duplicate record IDs to their review-group ID."""
    signatures: dict[tuple[str, ...], list[Measurement]] = defaultdict(list)
    for item in measurements:
        signatures[
            (
                item.source.casefold(),
                item.telescope.casefold(),
                _normal_float(item.frequency),
                _normal_float(item.flux),
                "upper_limit" if item.upper_limit else "detection",
                item.calendar_date.isoformat() if item.calendar_date else "NA",
                item.beam.raw,
            )
        ].append(item)

    groups: dict[str, str] = {}
    for signature, candidates in signatures.items():
        if len(candidates) < 2:
            continue
        group_id = _short_hash("duplicate", signature)
        for candidate in candidates:
            groups[candidate.record_id] = group_id
    return groups


def _measurement_fields(
    item: Measurement,
    suffix: str,
    duplicate_groups: dict[str, str],
) -> dict[str, str]:
    return {
        f"record{suffix}_id": item.record_id,
        f"telescope{suffix}": item.telescope,
        f"product{suffix}": item.product,
        f"frequency{suffix}_ghz": _output_float(item.frequency),
        f"flux{suffix}_density_mjy": _output_float(item.flux),
        f"flux{suffix}_status": "upper_limit_3sigma" if item.upper_limit else "detection",
        f"flux{suffix}_measurement_sigma_mjy": _output_float(item.formal_flux_sigma),
        f"date{suffix}_label": item.date_label,
        f"beam{suffix}_fwhm_input_arcsec": item.beam.raw,
        f"beam{suffix}_geometric_mean_fwhm_arcsec": _output_float(item.beam.geometric_mean_fwhm),
        f"beam{suffix}_range_min_arcsec": _output_float(item.beam.range_min),
        f"beam{suffix}_range_max_arcsec": _output_float(item.beam.range_max),
        f"record{suffix}_likely_duplicate_group": duplicate_groups.get(item.record_id, "NA"),
    }


def pair_row(pair: Pair, duplicate_groups: dict[str, str]) -> dict[str, str]:
    first_duplicate_group = duplicate_groups.get(pair.first.record_id, "")
    second_duplicate_group = duplicate_groups.get(pair.second.record_id, "")
    pair_is_duplicate = bool(
        first_duplicate_group and first_duplicate_group == second_duplicate_group
    )
    temporal_interpretation = {
        "same_day_repeat": "same_day_not_temporal_variability",
        "different_day": "different_day_screening_candidate_only",
        "date_unknown": "date_unknown_not_temporally_interpretable",
    }[pair.comparison_class]
    row = {
        "pair_id": pair.pair_id,
        "source": pair.first.source,
        "comparison_class": pair.comparison_class,
        "pair_is_likely_duplicate": "yes" if pair_is_duplicate else "no",
        "selected_for_matched_pair_screening": "yes" if pair.selected else "no",
        "selection_reason": pair.selection_reason,
        "temporal_interpretation": temporal_interpretation,
        "frequency_symmetric_difference_pct": _output_float(pair.frequency_difference_pct),
        "beam_geometric_mean_symmetric_difference_pct": _output_float(pair.beam_difference_pct),
        "flux_contrast_relation": pair.flux_contrast.relation,
        "flux_symmetric_contrast_pct": _output_float(pair.flux_contrast.percent),
        "flux_symmetric_contrast_formal_independent_sigma_pct": _output_float(
            pair.flux_contrast.formal_independent_sigma_percent
        ),
        "date_separation_days": "NA" if pair.date_separation_days is None else str(pair.date_separation_days),
        "date_separation_years": (
            "NA"
            if pair.date_separation_days is None
            else _output_float(pair.date_separation_days / DAYS_PER_YEAR)
        ),
    }
    row.update(_measurement_fields(pair.first, "1", duplicate_groups))
    row.update(_measurement_fields(pair.second, "2", duplicate_groups))
    return row


def write_pair_table(
    pairs: Iterable[Pair],
    path: Path,
    duplicate_groups: dict[str, str],
) -> int:
    rows = [pair_row(pair, duplicate_groups) for pair in pairs]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PAIR_COLUMNS, dialect="excel-tab", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def write_duplicate_report(
    measurements: Sequence[Measurement],
    duplicate_groups: dict[str, str],
    path: Path,
) -> int:
    candidates = [item for item in measurements if item.record_id in duplicate_groups]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=DUPLICATE_COLUMNS, dialect="excel-tab", lineterminator="\n")
        writer.writeheader()
        for item in candidates:
            writer.writerow(
                {
                    "likely_duplicate_group": duplicate_groups[item.record_id],
                    "source": item.source,
                    "reason": "same_source_instrument_frequency_flux_calendar_date_and_beam",
                    "record_id": item.record_id,
                    "telescope": item.telescope,
                    "product": item.product,
                    "frequency_ghz": _output_float(item.frequency),
                    "flux_density_mjy": _output_float(item.flux),
                    "flux_status": "upper_limit_3sigma" if item.upper_limit else "detection",
                    "date_label": item.date_label,
                    "beam_fwhm_input_arcsec": item.beam.raw,
                    "fit_snr": _output_float(item.fit_snr),
                }
            )
    return len(candidates)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--input", type=Path, default=None)
    parser.add_argument("--source-roles", type=Path, default=None)
    parser.add_argument("--variability-output", type=Path, default=None)
    parser.add_argument("--correlations-output", type=Path, default=None)
    parser.add_argument("--duplicates-output", type=Path, default=None)
    parser.add_argument("--schema-output", type=Path, default=None)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    base_dir = args.base_dir.resolve()
    machine_dir = base_dir / "machinetables"
    input_path = args.input or base_dir / "fitsumfiles" / "final" / "fitsummary_final_withextras.txt"
    source_roles_path = args.source_roles or base_dir / "continuum_source_roles.csv"
    variability_output = args.variability_output or machine_dir / "variability_machine.txt"
    correlations_output = args.correlations_output or machine_dir / "correlations_machine.txt"
    duplicates_output = args.duplicates_output or machine_dir / "variability_likely_duplicates.txt"
    schema_output = args.schema_output or machine_dir / "variability_tables_schema.txt"

    all_measurements = parse_fit_summary(input_path)
    source_roles = load_source_roles(source_roles_path)
    measurements = select_sample_measurements(all_measurements, source_roles)
    duplicate_groups = likely_duplicate_groups(measurements)
    pairs = build_pairs(measurements)
    selected_pairs = [pair for pair in pairs if pair.selected]

    correlation_count = write_pair_table(pairs, correlations_output, duplicate_groups)
    selected_count = write_pair_table(selected_pairs, variability_output, duplicate_groups)
    duplicate_record_count = write_duplicate_report(measurements, duplicate_groups, duplicates_output)
    schema_output.parent.mkdir(parents=True, exist_ok=True)
    schema_output.write_text(SCHEMA_TEXT, encoding="utf-8")

    same_day_count = sum(pair.comparison_class == "same_day_repeat" for pair in selected_pairs)
    geometry_matched_censored_count = sum(
        pair.flux_contrast.relation != "exact"
        and passes_selected_thresholds(pair.frequency_difference_pct, pair.beam_difference_pct)[0]
        for pair in pairs
    )
    print(f"Read {len(all_measurements)} continuum measurements from {input_path}")
    print(
        f"Retained {len(measurements)} measurements for the 21-source analysis sample "
        f"under {source_roles_path}"
    )
    print(f"Wrote {correlation_count} within-source pairs to {correlations_output}")
    print(f"Wrote {selected_count} matched-pair screening rows to {variability_output}")
    print(
        f"Screening rows include {same_day_count} explicitly non-temporal same-day "
        "and zero censored comparisons"
    )
    print(
        f"Retained {geometry_matched_censored_count} geometry-matched censored pairs only in "
        f"{correlations_output}"
    )
    print(
        f"Flagged {duplicate_record_count} records in {len(set(duplicate_groups.values()))} "
        f"likely-duplicate group(s) in {duplicates_output}"
    )
    print(f"Wrote schema to {schema_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
