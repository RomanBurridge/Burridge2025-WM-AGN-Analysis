"""Build six canonical manuscript contamination tables deterministically.

An additional NGC 1068 central-only table is retained as a labeled diagnostic.

The calculation is intentionally independent of the publication cells in
``get_contams.ipynb`` and ``Contam_variability.ipynb``.  Those notebooks are
retained as historical records, but this module is the publication path.

Candidate eligibility is defined relative to the fixed 230.0 GHz target:

* thermal dust: input frequency is strictly greater than 1.1 * 230 GHz;
* extended jet: input frequency is strictly less than 0.9 * 230 GHz;
* values exactly on either boundary are ineligible.

All eligible measurements are projected to 230 GHz, and the smallest
projected value is selected.  A selected measured flux gives the maximum
modeled contribution obtained if all of the measured input flux came from
that component.  An input upper limit remains an upper limit and is never
assigned a fabricated symmetric error.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path
from typing import Any, Iterable

from s230_reference import (
    build_s230_reference_rows,
    load_fitsummary_records,
    portable_reference,
    select_chosen_frequency_records,
)


TARGET_FREQUENCY_GHZ = 230.0
HIGH_FREQUENCY_RATIO = 1.1
LOW_FREQUENCY_RATIO = 0.9
# Store the decimal boundary values explicitly so binary evaluation of
# ``1.1 * 230.0`` cannot move the strict 253.0 GHz boundary by one ULP.
HIGH_FREQUENCY_BOUNDARY_GHZ = 253.0
LOW_FREQUENCY_BOUNDARY_GHZ = 207.0
DUST_SPECTRAL_EXPONENT = 4.5
JET_SPECTRAL_EXPONENT = -0.5
MODIFIED_BLACKBODY_BETA = 1.5
MODIFIED_BLACKBODY_TEMPERATURE_MIN_K = 20
MODIFIED_BLACKBODY_TEMPERATURE_MAX_K = 2000

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_SOURCE_TEXT = SCRIPT_DIR / "fitsumfiles/final/fitsummary_final_withextras.txt"
DEFAULT_LITERATURE_INPUTS = SCRIPT_DIR / "contamination_literature_inputs.csv"
DEFAULT_POLICY = SCRIPT_DIR / "contamination_adoption_policy.csv"
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "machinetables/contamination"

INPUT_STATUSES = {"detection", "upper_limit"}
ADOPTION_STATUSES = {"adopted", "not_adopted", "not_available", "pending_choice"}


def _optional_float(value: Any) -> float | None:
    if value in (None, "", "NA"):
        return None
    return float(value)


def _clean_text(value: Any) -> str:
    return str(value).strip()


def classify_candidate(frequency_ghz: float) -> str | None:
    """Return the eligible contaminant class using strict boundaries."""

    frequency_ghz = float(frequency_ghz)
    if frequency_ghz > HIGH_FREQUENCY_BOUNDARY_GHZ:
        return "dust"
    if frequency_ghz < LOW_FREQUENCY_BOUNDARY_GHZ:
        return "jet"
    return None


def projection_scale(
    component: str,
    input_frequency_ghz: float,
    input_beam_arcsec: float,
    target_beam_arcsec: float,
) -> float:
    """Return the exact multiplicative projection factor at 230 GHz."""

    if component == "dust":
        return (
            (TARGET_FREQUENCY_GHZ / float(input_frequency_ghz))
            ** DUST_SPECTRAL_EXPONENT
            * (float(target_beam_arcsec) / float(input_beam_arcsec)) ** 2
        )
    if component == "jet":
        # No beam scaling is applied to the extended-jet prescription.
        return (
            TARGET_FREQUENCY_GHZ / float(input_frequency_ghz)
        ) ** JET_SPECTRAL_EXPONENT
    raise ValueError(f"Unknown contamination component: {component!r}")


def _validate_provenance(record: dict[str, Any]) -> None:
    required = (
        "record_id",
        "source",
        "component",
        "scenario",
        "input_status",
        "telescope",
        "date",
        "reference",
    )
    missing = [name for name in required if not _clean_text(record.get(name, ""))]
    if missing:
        raise ValueError(
            f"Record {record.get('record_id', '<unknown>')} lacks explicit "
            f"provenance/status fields: {', '.join(missing)}"
        )
    if record["input_status"] not in INPUT_STATUSES:
        raise ValueError(
            f"Record {record['record_id']} has invalid input status "
            f"{record['input_status']!r}"
        )
    if record["component"] not in {"dust", "jet"}:
        raise ValueError(f"Record {record['record_id']} has invalid component")
    if float(record["frequency_ghz"]) <= 0 or float(record["beam_arcsec"]) <= 0:
        raise ValueError(f"Record {record['record_id']} has non-positive metadata")
    if float(record["flux_mjy"]) <= 0:
        raise ValueError(f"Record {record['record_id']} has non-positive flux")
    sigma = record.get("sigma_mjy")
    if record["input_status"] == "detection" and sigma is None:
        raise ValueError(f"Detection {record['record_id']} lacks a formal uncertainty")
    if record["input_status"] == "upper_limit" and sigma is not None:
        raise ValueError(
            f"Upper limit {record['record_id']} must not carry a symmetric uncertainty"
        )


def _project_candidate(record: dict[str, Any], target_beam_arcsec: float) -> dict[str, Any]:
    _validate_provenance(record)
    eligible_component = classify_candidate(float(record["frequency_ghz"]))
    if eligible_component != record["component"]:
        raise ValueError(
            f"Record {record['record_id']} is labeled {record['component']!r}, "
            f"but its frequency gives eligibility {eligible_component!r}"
        )

    scale = projection_scale(
        record["component"],
        float(record["frequency_ghz"]),
        float(record["beam_arcsec"]),
        target_beam_arcsec,
    )
    sigma = record.get("sigma_mjy")
    constraint_status = (
        "input_upper_limit"
        if record["input_status"] == "upper_limit"
        else "attribution_upper_bound"
    )
    return {
        **record,
        "target_frequency_ghz": TARGET_FREQUENCY_GHZ,
        "target_beam_arcsec": target_beam_arcsec,
        "scale_factor": scale,
        "projected_flux_mjy": float(record["flux_mjy"]) * scale,
        "projected_sigma_mjy": None if sigma is None else float(sigma) * scale,
        "constraint_status": constraint_status,
    }


def _anchor_signature(record: dict[str, Any]) -> tuple[Any, ...]:
    return (
        record["Source"],
        record["telescope"],
        float(record["freq_ghz"]),
        float(record["flux_mjy"]),
        record["date"],
        record["beam_arcsec"],
        _clean_text(record["reference"]),
    )


def build_continuum_rows(
    source_text_path: str | Path,
    sources: Iterable[str],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, dict[str, Any]]]:
    """Return continuum rows, raw records, and selected reference measurements."""

    records = load_fitsummary_records(source_text_path)
    physical_line_numbers = [
        line_number
        for line_number, line in enumerate(
            Path(source_text_path).read_text().splitlines(), start=1
        )
        if line.strip()
    ]
    if len(physical_line_numbers) != len(records):
        raise ValueError(
            "Parsed continuum records cannot be mapped to physical source-file lines"
        )
    indexed_records: list[dict[str, Any]] = []
    signature_to_id: dict[tuple[Any, ...], str] = {}
    for line_number, record in zip(physical_line_numbers, records, strict=True):
        indexed = {**record, "archive_line_number": line_number}
        indexed_records.append(indexed)
        signature_to_id[_anchor_signature(record)] = f"archive-L{line_number}"

    anchors = select_chosen_frequency_records(records)
    s230_rows = {
        row["Source"]: row for row in build_s230_reference_rows(source_text_path)
    }
    continuum_rows: list[dict[str, Any]] = []
    selected_anchors: dict[str, dict[str, Any]] = {}

    for source in sorted(set(sources), key=str.upper):
        if source not in anchors or source not in s230_rows:
            raise ValueError(
                f"No 200--400 GHz continuum reference measurement exists for {source}"
            )
        anchor = anchors[source]
        selected_anchors[source] = anchor
        s230 = s230_rows[source]
        record_id = signature_to_id[_anchor_signature(anchor)]
        continuum_row = {
            "source": source,
            "record_id": record_id,
            "input_status": anchor["flux_status"],
            "telescope": _clean_text(anchor["telescope"]),
            "date": _clean_text(anchor["date"]),
            "input_frequency_ghz": float(anchor["freq_ghz"]),
            "input_beam_arcsec": float(anchor["beam_arcsec"]),
            "input_flux_mjy": float(anchor["flux_mjy"]),
            "input_sigma_mjy": anchor["sigma_mjy"],
            "reference": _clean_text(s230["nu0_reference"]),
            "target_frequency_ghz": TARGET_FREQUENCY_GHZ,
            "s230_mjy": float(s230["S_230_mJy"]),
            "sigma_s230_mjy": s230["sigma_S_230_mJy"],
            "s230_model_min_mjy": float(s230["S_230_model_min_mJy"]),
            "s230_model_max_mjy": float(s230["S_230_model_max_mJy"]),
            "fiducial_spectral_index": float(s230["alpha_fid"]),
            "spectral_index_range": s230["alpha_model_range"],
        }
        for field in ("record_id", "input_status", "telescope", "date", "reference"):
            if not _clean_text(continuum_row[field]):
                raise ValueError(
                    f"Continuum reference measurement for {source} lacks {field}"
                )
        if continuum_row["input_status"] not in INPUT_STATUSES:
            raise ValueError(
                f"Continuum reference measurement for {source} has invalid status"
            )
        continuum_rows.append(continuum_row)
    return continuum_rows, indexed_records, selected_anchors


def build_archive_candidates(
    indexed_records: Iterable[dict[str, Any]],
    selected_anchors: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    for record in indexed_records:
        source = _clean_text(record["Source"])
        if source not in selected_anchors:
            continue
        component = classify_candidate(float(record["freq_ghz"]))
        if component is None or record["beam_arcsec"] is None:
            continue
        input_status = "upper_limit" if record["is_upper_limit"] else "detection"
        candidate = {
            "record_id": f"archive-L{record['archive_line_number']}",
            "origin": "continuum_archive",
            "source": source,
            "component": component,
            "scenario": "all",
            "input_status": input_status,
            "telescope": _clean_text(record["telescope"]),
            "date": _clean_text(record["date"]),
            "frequency_ghz": float(record["freq_ghz"]),
            "beam_arcsec": float(record["beam_arcsec"]),
            "flux_mjy": float(record["flux_mjy"]),
            "sigma_mjy": None if input_status == "upper_limit" else record["sigma_mjy"],
            "reference": portable_reference(record),
            "citation_key": "",
            "components": "archive measurement",
            "provenance_note": (
                "Exact row from fitsummary_final_withextras.txt; the physical "
                "source-file line number is encoded in record_id."
            ),
        }
        target_beam = float(selected_anchors[source]["beam_arcsec"])
        candidates.append(_project_candidate(candidate, target_beam))
    return candidates


def load_literature_candidates(
    path: str | Path,
    selected_anchors: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    candidates: list[dict[str, Any]] = []
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle):
            source = _clean_text(row["source"])
            if source not in selected_anchors:
                raise ValueError(f"Literature input has unconfigured source {source!r}")
            candidate = {
                "record_id": _clean_text(row["record_id"]),
                "origin": "literature_input",
                "source": source,
                "component": _clean_text(row["component"]),
                "scenario": _clean_text(row["scenario"]),
                "input_status": _clean_text(row["input_status"]),
                "telescope": _clean_text(row["telescope"]),
                "date": _clean_text(row["date"]),
                "frequency_ghz": float(row["frequency_ghz"]),
                "beam_arcsec": float(row["beam_arcsec"]),
                "flux_mjy": float(row["flux_mjy"]),
                "sigma_mjy": _optional_float(row["sigma_mjy"]),
                "reference": _clean_text(row["reference"]),
                "citation_key": _clean_text(row["citation_key"]),
                "components": _clean_text(row["components"]),
                "provenance_note": _clean_text(row["provenance_note"]),
            }
            target_beam = float(selected_anchors[source]["beam_arcsec"])
            candidates.append(_project_candidate(candidate, target_beam))
    return candidates


def load_policy(path: str | Path) -> list[dict[str, str]]:
    policies: list[dict[str, str]] = []
    keys: set[tuple[str, str, str]] = set()
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle):
            policy = {name: _clean_text(value) for name, value in row.items()}
            key = (policy["source"], policy["component"], policy["scenario"])
            if key in keys:
                raise ValueError(f"Duplicate adoption-policy row: {key}")
            keys.add(key)
            if policy["component"] not in {"dust", "jet"}:
                raise ValueError(f"Invalid component in policy row {key}")
            if policy["adoption_status"] not in ADOPTION_STATUSES:
                raise ValueError(f"Invalid adoption status in policy row {key}")
            if not policy["adoption_reason"]:
                raise ValueError(f"Policy row {key} lacks an explicit reason")
            policies.append(policy)
    return policies


def select_candidates(
    candidates: Iterable[dict[str, Any]],
    policies: Iterable[dict[str, str]],
) -> list[dict[str, Any]]:
    candidates = list(candidates)
    selections: list[dict[str, Any]] = []
    for policy in policies:
        compatible = [
            candidate
            for candidate in candidates
            if candidate["source"] == policy["source"]
            and candidate["component"] == policy["component"]
            and candidate["scenario"] in {"all", policy["scenario"]}
        ]
        compatible.sort(key=lambda row: (row["projected_flux_mjy"], row["record_id"]))
        selected = compatible[0] if compatible else None
        if selected is None and policy["adoption_status"] != "not_available":
            raise ValueError(
                "No eligible candidate for non-empty policy row "
                f"{policy['source']}/{policy['component']}/{policy['scenario']}"
            )
        if selected is not None and policy["adoption_status"] == "not_available":
            raise ValueError(
                f"Policy says not available but a candidate exists for {policy['source']} "
                f"{policy['component']} ({policy['scenario']})"
            )
        selections.append(
            {
                **policy,
                "candidate_count": len(compatible),
                "selection": selected,
            }
        )
    return selections


def modified_blackbody_scale(
    input_frequency_ghz: float,
    input_beam_arcsec: float,
    target_beam_arcsec: float,
    temperature_k: float,
    beta: float = MODIFIED_BLACKBODY_BETA,
) -> float:
    """Return optically thin modified-blackbody scaling to 230 GHz."""

    planck_h = 6.62607015e-34
    boltzmann_k = 1.380649e-23
    input_hz = float(input_frequency_ghz) * 1e9
    target_hz = TARGET_FREQUENCY_GHZ * 1e9
    planck_ratio = math.expm1(planck_h * input_hz / (boltzmann_k * temperature_k)) / math.expm1(
        planck_h * target_hz / (boltzmann_k * temperature_k)
    )
    return (
        (TARGET_FREQUENCY_GHZ / float(input_frequency_ghz)) ** (3.0 + beta)
        * planck_ratio
        * (float(target_beam_arcsec) / float(input_beam_arcsec)) ** 2
    )


def build_modified_blackbody_summary(
    selections: Iterable[dict[str, Any]],
    continuum_rows: Iterable[dict[str, Any]],
) -> list[dict[str, Any]]:
    continuum = {row["source"]: row for row in continuum_rows}
    unique_dust: dict[tuple[str, str], dict[str, Any]] = {}
    for item in selections:
        selected = item["selection"]
        if selected is None or item["component"] != "dust":
            continue
        unique_dust[(item["source"], selected["record_id"])] = selected

    summaries: list[dict[str, Any]] = []
    for (source, _), selected in sorted(unique_dust.items()):
        values: list[tuple[int, float, float | None]] = []
        for temperature_k in range(
            MODIFIED_BLACKBODY_TEMPERATURE_MIN_K,
            MODIFIED_BLACKBODY_TEMPERATURE_MAX_K + 1,
        ):
            scale = modified_blackbody_scale(
                selected["frequency_ghz"],
                selected["beam_arcsec"],
                selected["target_beam_arcsec"],
                temperature_k,
            )
            sigma = selected["sigma_mjy"]
            values.append(
                (
                    temperature_k,
                    selected["flux_mjy"] * scale,
                    None if sigma is None else sigma * scale,
                )
            )
        minimum = min(values, key=lambda row: row[1])
        maximum = max(values, key=lambda row: row[1])
        summaries.append(
            {
                "source": source,
                "record_id": selected["record_id"],
                "beta": MODIFIED_BLACKBODY_BETA,
                "temperature_min_k": MODIFIED_BLACKBODY_TEMPERATURE_MIN_K,
                "temperature_max_k": MODIFIED_BLACKBODY_TEMPERATURE_MAX_K,
                "minimum_flux_mjy": minimum[1],
                "minimum_flux_temperature_k": minimum[0],
                "minimum_sigma_mjy": minimum[2],
                "maximum_flux_mjy": maximum[1],
                "maximum_flux_temperature_k": maximum[0],
                "maximum_sigma_mjy": maximum[2],
                "simple_power_law_flux_mjy": selected["projected_flux_mjy"],
                "maximum_fraction_s230_percent": (
                    100.0 * maximum[1] / continuum[source]["s230_mjy"]
                ),
            }
        )
    return summaries


def _flatten_selection(item: dict[str, Any]) -> dict[str, Any]:
    base = {
        "source": item["source"],
        "component": item["component"],
        "scenario": item["scenario"],
        "adoption_status": item["adoption_status"],
        "adoption_reason": item["adoption_reason"],
        "candidate_count": item["candidate_count"],
    }
    selected = item["selection"]
    if selected is None:
        return base
    candidate_fields = {**selected, "candidate_scenario": selected["scenario"]}
    candidate_fields.pop("scenario")
    return {**base, **candidate_fields}


def _csv_value(value: Any) -> Any:
    return "" if value is None else value


def write_csv(rows: list[dict[str, Any]], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(row.get(key)) for key in fieldnames})
    return path


def _latex_escape(value: Any) -> str:
    text = _clean_text(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "_": r"\_",
        "#": r"\#",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def _latex_source(source: str) -> str:
    if source.startswith("NGC"):
        return source.replace("NGC", "NGC~", 1)
    return source


def _latex_telescope(value: Any) -> str:
    """Use a compact publication label while retaining full machine provenance."""

    text = _clean_text(value)
    if text.startswith("HSA ("):
        return "HSA"
    return _latex_escape(text)


def _latex_label_source(source: str) -> str:
    return "circinus" if source == "Circinus" else source


def _publication_date(value: Any) -> str:
    """Remove the internal same-day row suffix from a publication date."""

    return _clean_text(value).partition("_")[0]


def _format_number(value: float, significant_digits: int = 6) -> str:
    decimal_value = Decimal(str(value))
    if not decimal_value.is_finite() or decimal_value == 0:
        return str(decimal_value)
    exponent = decimal_value.copy_abs().adjusted()
    quantum = Decimal("1").scaleb(exponent - significant_digits + 1)
    rounded = decimal_value.quantize(quantum, rounding=ROUND_HALF_UP)
    rounded_exponent = rounded.copy_abs().adjusted()
    if rounded_exponent < -4 or rounded_exponent >= significant_digits:
        return f"{rounded:.{significant_digits - 1}e}"
    decimal_places = max(0, significant_digits - 1 - rounded_exponent)
    text = f"{rounded:.{decimal_places}f}"
    return text.rstrip("0").rstrip(".") if "." in text else text


def _format_latex_number(value: float, significant_digits: int = 6) -> str:
    """Format a number for math mode without Python's literal ``e`` notation."""

    text = _format_number(value, significant_digits)
    if "e" not in text:
        return text
    mantissa, exponent = text.split("e", maxsplit=1)
    return rf"{mantissa}\times10^{{{int(exponent)}}}"


def _decimal_places_for_significant_digits(value: float, significant_digits: int) -> int:
    """Return fixed-point decimal places needed for ``significant_digits``.

    Publication cells use this only for measured values with a reported
    uncertainty.  Computational and machine-readable products retain the
    unrounded floating-point values.
    """

    decimal_value = Decimal(str(value)).copy_abs()
    if not decimal_value.is_finite() or decimal_value <= 0:
        raise ValueError("A finite, positive value is required for display rounding")
    return significant_digits - 1 - decimal_value.adjusted()


def _format_fixed(value: float, decimal_places: int) -> str:
    """Format at an explicit decimal place, preserving meaningful zeros."""

    decimal_value = Decimal(str(value))
    quantum = Decimal("1").scaleb(-decimal_places)
    rounded = decimal_value.quantize(quantum, rounding=ROUND_HALF_UP)
    if decimal_places >= 0:
        return f"{rounded:.{decimal_places}f}"
    return f"{rounded:.0f}"


def _format_value_and_uncertainty(
    value: float,
    uncertainty: float,
    uncertainty_significant_digits: int = 2,
) -> tuple[str, str, int]:
    """Round an uncertainty and its central value for publication display."""

    decimal_places = _decimal_places_for_significant_digits(
        uncertainty, uncertainty_significant_digits
    )
    return (
        _format_fixed(value, decimal_places),
        _format_fixed(uncertainty, decimal_places),
        decimal_places,
    )


def _format_value_with_optional_uncertainty(
    value: float,
    uncertainty: float | None,
) -> str:
    """Render a measured value without inventing an unavailable uncertainty."""

    if uncertainty is None:
        return _format_latex_number(value, 3)
    value_text, uncertainty_text, _ = _format_value_and_uncertainty(
        value, uncertainty
    )
    return rf"{value_text} \pm {uncertainty_text}"


def _format_significant_fixed(value: float, significant_digits: int = 3) -> str:
    """Render a derived value with explicit significant figures and zeros."""

    return _format_fixed(
        value,
        _decimal_places_for_significant_digits(value, significant_digits),
    )


def _format_metadata_number(value: float) -> str:
    """Show observational metadata compactly without changing calculations."""

    return _format_number(value, 3)


def _latex_flux(selected: dict[str, Any]) -> str:
    sigma = selected["projected_sigma_mjy"]
    if sigma is None:
        return rf"$< {_format_latex_number(selected['projected_flux_mjy'], 3)}$"
    flux, sigma_text, _ = _format_value_and_uncertainty(
        selected["projected_flux_mjy"], sigma
    )
    return rf"$\leq ({flux} \pm {sigma_text})$"


def _latex_filename(source: str, scenario: str, canonical: bool = False) -> str:
    if source == "NGC1068":
        if canonical:
            return "contam_NGC1068.tex"
        return f"contam_NGC1068_{scenario}.tex"
    return f"contam_{source}.tex"


def write_latex_tables(
    continuum_rows: list[dict[str, Any]],
    selections: list[dict[str, Any]],
    output_dir: Path,
) -> list[Path]:
    continuum = {row["source"]: row for row in continuum_rows}
    groups: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for item in selections:
        groups.setdefault((item["source"], item["scenario"]), []).append(item)

    latex_dir = output_dir / "latex"
    latex_dir.mkdir(parents=True, exist_ok=True)
    # v1.0.9 resolves the former summed-components scenario into the canonical
    # NGC 1068 table; remove the obsolete scenario-named product on regeneration.
    obsolete_ngc1068 = latex_dir / "contam_NGC1068_summed_components.tex"
    if obsolete_ngc1068.exists():
        obsolete_ngc1068.unlink()
    written: list[Path] = []
    for (source, scenario), items in sorted(groups.items()):
        cont = continuum[source]
        s230_text = _format_value_with_optional_uncertainty(
            cont["s230_mjy"], cont["sigma_s230_mjy"]
        )
        source_tex = _latex_source(source)
        canonical_ngc1068 = (
            source == "NGC1068"
            and scenario == "summed_components"
            and all(item["adoption_status"] == "adopted" for item in items)
        )
        scenario_caption = ""
        if source == "NGC1068":
            if not canonical_ngc1068:
                scenario_caption = " (central-only diagnostic; not adopted)"
        lines = [
            r"\begin{deluxetable*}{lllllll}",
            rf"\tablecaption{{{source_tex}: 230 GHz Continuum Estimates{scenario_caption}}}",
            (
                rf"\label{{tab:{_latex_label_source(source)}_flux}}"
                if scenario == "standard" or canonical_ngc1068
                else rf"\label{{tab:{source}_flux_{scenario}}}"
            ),
            r"\tablehead{\colhead{Source} & \colhead{Emission} & \colhead{Telescope} & \colhead{$S_{230}$ [mJy]} & \colhead{Input Freq. [GHz]} & \colhead{Beam [\arcsec]} & \colhead{Date}}",
            r"\startdata",
            (
                f"{source_tex} & $S_{{230}}$ & {_latex_telescope(cont['telescope'])} & "
                f"${s230_text}$ "
                f"(${_format_significant_fixed(cont['s230_model_min_mjy'])}$--"
                f"${_format_significant_fixed(cont['s230_model_max_mjy'])}$) & "
                f"{_format_number(cont['input_frequency_ghz'])} & "
                f"{_format_metadata_number(cont['input_beam_arcsec'])} & "
                f"{_latex_escape(_publication_date(cont['date']))} \\\\"
            ),
        ]
        included_components: list[str] = []
        for item in sorted(items, key=lambda row: row["component"]):
            selected = item["selection"]
            include_diagnostic = (
                source == "NGC1068"
                and scenario == "central_only"
                and item["component"] == "jet"
                and item["adoption_status"] == "not_adopted"
            )
            if selected is None or (
                item["adoption_status"] == "not_adopted" and not include_diagnostic
            ):
                continue
            included_components.append(item["component"])
            emission = "Thermal dust" if item["component"] == "dust" else "Extended jet"
            if item["adoption_status"] != "adopted":
                emission += f" [{item['adoption_status'].replace('_', ' ')}]"
            lines.append(
                f"{source_tex} & {_latex_escape(emission)} & "
                f"{_latex_telescope(selected['telescope'])} & {_latex_flux(selected)} & "
                f"{_format_number(selected['frequency_ghz'])} & "
                f"{_format_metadata_number(selected['beam_arcsec'])} & "
                f"{_latex_escape(_publication_date(selected['date']))} \\\\"
            )
        if included_components == ["dust", "jet"]:
            contribution_phrase = (
                "The dust and jet rows give maximum modeled contributions"
            )
        elif included_components == ["dust"]:
            contribution_phrase = (
                "The dust row gives the maximum modeled contribution"
            )
        elif included_components == ["jet"]:
            contribution_phrase = (
                "The jet row gives the maximum modeled contribution"
            )
        else:
            raise ValueError(
                f"{source}/{scenario}: publication table has no contamination row"
            )
        comments = (
            "The $S_{230}$ value uses $\\alpha=-0.7$; parentheses give the range "
            "for $-1.0\\leq\\alpha\\leq-0.4$. "
            f"{contribution_phrase} at 230 GHz, "
            "calculated as described in Section~\\ref{subsec:contam}."
        )
        if source == "NGC1068":
            if canonical_ngc1068:
                comments += (
                    " NGC~1068 adopts the total 21.9 GHz HSA S1 model "
                    "($7\\pm1$ mJy; six Gaussian components and two point sources) "
                    "from \\citet{Gallimore2024}; HSA combines the VLBA, GBT, and "
                    "phased VLA."
                )
            else:
                comments += (
                    " This diagnostic uses the 21.9 GHz image peak at NGC~1068* "
                    "from \\citet{Gallimore2024}; HSA combines the VLBA, GBT, and "
                    "phased VLA."
                )
        if source == "NGC4258":
            comments += (
                " The dust estimate uses the SCUBA-JCMT measurement of "
                "\\citet{Doi2013}."
            )
        if source == "NGC3079":
            comments += (
                " The adopted 229.66 GHz NOEMA measurement is reported without a "
                "flux-density uncertainty. The jet "
                "estimate uses the VLBA+VLA measurement of \\citet{Sawada2000}."
            )
        if source == "Circinus":
            comments += (
                " The jet estimate uses the $70\\pm7$ mJy beam$^{-1}$ central "
                "peak measured at 2.37 GHz with ATCA by \\citet{Elmouttie1998}."
            )
        lines.extend(
            [
                r"\enddata",
                rf"\tablecomments{{{comments}}}",
                r"\end{deluxetable*}",
                "",
            ]
        )
        path = latex_dir / _latex_filename(source, scenario, canonical_ngc1068)
        path.write_text("\n".join(lines), encoding="utf-8")
        written.append(path)
    return written


def build_results(
    source_text_path: str | Path = DEFAULT_SOURCE_TEXT,
    literature_inputs_path: str | Path = DEFAULT_LITERATURE_INPUTS,
    policy_path: str | Path = DEFAULT_POLICY,
) -> dict[str, Any]:
    policies = load_policy(policy_path)
    sources = {row["source"] for row in policies}
    continuum_rows, indexed_records, selected_anchors = build_continuum_rows(
        source_text_path, sources
    )
    archive_candidates = build_archive_candidates(indexed_records, selected_anchors)
    literature_candidates = load_literature_candidates(
        literature_inputs_path, selected_anchors
    )
    candidates = sorted(
        archive_candidates + literature_candidates,
        key=lambda row: (
            row["source"],
            row["component"],
            row["projected_flux_mjy"],
            row["record_id"],
        ),
    )
    record_ids = [row["record_id"] for row in candidates]
    if len(record_ids) != len(set(record_ids)):
        raise ValueError("Candidate record IDs are not unique")
    selections = select_candidates(candidates, policies)
    s230_by_source = {row["source"]: row["s230_mjy"] for row in continuum_rows}
    for item in selections:
        selected = item["selection"]
        if selected is not None:
            item["selection"] = {
                **selected,
                "projected_fraction_s230_percent": (
                    100.0
                    * selected["projected_flux_mjy"]
                    / s230_by_source[item["source"]]
                ),
            }
    blackbody = build_modified_blackbody_summary(selections, continuum_rows)
    return {
        "method": {
            "target_frequency_ghz": TARGET_FREQUENCY_GHZ,
            "high_candidate_rule": "frequency_ghz > 1.1 * 230.0 (strict)",
            "low_candidate_rule": "frequency_ghz < 0.9 * 230.0 (strict)",
            "boundary_values_are_eligible": False,
            "dust_formula": "F230 = Fhi * (230/nu_hi)^4.5 * (beam_target/beam_hi)^2",
            "jet_formula": "F230 = Flo * (230/nu_lo)^(-0.5); no beam scaling",
            "formal_error_rule": (
                "multiply the adopted image scaled-MAD or reported literature "
                "input sigma by the same projection factor"
            ),
            "selection_rule": "minimum projected flux among eligible compatible candidates",
            "ngc1068_final_adoption": "summed_components",
            "ngc1068_choice": (
                "summed_components adopted from the source-reported total HSA "
                "S1 model flux; central_only retained as a non-adopted diagnostic"
            ),
        },
        "continuum": continuum_rows,
        "candidates": candidates,
        "selections": selections,
        "modified_blackbody": blackbody,
    }


def write_results(results: dict[str, Any], output_dir: str | Path) -> list[Path]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    paths.append(write_csv(results["continuum"], output_dir / "continuum_anchors.csv"))
    paths.append(write_csv(results["candidates"], output_dir / "contamination_candidates.csv"))
    flat_selections = [_flatten_selection(item) for item in results["selections"]]
    paths.append(write_csv(flat_selections, output_dir / "contamination_selected.csv"))
    paths.append(
        write_csv(
            results["modified_blackbody"],
            output_dir / "modified_blackbody_20_2000K.csv",
        )
    )
    json_path = output_dir / "contamination_results.json"
    json_path.write_text(
        json.dumps(results, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    paths.append(json_path)
    paths.extend(
        write_latex_tables(
            results["continuum"], results["selections"], output_dir
        )
    )
    return paths


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate deterministic 230 GHz contamination products."
    )
    parser.add_argument("--source-text", type=Path, default=DEFAULT_SOURCE_TEXT)
    parser.add_argument(
        "--literature-inputs", type=Path, default=DEFAULT_LITERATURE_INPUTS
    )
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    results = build_results(args.source_text, args.literature_inputs, args.policy)
    paths = write_results(results, args.output_dir)
    print(f"Generated {len(paths)} deterministic contamination products in {args.output_dir}")
    print("NGC 1068 uses the adopted Gallimore et al. (2024) summed-component result.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
