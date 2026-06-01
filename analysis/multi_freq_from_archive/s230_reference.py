from __future__ import annotations

import argparse
import ast
import csv
import math
from pathlib import Path
from typing import Iterable


FREQ_REF_GHZ = 230.0
ALPHA_FID = -0.7
ALPHA_MINUS_1P0 = -1.0
ALPHA_MINUS_0P4 = -0.4
ALPHA_MODEL_RANGE = "[-1.0, -0.4]"
FAR_NU0_RELATIVE_THRESHOLD = 0.3


def linear_to_log_uncertainty(lin_x: float, lin_delta_x: float) -> float:
    return lin_delta_x / (lin_x * math.log(10))


def _is_upper_limit(flux_field: str) -> bool:
    return "*" in str(flux_field)


def _safe_float(value):
    if value in ("", "NA", None):
        return None
    return float(value)


def _parse_beam(value):
    text = str(value).strip()
    if "-" in text or text in ("", "NA"):
        return None
    return float(text)


def _parse_snr(value):
    if value in ("", "NA", None):
        return None
    return float(str(value).splitlines()[0])


def _parse_record(line: str) -> dict[str, object]:
    row = ast.literal_eval(line.strip())
    flux_field = str(row[3]).strip()
    flux_mjy = float(flux_field.split("*")[0])
    is_upper_limit = _is_upper_limit(flux_field)
    snr = _parse_snr(row[6])
    sigma_mjy = None
    if snr not in (None, 0):
        sigma_mjy = flux_mjy / snr

    return {
        "Source": row[0],
        "telescope": row[1],
        "freq_ghz": float(row[2]),
        "flux_field": flux_field,
        "flux_mjy": flux_mjy,
        "is_upper_limit": is_upper_limit,
        "date": row[4],
        "beam_arcsec": _parse_beam(row[5]),
        "snr": snr,
        "sigma_mjy": sigma_mjy,
        "reference": row[7],
    }


def load_fitsummary_records(source_text_path: str | Path) -> list[dict[str, object]]:
    path = Path(source_text_path)
    records: list[dict[str, object]] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            records.append(_parse_record(line))
    return records


def select_chosen_frequency_records(
    records: Iterable[dict[str, object]],
    freq_min_ghz: float = 200.0,
    freq_max_ghz: float = 400.0,
) -> dict[str, dict[str, object]]:
    grouped: dict[str, dict[str, list[dict[str, object]]]] = {}

    for record in records:
        freq_ghz = float(record["freq_ghz"])
        beam_arcsec = record["beam_arcsec"]
        if beam_arcsec is None:
            continue
        if not (freq_min_ghz < freq_ghz < freq_max_ghz):
            continue
        source = str(record["Source"])
        grouped.setdefault(source, {"detections": [], "upper_limits": []})
        if bool(record["is_upper_limit"]):
            grouped[source]["upper_limits"].append(record)
        else:
            grouped[source]["detections"].append(record)

    chosen: dict[str, dict[str, object]] = {}
    for source, candidates in grouped.items():
        if candidates["detections"]:
            best = min(candidates["detections"], key=lambda rec: float(rec["beam_arcsec"]))
            chosen[source] = {**best, "flux_status": "detection"}
        elif candidates["upper_limits"]:
            best = min(candidates["upper_limits"], key=lambda rec: float(rec["beam_arcsec"]))
            chosen[source] = {**best, "flux_status": "upper_limit"}
    return chosen


def _compute_s230_fields(record: dict[str, object]) -> dict[str, object]:
    flux_status = str(record.get("flux_status", "detection"))
    nu0_ghz = float(record["freq_ghz"])
    s_nu0_mjy = float(record["flux_mjy"])
    sigma_s_nu0_mjy = None if flux_status == "upper_limit" else record["sigma_mjy"]

    fid_scale = (FREQ_REF_GHZ / nu0_ghz) ** ALPHA_FID
    s_230_mjy = s_nu0_mjy * fid_scale
    sigma_s_230_mjy = None
    if sigma_s_nu0_mjy is not None:
        sigma_s_230_mjy = float(sigma_s_nu0_mjy) * fid_scale

    s_230_alpha_minus_1p0_mjy = s_nu0_mjy * (FREQ_REF_GHZ / nu0_ghz) ** ALPHA_MINUS_1P0
    s_230_alpha_minus_0p4_mjy = s_nu0_mjy * (FREQ_REF_GHZ / nu0_ghz) ** ALPHA_MINUS_0P4
    s_230_model_min_mjy = min(s_230_alpha_minus_1p0_mjy, s_230_alpha_minus_0p4_mjy)
    s_230_model_max_mjy = max(s_230_alpha_minus_1p0_mjy, s_230_alpha_minus_0p4_mjy)

    percent_change_to_230 = 100.0 * (s_230_mjy - s_nu0_mjy) / s_nu0_mjy
    model_range_percent_width = 100.0 * (s_230_model_max_mjy - s_230_model_min_mjy) / s_230_mjy

    return {
        "Source": record["Source"],
        "flux_status": flux_status,
        "nu0_GHz": nu0_ghz,
        "S_nu0_mJy": s_nu0_mjy,
        "sigma_S_nu0_mJy": sigma_s_nu0_mjy,
        "alpha_fid": ALPHA_FID,
        "S_230_mJy": s_230_mjy,
        "sigma_S_230_mJy": sigma_s_230_mjy,
        "alpha_model_range": ALPHA_MODEL_RANGE,
        "S_230_alpha_minus_1p0_mJy": s_230_alpha_minus_1p0_mjy,
        "S_230_alpha_minus_0p4_mJy": s_230_alpha_minus_0p4_mjy,
        "S_230_model_min_mJy": s_230_model_min_mjy,
        "S_230_model_max_mJy": s_230_model_max_mjy,
        "percent_change_to_230": percent_change_to_230,
        "model_range_percent_width": model_range_percent_width,
        "nu0_telescope": record["telescope"],
        "nu0_date": record["date"],
        "nu0_beam_arcsec": record["beam_arcsec"],
    }


def build_s230_reference_rows(
    source_text_path: str | Path,
    output_csv_path: str | Path | None = None,
    verbose: bool = False,
) -> list[dict[str, object]]:
    records = load_fitsummary_records(source_text_path)
    chosen = select_chosen_frequency_records(records)
    rows = [_compute_s230_fields(chosen[source]) for source in sorted(chosen, key=str.upper)]

    if output_csv_path is not None:
        write_s230_reference_csv(rows, output_csv_path)

    if verbose:
        print_s230_reference_table(rows)
        print_diagnostic_warnings(rows)

    return rows


def _format_csv_value(value):
    if value is None:
        return "NA"
    return value


def write_s230_reference_csv(rows: list[dict[str, object]], output_csv_path: str | Path) -> Path:
    path = Path(output_csv_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "Source",
        "flux_status",
        "nu0_GHz",
        "S_nu0_mJy",
        "sigma_S_nu0_mJy",
        "alpha_fid",
        "S_230_mJy",
        "sigma_S_230_mJy",
        "alpha_model_range",
        "S_230_alpha_minus_1p0_mJy",
        "S_230_alpha_minus_0p4_mJy",
        "S_230_model_min_mJy",
        "S_230_model_max_mJy",
        "percent_change_to_230",
        "model_range_percent_width",
    ]

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: _format_csv_value(row.get(name)) for name in fieldnames})

    return path


def _fmt(value, precision: int = 2) -> str:
    if value is None:
        return "NA"
    return f"{float(value):.{precision}f}"


def print_s230_reference_table(rows: list[dict[str, object]]) -> None:
    print("S230 chosen-frequency extrapolation comparison")
    header = (
        f"{'Source':<14} "
        f"{'status':<12} "
        f"{'nu0_GHz':>8} "
        f"{'S_nu0':>10} "
        f"{'sigma_nu0':>10} "
        f"{'S_230':>10} "
        f"{'sigma_230':>10} "
        f"{'%chg230':>10} "
        f"{'model_min':>10} "
        f"{'model_max':>10} "
        f"{'model_w%':>10}"
    )
    print(header)
    print("-" * len(header))

    for row in rows:
        print(
            f"{str(row['Source']):<14} "
            f"{str(row['flux_status']):<12} "
            f"{_fmt(row['nu0_GHz'], 2):>8} "
            f"{_fmt(row['S_nu0_mJy'], 2):>10} "
            f"{_fmt(row['sigma_S_nu0_mJy'], 2):>10} "
            f"{_fmt(row['S_230_mJy'], 2):>10} "
            f"{_fmt(row['sigma_S_230_mJy'], 2):>10} "
            f"{_fmt(row['percent_change_to_230'], 2):>10} "
            f"{_fmt(row['S_230_model_min_mJy'], 2):>10} "
            f"{_fmt(row['S_230_model_max_mJy'], 2):>10} "
            f"{_fmt(row['model_range_percent_width'], 2):>10}"
        )


def print_diagnostic_warnings(rows: list[dict[str, object]]) -> None:
    warnings: list[str] = []

    for row in rows:
        nu0_ghz = float(row["nu0_GHz"])
        if abs(nu0_ghz - FREQ_REF_GHZ) / FREQ_REF_GHZ > FAR_NU0_RELATIVE_THRESHOLD:
            warnings.append(
                f"WARNING: {row['Source']} has nu0={nu0_ghz:.2f} GHz, "
                "which is more than 30% from 230 GHz; the spectral-index model bracket is more important."
            )

    if any(str(row["Source"]).upper() == "NGC4258" for row in rows):
        warnings.append(
            "NOTE: NGC4258 uses the selected 200--400 GHz anchor record in this table. "
            "Check any manuscript prose against the adopted anchor table before submission."
        )

    if not warnings:
        print("No S230 diagnostic warnings.")
        return

    print("\nS230 diagnostic warnings")
    for warning in warnings:
        print(warning)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a shared 230 GHz chosen-frequency extrapolation comparison table."
    )
    parser.add_argument(
        "--source-text",
        default="fitsumfiles/final/fitsummary_final_withextras.txt",
        help="Path to the final continuum summary text file.",
    )
    parser.add_argument(
        "--output",
        default="machinetables/S230_extrapolation_comparison.csv",
        help="Path to the output comparison CSV.",
    )
    args = parser.parse_args()

    rows = build_s230_reference_rows(
        source_text_path=args.source_text,
        output_csv_path=args.output,
        verbose=False,
    )
    output_path = Path(args.output).resolve()
    print_s230_reference_table(rows)
    print_diagnostic_warnings(rows)
    print(f"\nWrote {output_path}")


if __name__ == "__main__":
    main()
