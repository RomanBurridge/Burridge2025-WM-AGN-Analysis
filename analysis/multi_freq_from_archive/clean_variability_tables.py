"""Repair variability machine tables generated with shifted date strings.

The original variability notebook converted dates like ``2016/07/02_1`` to a
decimal year without stripping the ``_1`` suffix. Python accepts underscores in
numeric strings, so ``float("02_1")`` becomes ``21.0``. This script maps the
printed table rows back to the fit-summary rows and restores the true dates.
It also removes warning lines that were captured by stdout redirection.
"""

from __future__ import annotations

import argparse
import ast
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ARCHIVE_DIR = ROOT / "analysis" / "multi_freq_from_archive"
FIT_SUMMARY = ARCHIVE_DIR / "fitsumfiles" / "final" / "fitsummary_final_withextras.txt"
MACHINE_DIR = ARCHIVE_DIR / "machinetables"
FIT_SUMMARY_MACHINE = MACHINE_DIR / "fitsummary_machine.txt"
TABLES = (
    MACHINE_DIR / "variability_machine.txt",
    MACHINE_DIR / "correlations_machine.txt",
)

SPACINGS = [13, 8, 8, 10, 9, 9, 10, 10, 10, 15, 15, 15, 10]


@dataclass(frozen=True)
class FitRecord:
    source: str
    source_key: str
    freq: float
    flux: float
    upper_limit: bool
    beam: float | None
    clean_date: str
    buggy_date: str


def buggy_decimal_year_date(raw_date: str) -> str:
    """Reproduce the date bug from Variability_Search.ipynb."""
    if raw_date == "NA":
        return "NA"

    year_s, month_s, day_s = raw_date.split("/")
    decimal_year = (
        float(year_s)
        + (float(month_s) - 1.0) / 12.0
        + (float(day_s) - 1.0) / 365.0
    )

    year = int(decimal_year)
    year_remainder = decimal_year - year
    dec_month = year_remainder * 12.0 + 1.0
    month = int(dec_month)
    month_remainder = dec_month - month
    day = int(round(month_remainder / 12.0 * 365.0 + 1.0))
    return f"{year}/{month:02}/{day:02}"


def clean_date(raw_date: str) -> str:
    return raw_date.split("_", 1)[0]


def parse_number(value: str) -> float | None:
    value = value.strip().strip("[]").replace("*", "")
    if value == "NA":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def iter_fit_summary_rows(path: Path):
    """Yield source, freq, flux, date, beam rows from either summary format."""
    for raw_line in path.read_text().splitlines():
        if not raw_line.strip():
            continue
        if raw_line.lstrip().startswith("["):
            row = ast.literal_eval(raw_line)
            source, _telescope, freq, flux, date, beam, *_rest = row
            yield source, freq, flux, date, beam
        else:
            parts = raw_line.split(maxsplit=7)
            if len(parts) < 6:
                continue
            source, _telescope, freq, flux, date, beam = parts[:6]
            yield source, freq, flux, date, beam


def load_fit_records(path: Path = FIT_SUMMARY) -> dict[str, list[FitRecord]]:
    if not path.exists():
        path = FIT_SUMMARY_MACHINE
    records: dict[str, list[FitRecord]] = {}
    for source, freq, flux, date, beam in iter_fit_summary_rows(path):
        flux_s = str(flux)
        beam_s = str(beam)
        record = FitRecord(
            source=source,
            source_key=source.upper(),
            freq=float(freq),
            flux=float(flux_s.rstrip("*")),
            upper_limit=flux_s.endswith("*"),
            beam=None if "-" in beam_s else float(beam_s),
            clean_date=clean_date(str(date)),
            buggy_date=buggy_decimal_year_date(str(date)),
        )
        records.setdefault(record.source_key, []).append(record)
    return records


def score_record(
    record: FitRecord,
    freq_s: str,
    beam_s: str,
    flux_s: str,
    printed_date: str,
) -> float:
    score = 0.0
    if printed_date != "NA":
        score += 0.0 if printed_date in {record.buggy_date, record.clean_date} else 1000.0

    freq = parse_number(freq_s)
    if freq is not None:
        # The notebook rounds displayed frequencies heavily; use frequency as a
        # weak discriminator after the buggy printed date, beam, and flux.
        freq_tol = 0.06 if freq < 10 else 0.55
        score += 0.2 * abs(record.freq - freq) / freq_tol

    beam = parse_number(beam_s)
    if beam is not None and record.beam is not None:
        score += abs(record.beam - beam) / max(0.003, 0.05 * max(abs(beam), 1e-9))

    flux = parse_number(flux_s)
    if flux is not None:
        score += abs(record.flux - flux) / max(0.05, 0.06 * max(abs(flux), 1e-9))
        if flux_s.endswith("*") != record.upper_limit:
            score += 2.0

    return score


def match_clean_date(
    records: dict[str, list[FitRecord]],
    source: str,
    freq_s: str,
    beam_s: str,
    flux_s: str,
    printed_date: str,
) -> tuple[str, str]:
    if printed_date == "NA":
        return printed_date, "kept_na"

    candidates = records.get(source.upper(), [])
    if not candidates:
        return printed_date, "missing_source"

    ranked = sorted(
        (
            score_record(record, freq_s, beam_s, flux_s, printed_date),
            idx,
            record,
        )
        for idx, record in enumerate(candidates)
    )
    best_score, _idx, best = ranked[0]
    if best_score >= 1000.0:
        return printed_date, "unmatched"
    return best.clean_date, "fixed" if best.clean_date != printed_date else "unchanged"


def recalc_sep(date1: str, date2: str, old_sep: str) -> str:
    if date1 == "NA" or date2 == "NA":
        return "[NA]"
    try:
        d1 = datetime.strptime(date1, "%Y/%m/%d")
        d2 = datetime.strptime(date2, "%Y/%m/%d")
    except ValueError:
        return old_sep
    return f"[{abs((d2 - d1).days) / 365.25:.3f}]"


def format_row(parts: list[str]) -> str:
    return " ".join(part.ljust(width) for part, width in zip(parts, SPACINGS)).rstrip()


def clean_table(path: Path, records: dict[str, list[FitRecord]]) -> tuple[list[str], dict[str, int]]:
    stats = {
        "data_rows": 0,
        "warning_rows_removed": 0,
        "dates_fixed": 0,
        "dates_unchanged": 0,
        "dates_unmatched": 0,
        "dates_missing_source": 0,
        "dates_kept_na": 0,
        "rows_removed_unmatched": 0,
    }
    cleaned: list[str] = []

    for raw_line in path.read_text().splitlines():
        if not raw_line.strip():
            continue
        if raw_line.startswith("WARNING:"):
            stats["warning_rows_removed"] += 1
            continue

        parts = raw_line.split()
        if not parts:
            continue
        if parts[0] == "target":
            cleaned.append(raw_line.rstrip())
            continue
        if len(parts) < 13:
            cleaned.append(raw_line.rstrip())
            continue

        stats["data_rows"] += 1
        source = parts[0]
        fixed1, status1 = match_clean_date(
            records,
            source,
            freq_s=parts[1],
            beam_s=parts[4],
            flux_s=parts[7],
            printed_date=parts[10],
        )
        fixed2, status2 = match_clean_date(
            records,
            source,
            freq_s=parts[2],
            beam_s=parts[5],
            flux_s=parts[8],
            printed_date=parts[11],
        )

        for status in (status1, status2):
            key = f"dates_{status}"
            if key in stats:
                stats[key] += 1
        if status1 in {"unmatched", "missing_source"} or status2 in {
            "unmatched",
            "missing_source",
        }:
            stats["rows_removed_unmatched"] += 1
            continue
        parts[10] = fixed1
        parts[11] = fixed2
        parts[12] = recalc_sep(fixed1, fixed2, parts[12])
        cleaned.append(format_row(parts))

    return cleaned, stats


def validate_dates(lines: list[str]) -> list[str]:
    bad: list[str] = []
    for line_number, line in enumerate(lines, start=1):
        parts = line.split()
        if len(parts) < 13 or parts[0] == "target":
            continue
        for idx in (10, 11):
            date = parts[idx]
            if date == "NA":
                continue
            try:
                datetime.strptime(date, "%Y/%m/%d")
            except ValueError:
                bad.append(f"line {line_number}: {date}")
    return bad


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="overwrite the existing machine tables instead of writing *_clean.txt files",
    )
    args = parser.parse_args()

    records = load_fit_records()
    report_lines = ["Variability machine-table cleaning report", ""]

    for table in TABLES:
        cleaned, stats = clean_table(table, records)
        bad_dates = validate_dates(cleaned)
        if bad_dates:
            report_lines.append(f"{table.name}: validation failed")
            report_lines.extend(f"  {item}" for item in bad_dates)
            continue

        output = table if args.in_place else table.with_name(f"{table.stem}_clean{table.suffix}")
        output.write_text("\n".join(cleaned) + "\n")
        report_lines.append(f"{table.name} -> {output.name}")
        for key, value in stats.items():
            report_lines.append(f"  {key}: {value}")
        report_lines.append("  date_validation: passed")
        report_lines.append("")

    report_path = MACHINE_DIR / "variability_table_cleaning_report.txt"
    report_path.write_text("\n".join(report_lines).rstrip() + "\n")
    print(report_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
