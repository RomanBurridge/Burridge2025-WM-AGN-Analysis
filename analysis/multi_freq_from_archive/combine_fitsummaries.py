#!/usr/bin/env python3

"""Combine per-run fixed-to-chosen continuum summaries for post-processing.

This helper builds the legacy input file expected by
`post_processing_finalcontinuum.ipynb`:

    fitsumfiles/sorted/fitsummary_fixtochosen_sorted.txt

It prefers already-sorted per-run files when they exist and falls back to the
corresponding `witherrors` summaries when a run crashed before the later
noerrors/sorted stages.
"""

from __future__ import annotations

import argparse
import ast
import re
from collections import defaultdict
from pathlib import Path


SORTED_SUFFIX = "_fitsummary_fixedtochosen_sorted.txt"
WITHERRORS_SUFFIX = "_fitsummary_fixedtochosen_witherrors.txt"
OUTPUT_NAME = "fitsummary_fixtochosen_sorted.txt"

# Manual row exclusions for the downstream continuum summary.
# For NGC4258 we keep the 2023-03-16 SMA detection and drop the older
# 2022-01-31 SMA upper-limit entry from the downstream combined list.
EXCLUDED_ROWS = {
    ("NGC4258", "SMA", "2022/01/31_1"),
}


def summary_prefix(path: Path, suffix: str) -> str:
    return path.name[: -len(suffix)]


def date_key(value: str) -> float:
    value = str(value).strip().strip(",")
    if value.upper() == "NA" or not value:
        return -1.0
    parts = value.split("/")
    if len(parts) != 3:
        return -1.0
    try:
        year = float(parts[0])
        month = float(parts[1])
        day_part = parts[2].split("_")[0]
        day = float(day_part)
    except ValueError:
        return -1.0
    return year + month / 12.0 + day / 365.0


def float_key(value: str) -> float:
    value = str(value).strip()
    if value.upper() == "NA" or not value:
        return -1.0
    return float(value)


def parse_sorted_file(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line or not line.startswith("["):
            continue
        row = ast.literal_eval(line)
        if len(row) < 8:
            continue
        rows.append([str(item) for item in row[:8]])
    return rows


def parse_witherrors_file(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("target"):
            continue
        if line == "images/chosenfreq":
            continue
        if "Error" in line or "pixels" in line:
            continue
        if "did not converge" in line:
            continue
        if "on boundry" in line.lower():
            continue
        if "imfit does not work" in line:
            continue

        parts = [part.strip() for part in re.split(r"\s{2,}", line) if part.strip()]

        # Bare source-name marker lines appear between data blocks.
        if len(parts) == 1:
            continue

        if len(parts) < 8:
            continue

        if len(parts) > 8:
            parts = parts[:7] + [" ".join(parts[7:])]

        rows.append(parts[:8])
    return rows


def parse_bash_array(script_text: str, name: str) -> list[str]:
    match = re.search(rf"{re.escape(name)}=\(([^)]*)\)", script_text)
    if not match:
        return []
    return re.findall(r"[^\s()]+", match.group(1))


def get_target_prefixes(base_dir: Path) -> set[str]:
    script_path = base_dir / "run_continuum_fits_clean.sh"
    if not script_path.exists():
        return set()

    script_text = script_path.read_text()
    dates = parse_bash_array(script_text, "datearray")
    flags = parse_bash_array(script_text, "flagarray")
    if not dates or not flags:
        return set()

    return {f"{date}_{flag}" for date in dates for flag in flags}


def collect_summary_files(summary_dirs: list[Path]) -> tuple[dict[str, Path], dict[str, Path]]:
    sorted_files: dict[str, Path] = {}
    witherrors_files: dict[str, Path] = {}

    for directory in summary_dirs:
        sorted_dir = directory / "fitsumfiles" / "sorted"
        witherrors_dir = directory / "fitsumfiles" / "witherrors"

        if sorted_dir.exists():
            for path in sorted_dir.glob(f"*{SORTED_SUFFIX}"):
                if path.name == OUTPUT_NAME:
                    continue
                sorted_files.setdefault(summary_prefix(path, SORTED_SUFFIX), path)

        if witherrors_dir.exists():
            for path in witherrors_dir.glob(f"*{WITHERRORS_SUFFIX}"):
                witherrors_files.setdefault(summary_prefix(path, WITHERRORS_SUFFIX), path)

    return sorted_files, witherrors_files


def collect_rows(base_dir: Path) -> tuple[list[list[str]], list[str]]:
    summary_dirs = [base_dir]

    moved_output_dir = (
        base_dir.parent / "sma_calibration" / "nitercal_0"
        if base_dir.name == "multi_freq_from_archive"
        else None
    )
    if moved_output_dir is not None and moved_output_dir.exists():
        summary_dirs.append(moved_output_dir)

    sorted_files, witherrors_files = collect_summary_files(summary_dirs)

    target_prefixes = get_target_prefixes(base_dir)
    prefixes = sorted(set(sorted_files) | set(witherrors_files))
    if target_prefixes:
        prefixes = [prefix for prefix in prefixes if prefix in target_prefixes]

    rows: list[list[str]] = []
    sources_used: list[str] = []

    for prefix in prefixes:
        if prefix in sorted_files:
            rows.extend(parse_sorted_file(sorted_files[prefix]))
            sources_used.append(f"sorted:{sorted_files[prefix].name}")
        elif prefix in witherrors_files:
            rows.extend(parse_witherrors_file(witherrors_files[prefix]))
            sources_used.append(f"witherrors:{witherrors_files[prefix].name}")

    return rows, sources_used


def deduplicate_rows(rows: list[list[str]]) -> list[list[str]]:
    seen: set[tuple[str, ...]] = set()
    unique_rows: list[list[str]] = []
    for row in rows:
        clean = [str(item).rstrip("\n") for item in row[:8]]
        if (clean[0], clean[1], clean[4]) in EXCLUDED_ROWS:
            continue
        key = tuple(clean)
        if key in seen:
            continue
        seen.add(key)
        unique_rows.append(clean)
    return unique_rows


def write_combined(rows: list[list[str]], output_path: Path) -> None:
    grouped: dict[str, list[list[str]]] = defaultdict(list)
    for row in deduplicate_rows(rows):
        clean = [str(item) for item in row[:8]]
        if not clean[7].endswith("\n"):
            clean[7] = clean[7] + "\n"
        grouped[clean[0]].append(clean)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w") as handle:
        for source in sorted(grouped):
            source_rows = sorted(
                grouped[source],
                key=lambda row: (
                    float_key(row[2]),
                    date_key(row[4]),
                    float_key(row[5]),
                ),
            )
            for row in source_rows:
                handle.write(f"{row!r}\n")
            handle.write("\n\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base-dir",
        default=Path(__file__).resolve().parent,
        type=Path,
        help="Directory containing fitsumfiles/",
    )
    parser.add_argument(
        "--output",
        default=None,
        type=Path,
        help="Override output path. Defaults to fitsumfiles/sorted/fitsummary_fixtochosen_sorted.txt under --base-dir.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only report which files would be combined.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    base_dir = args.base_dir.resolve()
    output_path = (
        args.output.resolve()
        if args.output is not None
        else base_dir / "fitsumfiles" / "sorted" / OUTPUT_NAME
    )

    rows, sources_used = collect_rows(base_dir)
    if not rows:
        parser.error("No fixed-to-chosen summary files were found to combine.")

    print("Combining inputs:")
    for source in sources_used:
        print(f"  - {source}")
    print(f"Rows collected: {len(rows)}")

    if args.dry_run:
        return 0

    write_combined(rows, output_path)
    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
