#!/usr/bin/env python3

"""Build the continuum summary tables from the combined fit summary.

This script replaces the fragile notebook-only table assembly logic.  It keeps
the literature extras in one place, checks that every extra is inserted, and
generates both the repo-local fixed-width table and an AAS/CDS MRT version.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


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


EXTRAS = [
    ExtraPoint("NGC4258", "SCUBA-JCMT", "347", "93.7", "29.3", "15", "NA", "Doi et al, 2005, MNRAS, 363:2"),
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
        rows.append([str(item).strip().replace("\\n", "") for item in row[:8]])
    return rows


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
    if "*" in extra.flux_mjy:
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
    if telescope == "SMA" and raw.startswith("/"):
        return "SMA_this_work"
    if raw.startswith(("ALMA/", "VLA/")):
        return raw.split("/")[-1]
    if raw.startswith("/"):
        return Path(raw).name
    return raw


def write_fixed_width_machine(rows: list[list[str]], path: Path) -> None:
    widths = [15, 15, 19, 16, 15, 19, 15, 112]
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
    return f"{float(snr):7.3f}"


def write_aas_mrt(rows: list[list[str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = """Title: When the Shadow Meets Its Measure: Assessing the Feasibility of
       Submillimeter Black Hole Shadow Imaging in Megamaser Disk AGN
Authors: Burridge R.N., Bower G.C.
Table: Continuum Brightness Parameters and Archive Product Identifiers
================================================================================
Byte-by-byte Description of file: datafile2.txt
--------------------------------------------------------------------------------
   Bytes Format Units      Label  Explanations
--------------------------------------------------------------------------------
   1- 12 A12    ---        Name   Source name
  14- 23 A10    ---        Tel    Telescope identifier
  25- 30 F6.2   GHz        Freq   Frequency
  32- 38 F7.2   mJy/beam   Flux   Flux density
  40- 40 A1     ---      f_Flux   [*] Flag on Flux (1)
  42- 53 A12    ---        Date   Observation date; YYYY/MM/DD_N (2)
  55- 62 F8.5   arcsec     FWHM   Beam Full-Width at Half-Maximum size
  64- 64 A1     ---      f_FWHM   [*] Flag on FWHM (3)
  66- 72 A7     ---        SNR    Signal-to-Noise
  74-180 A107   ---        File   File identifier
--------------------------------------------------------------------------------
Note (1):
    * = denotes a non-detection.
Note (2): If the number following the underscore in the date (i.e., the
          \"N\" value) is not 1, it denotes an additional observation taken for
          the same source, at the same frequency, on the same date.
Note (3):
    * = The beam size is 0.19-0.50 arcseconds.
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
                aas_snr(row[6]).ljust(7),
                display_file_label(row).ljust(107),
            ]
            line = " ".join(fields)
            if len(line) != 180:
                raise AssertionError(f"AAS MRT row has length {len(line)}, expected 180: {line!r}")
            handle.write(line + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--sorted-input", type=Path, default=None)
    parser.add_argument("--final-output", type=Path, default=None)
    parser.add_argument("--machine-output", type=Path, default=None)
    parser.add_argument("--aas-output", type=Path, default=None)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    base_dir = args.base_dir.resolve()
    sorted_input = args.sorted_input or base_dir / "fitsumfiles" / "sorted" / "fitsummary_fixtochosen_sorted.txt"
    final_output = args.final_output or base_dir / "fitsumfiles" / "final" / "fitsummary_final_withextras.txt"
    machine_output = args.machine_output or base_dir / "machinetables" / "fitsummary_machine.txt"
    aas_output = args.aas_output or base_dir / "machinetables" / "fitsummary_aas_mrt.txt"

    sorted_rows = parse_python_rows(sorted_input)
    final_rows = build_final_rows(sorted_rows)

    write_final_rows(final_rows, final_output)
    write_fixed_width_machine(final_rows, machine_output)
    write_aas_mrt(final_rows, aas_output)

    print(f"Read {len(sorted_rows)} sorted rows")
    print(f"Inserted {len(EXTRAS)} literature extras")
    print(f"Wrote {len(final_rows)} final rows to {final_output}")
    print(f"Wrote fixed-width table to {machine_output}")
    print(f"Wrote AAS MRT table to {aas_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
