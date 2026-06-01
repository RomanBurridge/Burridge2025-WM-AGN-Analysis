#!/usr/bin/env python3
"""Build the legacy Greene_params.txt file from BHS_params.csv.

This file is kept only as a compatibility bridge for any older notebook that
still expects ``analysis/CountingBHS_and_Spectra/Greene_params/Greene_params.txt``.
The source of truth is now ``BH_Parameters.ipynb``, which exports
``analysis/CountingBHS_and_Spectra/BHS_params.csv``.
"""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "analysis" / "CountingBHS_and_Spectra" / "BHS_params.csv"
OUTPUT = (
    ROOT
    / "analysis"
    / "CountingBHS_and_Spectra"
    / "Greene_params"
    / "Greene_params.txt"
)


def _format_value(value: str) -> str:
    if value in ("", "NA"):
        return "NA"
    return f"{float(value):.5f}"


def build_rows() -> list[str]:
    with INPUT.open(newline="") as fh:
        reader = csv.DictReader(fh)
        rows = []
        for row in reader:
            source = row["source_key"]
            log_mass = _format_value(row["log_mass"])
            log_mass_unc = _format_value(row["log_mass_unc"])
            rows.append(f"{source} x x {log_mass} x {log_mass_unc}")
    return rows


def main() -> None:
    if not INPUT.exists():
        raise FileNotFoundError(f"BHS parameter table not found: {INPUT}")

    rows = build_rows()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text("\n".join(rows) + "\n")
    print(f"Wrote {len(rows)} rows to {OUTPUT}")


if __name__ == "__main__":
    main()
