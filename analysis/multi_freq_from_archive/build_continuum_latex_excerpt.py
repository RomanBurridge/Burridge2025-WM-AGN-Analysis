#!/usr/bin/env python3
"""Build the ten-row printed excerpt of the continuum machine table.

The full machine-readable table retains the canonical numerical fields.  This
builder selects the first ten Circinus rows used in the manuscript and applies
an explicit compact display policy, preventing the printed excerpt from
drifting when a canonical measurement is corrected.
"""

from __future__ import annotations

import argparse
import math
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

from s230_reference import load_fitsummary_records


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = SCRIPT_DIR / "fitsumfiles/final/fitsummary_final_withextras.txt"
DEFAULT_OUTPUT = (
    SCRIPT_DIR.parents[2]
    / "revision_package"
    / "working_clean"
    / "Continuum_from_Pipeline.tex"
)
EXCERPT_SOURCE = "Circinus"
EXCERPT_ROWS = 10
EXPECTED_SOURCE_ROWS = 53


def _finite_decimal(value: object) -> Decimal:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"cannot display a non-finite number: {value!r}")
    return Decimal(str(value))


def format_decimal_places(value: object, decimal_places: int) -> str:
    """Round half up at an explicit fixed-point decimal place."""

    number = _finite_decimal(value)
    quantum = Decimal("1").scaleb(-decimal_places)
    rounded = number.quantize(quantum, rounding=ROUND_HALF_UP)
    return f"{rounded:.{max(decimal_places, 0)}f}"


def format_significant_figures(value: object, significant_figures: int) -> str:
    """Round half up while retaining meaningful trailing zeros."""

    number = _finite_decimal(value)
    if number <= 0:
        raise ValueError(f"significant-figure display requires a positive value: {value!r}")
    decimal_places = significant_figures - 1 - number.adjusted()
    return format_decimal_places(value, decimal_places)


def excerpt_rows(input_path: Path = DEFAULT_INPUT) -> list[dict[str, object]]:
    records = load_fitsummary_records(input_path)
    rows = [row for row in records if row["Source"] == EXCERPT_SOURCE]
    if len(rows) != EXPECTED_SOURCE_ROWS:
        raise ValueError(
            f"expected {EXPECTED_SOURCE_ROWS} {EXCERPT_SOURCE} rows, found {len(rows)}"
        )
    selected = rows[:EXCERPT_ROWS]
    if len(selected) != EXCERPT_ROWS:
        raise ValueError(f"expected {EXCERPT_ROWS} printed rows")
    return selected


def render_table(input_path: Path = DEFAULT_INPUT) -> str:
    lines = [
        r"\begin{deluxetable*}{lllllll}{",
        r"\tablecaption {Continuum Flux-density Measurements}",
        r"\label{tab:mr1}",
        r"\tablehead{ ",
        r"\colhead{[1]} & \colhead{[2]} & \colhead{[3]}  & \colhead{[4]} & \colhead{[5]} & \colhead{[6]}  & \colhead{[7]} }",
        r"}",
        r"\startdata",
    ]
    for row in excerpt_rows(input_path):
        flux = format_decimal_places(row["flux_mjy"], 1)
        if bool(row["is_upper_limit"]):
            flux += "*"
        snr = "NA" if row["snr"] is None else format_decimal_places(row["snr"], 0)
        beam = row["beam_arcsec"]
        if beam is None:
            raise ValueError("printed image-based row lacks a numerical beam")
        date = str(row["date"]).replace("_", r"\_")
        lines.append(
            " & ".join(
                (
                    str(row["Source"]),
                    str(row["telescope"]),
                    format_decimal_places(row["freq_ghz"], 0),
                    flux,
                    date,
                    format_significant_figures(beam, 3),
                    snr,
                )
            )
            + r" \\"
        )
    lines.extend(
        [
            r"\enddata",
            r"\tablecomments{",
            r"[1]: Source Name",
            r"[2]: Telescope",
            r"[3]: Frequency [GHz]",
            r"[4]: Flux Density [mJy]",
            r"[5]: Date [yyyy/mm/dd\_x]",
            r"[6]: Beam Size FWHM [arcsec]",
            r"[7]: SNR",
            r"\\",
            (
                "Ten representative rows are shown; the full 501-row table is provided as "
                "\\texttt{tab2.txt}. It comprises 485 image-based "
                "measurements and 16 literature measurements. Image rows use the beam-shaped "
                "fits described in Section~\\ref{subsec:casa}; literature rows preserve the reported "
                "quantities. An asterisk marks an upper limit, and a date suffix distinguishes "
                "additional products from the same day. Full numerical precision remains in the "
                "machine-readable table.}"
            ),
            r"\end{deluxetable*}",
            "",
        ]
    )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    rendered = render_table(args.input)
    if args.check:
        if not args.output.is_file() or args.output.read_text(encoding="utf-8") != rendered:
            raise RuntimeError(f"stale continuum excerpt: {args.output}")
        print(f"verified {args.output}")
        return 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
