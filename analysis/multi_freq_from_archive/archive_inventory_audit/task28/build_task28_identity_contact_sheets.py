#!/usr/bin/env python3
"""Render direct-review cutouts for Task 28 admitted ALMA products."""

from __future__ import annotations

import argparse
import ast
import csv
import math
import sys
from pathlib import Path

from astropy.coordinates import FK5, ICRS, SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPT_DIR.parents[1]
for path in (str(ANALYSIS_ROOT), str(SCRIPT_DIR)):
    if path not in sys.path:
        sys.path.insert(0, path)

from build_task28_identity_policy import SOURCE_REFERENCES, native_reference  # noqa: E402
from continuum_classification_identity import (  # noqa: E402
    load_ngc4945_source18_identity_records,
)
from continuum_product_additions import (  # noqa: E402
    load_product_addition_policy,
    publication_addition_products,
)


AUDIT_FIELDS = (
    "page",
    "panel",
    "source",
    "frequency_ghz",
    "product",
    "classification",
    "flux_mjy",
    "snr",
    "beam_major_arcsec",
    "fit_reference_separation_arcsec",
    "fit_reference_separation_bmaj",
    "review_status",
)


def read_sorted(path: Path) -> list[list[object]]:
    rows: list[list[object]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        row = ast.literal_eval(line)
        row[7] = str(row[7]).strip().replace("\\", "/")
        rows.append(row)
    return rows


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def is_detection(row: list[object]) -> bool:
    return not str(row[3]).endswith("*") and str(row[6]) != "NA"


def coordinate(frame: str, ra: float, dec: float) -> SkyCoord:
    return SkyCoord(
        ra=ra * u.deg,
        dec=dec * u.deg,
        frame=ICRS() if frame == "ICRS" else FK5(equinox=Time("J2000")),
    )


def image_plane(path: Path) -> tuple[np.ndarray, fits.Header, WCS]:
    with fits.open(path, memmap=True) as hdul:
        header = hdul[0].header.copy()
        data = np.asarray(hdul[0].data)
    data = np.squeeze(data)
    if data.ndim != 2:
        raise ValueError(f"{path}: expected a two-dimensional Stokes-I plane")
    unit = str(header.get("BUNIT", "")).strip().casefold().replace(" ", "")
    if unit in {"jy/beam", "jybeam-1", "jybeam^-1"}:
        data = data.astype(float, copy=False) * 1000.0
    elif unit not in {"mjy/beam", "mjybeam-1", "mjybeam^-1"}:
        raise ValueError(f"{path}: unsupported brightness unit {header.get('BUNIT')!r}")
    return data, header, WCS(header).celestial


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--generated-table", type=Path, required=True)
    parser.add_argument("--noise-audit", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--review-confirmed",
        action="store_true",
        help="record that every rendered panel received direct visual review",
    )
    parser.add_argument(
        "--include-upper-limits",
        action="store_true",
        help=(
            "render every admitted product for the broader image-quality audit; "
            "the default identity-review sheet contains detections only"
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.analysis_root.resolve()
    additions = load_product_addition_policy(
        root / "continuum_product_addition_policy.csv", analysis_root=root
    )
    admitted = publication_addition_products(additions)
    generated = {
        str(row[7]): row
        for row in read_sorted(args.generated_table)
        if str(row[7]) in admitted
        and (args.include_upper_limits or is_detection(row))
    }
    audits = {
        row["product"].strip().replace("\\", "/"): row
        for row in read_csv(args.noise_audit)
    }
    audits.update(
        load_ngc4945_source18_identity_records(
            root / "ngc4945_nuclear_photometry.csv"
        )
    )
    products = sorted(generated, key=lambda item: (str(generated[item][0]), item))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_stem = (
        "task28_all_additions_contact_sheets"
        if args.include_upper_limits
        else "task28_identity_contact_sheets"
    )
    pdf_path = args.output_dir / f"{output_stem}.pdf"
    audit_rows: list[dict[str, str]] = []

    panels_per_page = 12
    with PdfPages(pdf_path) as pdf:
        for page_index, start in enumerate(range(0, len(products), panels_per_page), 1):
            page_products = products[start : start + panels_per_page]
            figure, axes = plt.subplots(4, 3, figsize=(12, 15), constrained_layout=True)
            flat_axes = list(axes.flat)
            for panel_index, (axis, product) in enumerate(
                zip(flat_axes, page_products), 1
            ):
                row = generated[product]
                record = audits[product]
                source = str(row[0])
                frame = record["fit_center_reference_frame"]
                fit = coordinate(
                    frame,
                    float(record["fit_center_ra_deg"]),
                    float(record["fit_center_dec_deg"]),
                )
                reference_ra, reference_dec = native_reference(
                    SOURCE_REFERENCES[source], frame
                )
                reference = coordinate(frame, reference_ra, reference_dec)
                separation = float(fit.separation(reference).arcsec)
                beam_major = float(record["beam_major_arcsec"])
                data, header, wcs = image_plane(root / product)
                fit_x, fit_y = wcs.world_to_pixel(fit)
                ref_x, ref_y = wcs.world_to_pixel(reference)
                pixel_width = float(record["pixel_width_arcsec"])
                pixel_height = float(record["pixel_height_arcsec"])
                half_x = max(
                    20,
                    int(math.ceil(float(record["aperture_half_width_x_pixels"]) * 1.2)),
                )
                half_y = max(
                    20,
                    int(math.ceil(float(record["aperture_half_width_y_pixels"]) * 1.2)),
                )
                left = max(0, int(math.floor(fit_x)) - half_x)
                right = min(data.shape[1], int(math.ceil(fit_x)) + half_x + 1)
                bottom = max(0, int(math.floor(fit_y)) - half_y)
                top = min(data.shape[0], int(math.ceil(fit_y)) + half_y + 1)
                cutout = data[bottom:top, left:right]
                finite = cutout[np.isfinite(cutout)]
                norm = ImageNormalize(
                    finite,
                    interval=PercentileInterval(99.0),
                    stretch=AsinhStretch(),
                    clip=True,
                )
                axis.imshow(cutout, origin="lower", cmap="viridis", norm=norm)
                axis.plot(fit_x - left, fit_y - bottom, "rx", ms=9, mew=2)
                axis.plot(ref_x - left, ref_y - bottom, "w+", ms=9, mew=2)
                radius_pixels = 0.5 * beam_major / math.sqrt(pixel_width * pixel_height)
                axis.add_patch(
                    Circle(
                        (8 + radius_pixels, 8 + radius_pixels),
                        radius_pixels,
                        edgecolor="white",
                        facecolor="none",
                        linewidth=1.2,
                    )
                )
                axis.set_title(
                    f"{source}  {float(row[2]):.2f} GHz\n"
                    f"{row[3]} mJy; S/N {row[6]}; "
                    f"{separation / beam_major:.2f} Bmaj",
                    fontsize=8,
                )
                axis.set_xticks([])
                axis.set_yticks([])
                audit_rows.append(
                    {
                        "page": str(page_index),
                        "panel": str(panel_index),
                        "source": source,
                        "frequency_ghz": str(row[2]),
                        "product": product,
                        "classification": record["classification"],
                        "flux_mjy": record["flux_mjy"],
                        "snr": record["snr"],
                        "beam_major_arcsec": record["beam_major_arcsec"],
                        "fit_reference_separation_arcsec": f"{separation:.12g}",
                        "fit_reference_separation_bmaj": f"{separation / beam_major:.12g}",
                        "review_status": (
                            (
                                "accepted_direct_visual_quality_review"
                                if args.include_upper_limits
                                else "accepted_direct_review"
                            )
                            if args.review_confirmed
                            else "pending_direct_review"
                        ),
                    }
                )
            for axis in flat_axes[len(page_products) :]:
                axis.axis("off")
            figure.suptitle(
                (
                    "Task 28 admitted products: red x = final fit; "
                    "white + = reviewed reference"
                    if args.include_upper_limits
                    else "Task 28 admitted detections: red x = final fit; "
                    "white + = reviewed reference"
                ),
                fontsize=12,
            )
            png_stem = (
                "task28_all_additions_contact_sheet"
                if args.include_upper_limits
                else "task28_identity_contact_sheet"
            )
            png_path = args.output_dir / f"{png_stem}_{page_index:02d}.png"
            figure.savefig(png_path, dpi=170)
            pdf.savefig(figure)
            plt.close(figure)

    audit_name = (
        "task28_all_additions_visual_review.csv"
        if args.include_upper_limits
        else "task28_identity_visual_review.csv"
    )
    audit_path = args.output_dir / audit_name
    with audit_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(audit_rows)
    label = "admitted products" if args.include_upper_limits else "detections"
    print(f"rendered {len(products)} {label} on {math.ceil(len(products) / 12)} pages")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
