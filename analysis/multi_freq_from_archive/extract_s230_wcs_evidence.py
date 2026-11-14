#!/usr/bin/env python3
"""Freeze the celestial WCS of the 19 selected S230 input images.

Run this script with the CASA Python interpreter.  The generated CSV is a
portable snapshot of the actual selected image metadata; the correction
release can therefore validate the astrometry audit without shipping the
large FITS/CASA image products.  ``--check`` re-reads the images and compares
the resulting bytes with the committed snapshot.
"""

from __future__ import annotations

import argparse
import csv
import io
import math
from pathlib import Path, PurePosixPath
from typing import Mapping, Sequence

from casatools import image


ANALYSIS_DIR = Path(__file__).resolve().parent
S230_TABLE = ANALYSIS_DIR / "machinetables/S230_extrapolation_comparison.csv"
RECOVERY_MANIFEST = ANALYSIS_DIR / "machinetables/continuum_recovery_manifest.csv"
OUTPUT = ANALYSIS_DIR / "s230_wcs_evidence.csv"

FIELDS = (
    "source",
    "selected_observation_id",
    "selected_reference",
    "recovery_record_id",
    "product_sha256",
    "product_sha256_algorithm",
    "wcs_evidence_status",
    "wcs_direction_reference_code",
    "fits_radesys_keyword",
    "fits_equinox_keyword",
    "fits_ctype1",
    "fits_ctype2",
    "wcs_projection",
    "wcs_axis1_name",
    "wcs_axis2_name",
    "wcs_axis1_unit",
    "wcs_axis2_unit",
    "wcs_crpix1_zero_based",
    "wcs_crpix2_zero_based",
    "wcs_crval1_deg",
    "wcs_crval2_deg",
    "wcs_cdelt1_deg_per_pixel",
    "wcs_cdelt2_deg_per_pixel",
    "wcs_pc11",
    "wcs_pc12",
    "wcs_pc21",
    "wcs_pc22",
    "wcs_transform_convention",
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _one(value: object) -> str:
    if isinstance(value, str):
        values = [value]
    else:
        values = list(value)  # type: ignore[arg-type]
    clean = [str(item).strip().upper() for item in values if str(item).strip()]
    if len(clean) != 1:
        raise ValueError(f"expected one direction reference code, found {clean!r}")
    if clean[0] not in {"ICRS", "J2000"}:
        raise ValueError(f"unsupported celestial reference code {clean[0]!r}")
    return clean[0]


def _number(value: object) -> str:
    numeric = float(value)
    if not math.isfinite(numeric):
        raise ValueError(f"non-finite WCS value {value!r}")
    if numeric == 0:
        numeric = 0.0
    return f"{numeric:.17g}"


def _resolve_product(archive_root: Path, reference: str) -> Path:
    relative = PurePosixPath(reference)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"non-portable selected reference: {reference!r}")
    if relative.parts[0] in {"ALMA", "VLA"}:
        candidate = archive_root / "analysis/multi_freq_from_archive" / Path(*relative.parts)
    elif relative.parts[:2] == ("analysis", "sma_calibration"):
        candidate = archive_root / Path(*relative.parts)
    else:
        raise ValueError(f"cannot resolve selected reference: {reference!r}")
    resolved = candidate.resolve()
    root = archive_root.resolve()
    if root not in resolved.parents:
        raise ValueError(f"selected reference escapes archive root: {reference!r}")
    return resolved


def _raw_fits_keywords(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {
            "fits_radesys_keyword": "NA_CASA_IMAGE_TABLE",
            "fits_equinox_keyword": "NA_CASA_IMAGE_TABLE",
            "fits_ctype1": "NA_CASA_IMAGE_TABLE",
            "fits_ctype2": "NA_CASA_IMAGE_TABLE",
        }
    from astropy.io import fits

    with fits.open(path, memmap=True) as hdus:
        header = hdus[0].header
        radesys = header.get("RADESYS", header.get("RADECSYS", "UNAVAILABLE"))
        equinox = header.get("EQUINOX", "UNAVAILABLE")
        return {
            "fits_radesys_keyword": str(radesys).strip() or "UNAVAILABLE",
            "fits_equinox_keyword": str(equinox).strip() or "UNAVAILABLE",
            "fits_ctype1": str(header.get("CTYPE1", "UNAVAILABLE")).strip(),
            "fits_ctype2": str(header.get("CTYPE2", "UNAVAILABLE")).strip(),
        }


def _extract_one(
    selected: Mapping[str, str],
    recovery: Mapping[str, str],
    archive_root: Path,
) -> dict[str, str]:
    reference = selected["nu0_reference"]
    if recovery["local_product_kind"] == "none":
        if recovery["recovery_status"] != "literature_only_no_local_product":
            raise ValueError(
                f"unexpected literature recovery status for {reference}: "
                f"{recovery['recovery_status']!r}"
            )
        sentinel = "NOT_APPLICABLE_LITERATURE_ONLY"
        return {
            "source": selected["Source"],
            "selected_observation_id": selected["nu0_observation_id"],
            "selected_reference": reference,
            "recovery_record_id": recovery["record_id"],
            "product_sha256": "NA",
            "product_sha256_algorithm": "not_applicable",
            "wcs_evidence_status": "literature_only_no_local_product",
            "wcs_direction_reference_code": sentinel,
            "fits_radesys_keyword": sentinel,
            "fits_equinox_keyword": sentinel,
            "fits_ctype1": sentinel,
            "fits_ctype2": sentinel,
            "wcs_projection": sentinel,
            "wcs_axis1_name": sentinel,
            "wcs_axis2_name": sentinel,
            "wcs_axis1_unit": sentinel,
            "wcs_axis2_unit": sentinel,
            "wcs_crpix1_zero_based": sentinel,
            "wcs_crpix2_zero_based": sentinel,
            "wcs_crval1_deg": sentinel,
            "wcs_crval2_deg": sentinel,
            "wcs_cdelt1_deg_per_pixel": sentinel,
            "wcs_cdelt2_deg_per_pixel": sentinel,
            "wcs_pc11": sentinel,
            "wcs_pc12": sentinel,
            "wcs_pc21": sentinel,
            "wcs_pc22": sentinel,
            "wcs_transform_convention": sentinel,
        }
    path = _resolve_product(archive_root, reference)
    if not path.exists():
        raise FileNotFoundError(path)

    tool = image()
    coordinates = None
    try:
        if not tool.open(str(path)):
            raise RuntimeError(f"CASA could not open {path}")
        coordinates = tool.coordsys()
        record = coordinates.torecord()
        code = _one(coordinates.referencecode(type="direction"))
    finally:
        if coordinates is not None:
            try:
                coordinates.done()
            except Exception:
                pass
        try:
            tool.close()
        except Exception:
            pass

    direction_keys = sorted(key for key in record if key.startswith("direction"))
    if len(direction_keys) != 1:
        raise ValueError(f"expected one direction coordinate in {path}")
    direction = record[direction_keys[0]]
    axes = list(direction["axes"])
    units = list(direction["units"])
    crpix = list(direction["crpix"])
    crval = list(direction["crval"])
    cdelt = list(direction["cdelt"])
    pc = direction["pc"]
    if len(axes) != 2 or len(units) != 2:
        raise ValueError(f"unexpected celestial axes in {path}: {axes!r}, {units!r}")
    if list(units) != ["rad", "rad"]:
        raise ValueError(f"unexpected celestial WCS units in {path}: {units!r}")
    if str(direction.get("system", "")).upper() != code:
        raise ValueError(f"contradictory celestial reference metadata in {path}")

    radians_to_degrees = 180.0 / math.pi
    raw = _raw_fits_keywords(path)
    return {
        "source": selected["Source"],
        "selected_observation_id": selected["nu0_observation_id"],
        "selected_reference": reference,
        "recovery_record_id": recovery["record_id"],
        "product_sha256": recovery["sha256"],
        "product_sha256_algorithm": recovery["sha256_algorithm"],
        "wcs_evidence_status": recovery["recovery_status"],
        "wcs_direction_reference_code": code,
        **raw,
        "wcs_projection": str(direction.get("projection", "UNAVAILABLE")),
        "wcs_axis1_name": str(axes[0]),
        "wcs_axis2_name": str(axes[1]),
        "wcs_axis1_unit": str(units[0]),
        "wcs_axis2_unit": str(units[1]),
        "wcs_crpix1_zero_based": _number(crpix[0]),
        "wcs_crpix2_zero_based": _number(crpix[1]),
        "wcs_crval1_deg": _number(float(crval[0]) * radians_to_degrees),
        "wcs_crval2_deg": _number(float(crval[1]) * radians_to_degrees),
        "wcs_cdelt1_deg_per_pixel": _number(float(cdelt[0]) * radians_to_degrees),
        "wcs_cdelt2_deg_per_pixel": _number(float(cdelt[1]) * radians_to_degrees),
        "wcs_pc11": _number(pc[0][0]),
        "wcs_pc12": _number(pc[0][1]),
        "wcs_pc21": _number(pc[1][0]),
        "wcs_pc22": _number(pc[1][1]),
        "wcs_transform_convention": (
            "CASA_direction0_zero_based_crpix;world_rad="
            "crval_rad+PC@(pixel-crpix)*cdelt_rad"
        ),
    }


def build_rows(archive_root: Path) -> list[dict[str, str]]:
    selected_rows = _read_csv(S230_TABLE)
    recovery_rows = _read_csv(RECOVERY_MANIFEST)
    recovery_by_reference = {row["product_or_citation"]: row for row in recovery_rows}
    if len(selected_rows) != 19:
        raise ValueError(f"expected 19 S230 rows, found {len(selected_rows)}")

    rows: list[dict[str, str]] = []
    for selected in selected_rows:
        reference = selected["nu0_reference"]
        if reference not in recovery_by_reference:
            raise ValueError(f"selected product has no recovery record: {reference}")
        recovery = recovery_by_reference[reference]
        if recovery["source"] != selected["Source"]:
            raise ValueError(f"source mismatch for selected product {reference}")
        rows.append(_extract_one(selected, recovery, archive_root))
    return sorted(rows, key=lambda row: row["source"])


def render(rows: Sequence[Mapping[str, str]]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--archive-root",
        type=Path,
        required=True,
        help="local Burridge2025-WM-AGN-Analysis project root",
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    content = render(build_rows(args.archive_root))
    if args.check:
        if not args.output.is_file() or args.output.read_text(encoding="utf-8") != content:
            raise SystemExit(f"stale S230 WCS evidence: {args.output}")
        print(f"S230 WCS evidence is current: 19/19 rows ({args.output})")
        return 0
    args.output.write_text(content, encoding="utf-8")
    print(f"Wrote 19/19 selected-image WCS rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
