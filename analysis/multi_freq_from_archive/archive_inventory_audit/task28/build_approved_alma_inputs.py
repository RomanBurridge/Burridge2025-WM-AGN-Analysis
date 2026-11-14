#!/usr/bin/env python3
"""Materialize and validate the 82 approved Task 28 ALMA inputs.

The source/MOUS inventory is frozen by ``alma_archive_mous_candidate_review.csv``.
For ordinary products this script selects the representative target continuum
image.  For the seven manually reviewed MFS deliveries it selects either the
native target Stokes-I image or the target IQUV image.  An IQUV parent is reduced
to a singleton Stokes-I plane without changing its spatial or frequency axes.

Every output retains its official ALMA URL, MOUS, parent SHA-256, output
SHA-256, FITS identity, and sky-coverage check in a generated validation CSV.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import math
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel


APPROVED_DISPOSITIONS = {
    "pre_cutoff_possible_omission",
    "post_cutoff_possible_addition",
    "no_aggregate_image_manual_spw_review",
}
MANUAL_DISPOSITION = "no_aggregate_image_manual_spw_review"
EXPECTED_COUNTS = {
    "pre_cutoff_possible_omission": 49,
    "post_cutoff_possible_addition": 26,
    "no_aggregate_image_manual_spw_review": 7,
}
EXPECTED_TOTAL = 82
# CHECKSUM/DATASUM comments are serialization metadata, not observing or
# processing timestamps.  Preserve the five reviewed products' original
# comment templates so an independent rebuild satisfies the unchanged
# whole-file admission hashes.  Keys are exact official-parent SHA-256s.
# Unknown parents use a fixed format epoch; they are NOT thereby admitted.
# The downstream product policy still requires the complete output SHA-256.
STOKES_I_CHECKSUM_FORMAT_EPOCH = "2026-08-26T00:00:00"
REVIEWED_STOKES_I_CHECKSUM_TIMES = {
    "9eaae5e87e0eac00c34ca90c14041129a56eda168ec16ca5398a48334257af72": "2026-08-25T19:49:29",
    "59440f136e1a2a151b06a8bf87dc5720c5d9a7ddfefdb49b22428e17cf31f777": "2026-08-25T19:49:29",
    "9f86b506dca89adc6967ec19c64b1242baa38443cbbf10a81ea38317f28e3ef2": "2026-08-25T19:55:02",
    "93c43a5e5d7092d394626ffef0cb08242f1ec2a9f5670c133c4ee3ef3c13b504": "2026-08-25T19:49:29",
    "56376e9e31e52f91ea485d815d14e4322e68d1f6fc31d30110a56f45731ab73a": "2026-08-25T19:49:30",
}
# The five 2014 multi-source deliveries use anonymous ``sourceN`` filenames.
# The archive-level target metadata cannot identify the matching image.  These
# exact choices were therefore resolved from each FITS OBJECT card and WCS:
# every selected image is centered within one arcsecond of NGC 4388.
SPLIT_CAL_SOURCE_OVERRIDES = {
    "uid://A002/X5a9a13/Xe5": "uid___A002_X828760_X533.ms.split.cal.image.continuum.source6.pbcor.fits",
    "uid://A002/X5a9a13/Xe9": "uid___A002_X828760_X398.ms.split.cal.image.continuum.source12.pbcor.fits",
    "uid://A002/X5a9a13/Xed": "uid___A002_X828760_X207.ms.split.cal.image.continuum.source5.pbcor.fits",
    "uid://A002/X5a9a13/Xf1": "uid___A002_X828760_X2c2.ms.split.cal.image.continuum.source9.pbcor.fits",
    "uid://A002/X5a9a13/Xf5": "uid___A002_X828760_X454.ms.split.cal.image.continuum.source13.pbcor.fits",
}
FIELDS = (
    "source",
    "member_ous_uid",
    "screening_disposition",
    "proposal_ids",
    "parent_product_basename",
    "parent_product_url",
    "output_product_basename",
    "output_relative_path",
    "transformation",
    "parent_sha256",
    "output_sha256",
    "output_size_bytes",
    "header_object",
    "header_telescope",
    "header_bunit",
    "header_date_obs",
    "header_frequency_ghz",
    "header_shape",
    "stokes_axis_length",
    "stokes_i_world_value",
    "source_inside_image",
    "source_to_image_center_arcsec",
    "validation_status",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def write_csv(path: Path, rows: list[dict[str, str]], fields: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def select_rows(review: list[dict[str, str]], products: list[dict[str, str]]) -> list[dict[str, str]]:
    product_urls = {
        (row["member_ous_uid"], row["product_basename"]): row["product_url"]
        for row in products
    }
    selected: list[dict[str, str]] = []
    for row in review:
        disposition = row["screening_disposition"]
        if disposition not in APPROVED_DISPOSITIONS:
            continue
        if disposition != MANUAL_DISPOSITION:
            basename = SPLIT_CAL_SOURCE_OVERRIDES.get(
                row["member_ous_uid"], row["representative_product_basename"]
            )
            url = product_urls.get(
                (row["member_ous_uid"], basename), row["representative_product_url"]
            )
            transformation = "none"
        else:
            eligible = row["eligible_product_basenames"].split(";")
            native_i = [
                name
                for name in eligible
                if name.endswith(".mfs.I.manual.image.pbcor.fits")
            ]
            target_iquv = [
                name
                for name in eligible
                if name.endswith(".mfs.IQUV.manual.pbcor.fits")
            ]
            if len(native_i) == 1:
                basename = native_i[0]
                transformation = "none_native_stokes_i"
            elif len(target_iquv) == 1:
                basename = target_iquv[0]
                transformation = "extract_singleton_stokes_i_from_iquv"
            else:
                raise RuntimeError(
                    f"manual MOUS does not have one approved target intensity product: "
                    f"{row['member_ous_uid']}"
                )
            try:
                url = product_urls[(row["member_ous_uid"], basename)]
            except KeyError as exc:
                raise RuntimeError(f"missing audited URL for {basename}") from exc
        if not basename or not url:
            raise RuntimeError(f"selected product is incomplete for {row['member_ous_uid']}")
        output_basename = basename
        if transformation == "extract_singleton_stokes_i_from_iquv":
            output_basename = basename.replace(
                ".mfs.IQUV.manual.pbcor.fits",
                ".mfs.I.stokesI.manual.pbcor.fits",
            )
        selected.append(
            {
                "source": row["source"],
                "member_ous_uid": row["member_ous_uid"],
                "screening_disposition": disposition,
                "proposal_ids": row["proposal_ids"],
                "parent_product_basename": basename,
                "parent_product_url": url,
                "output_product_basename": output_basename,
                "output_relative_path": f"ALMA/{row['source']}/{output_basename}",
                "transformation": transformation,
            }
        )
    counts = {key: 0 for key in EXPECTED_COUNTS}
    for row in selected:
        counts[row["screening_disposition"]] += 1
    if counts != EXPECTED_COUNTS or len(selected) != EXPECTED_TOTAL:
        raise RuntimeError(f"approved inventory changed: counts={counts}, total={len(selected)}")
    if len({row["member_ous_uid"] for row in selected}) != EXPECTED_TOTAL:
        raise RuntimeError("approved inventory contains duplicate MOUS identifiers")
    if len({row["output_relative_path"] for row in selected}) != EXPECTED_TOTAL:
        raise RuntimeError("approved inventory contains duplicate output paths")
    return sorted(selected, key=lambda row: row["output_relative_path"])


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    command = [
        "curl",
        "--fail",
        "--location",
        "--retry",
        "8",
        "--retry-all-errors",
        "--connect-timeout",
        "20",
        "--continue-at",
        "-",
        "--silent",
        "--show-error",
        "--output",
        str(partial),
        url,
    ]
    subprocess.run(command, check=True)
    partial.replace(destination)


def materialize_parent(
    row: dict[str, str],
    parent_root: Path,
    output_root: Path,
    reusable: dict[str, Path],
) -> Path:
    output = output_root / row["output_relative_path"]
    if row["transformation"].startswith("extract_"):
        parent = parent_root / row["source"] / row["parent_product_basename"]
    else:
        parent = output
    if parent.is_file():
        return parent
    reuse = reusable.get(row["parent_product_url"])
    if reuse is not None:
        parent.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(reuse, parent)
        return parent
    download(row["parent_product_url"], parent)
    return parent


def stokes_axis(header: fits.Header) -> int | None:
    for axis in range(1, int(header.get("NAXIS", 0)) + 1):
        if str(header.get(f"CTYPE{axis}", "")).strip().upper().startswith("STOKES"):
            return axis
    return None


def stokes_world_value(header: fits.Header, axis: int, index_zero_based: int) -> float:
    return float(header[f"CRVAL{axis}"]) + (
        index_zero_based + 1 - float(header.get(f"CRPIX{axis}", 1.0))
    ) * float(header.get(f"CDELT{axis}", 1.0))


def extract_stokes_i(parent: Path, output: Path) -> None:
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Stokes-I output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    parent_sha256 = sha256(parent)
    with fits.open(parent, memmap=True) as hdul:
        data = np.asarray(hdul[0].data)
        header = hdul[0].header.copy()
        axis = stokes_axis(header)
        if axis is None:
            raise RuntimeError(f"IQUV parent has no Stokes axis: {parent}")
        numpy_axis = data.ndim - axis
        values = [stokes_world_value(header, axis, index) for index in range(data.shape[numpy_axis])]
        matches = [index for index, value in enumerate(values) if math.isclose(value, 1.0, abs_tol=1e-8)]
        if len(matches) != 1:
            raise RuntimeError(f"IQUV parent lacks one Stokes-I plane: {parent}; values={values}")
        index = matches[0]
        slices = [slice(None)] * data.ndim
        slices[numpy_axis] = slice(index, index + 1)
        selected = np.array(data[tuple(slices)], copy=True)
        header[f"CRPIX{axis}"] = float(header.get(f"CRPIX{axis}", 1.0)) - index
        # Some CASA polarization exports store one beam per polarization in a
        # BEAMS extension instead of the primary BMAJ/BMIN/BPA cards.  Once a
        # single polarization is selected, promote that exact beam to the
        # primary header so ordinary CASA image fitting sees one unambiguous
        # restoring beam.
        if "BMAJ" not in header:
            beams = next((hdu for hdu in hdul[1:] if hdu.name == "BEAMS"), None)
            if beams is None or beams.data is None:
                raise RuntimeError(f"IQUV parent has no restoring beam: {parent}")
            matches = [
                beam
                for beam in beams.data
                if int(beam["CHAN"]) == 0 and int(beam["POL"]) == index
            ]
            if len(matches) != 1:
                raise RuntimeError(
                    f"IQUV parent has no unique Stokes-I restoring beam: {parent}"
                )
            beam = matches[0]
            units = {column.name: str(column.unit).strip().casefold() for column in beams.columns}
            major = float(beam["BMAJ"])
            minor = float(beam["BMIN"])
            if units["BMAJ"] == "arcsec":
                major /= 3600.0
            elif units["BMAJ"] != "deg":
                raise RuntimeError(f"unsupported BMAJ unit in {parent}: {units['BMAJ']}")
            if units["BMIN"] == "arcsec":
                minor /= 3600.0
            elif units["BMIN"] != "deg":
                raise RuntimeError(f"unsupported BMIN unit in {parent}: {units['BMIN']}")
            if units["BPA"] != "deg":
                raise RuntimeError(f"unsupported BPA unit in {parent}: {units['BPA']}")
            header["BMAJ"] = major
            header["BMIN"] = minor
            header["BPA"] = float(beam["BPA"])
        header.add_history("Task 28: retained only the Stokes-I plane from official ALMA IQUV parent")
        header.add_history(f"Parent SHA256: {parent_sha256}")
        hdu = fits.PrimaryHDU(data=selected, header=header)
        stamp = REVIEWED_STOKES_I_CHECKSUM_TIMES.get(
            parent_sha256, STOKES_I_CHECKSUM_FORMAT_EPOCH
        )
        hdu.add_datasum(when=f"data unit checksum updated {stamp}")
        hdu.add_checksum(when=f"HDU checksum updated {stamp}", override_datasum=True)
        # checksum=True would replace the deterministic comments with the
        # current wall clock.  Both real checksums were calculated above.
        hdu.writeto(output, overwrite=False, checksum=False)


def source_coordinates(path: Path) -> dict[str, SkyCoord]:
    raw = ast.literal_eval(path.read_text(encoding="utf-8"))
    return {
        item[0][0]: SkyCoord(item[1][1], item[1][2], unit=(u.hourangle, u.deg), frame="icrs")
        for item in raw
    }


def header_frequency_ghz(header: fits.Header) -> float:
    factors = {"hz": 1e-9, "khz": 1e-6, "mhz": 1e-3, "ghz": 1.0}
    for axis in range(1, int(header.get("NAXIS", 0)) + 1):
        if str(header.get(f"CTYPE{axis}", "")).strip().upper().startswith("FREQ"):
            unit = str(header.get(f"CUNIT{axis}", "Hz")).strip().lower()
            return float(header[f"CRVAL{axis}"]) * factors[unit]
    raise RuntimeError("FITS header has no frequency axis")


def validate(
    row: dict[str, str],
    parent: Path,
    output: Path,
    coordinates: dict[str, SkyCoord],
) -> dict[str, str]:
    with fits.open(output, memmap=True, checksum=True) as hdul:
        data = hdul[0].data
        header = hdul[0].header.copy()
        if data is None or data.ndim < 2:
            raise RuntimeError(f"invalid FITS image array: {output}")
        shape = data.shape
    celestial = WCS(header).celestial
    source = coordinates[row["source"]]
    x, y = skycoord_to_pixel(source, celestial, origin=0)
    inside = bool(np.isfinite(x) and np.isfinite(y) and 0 <= x < shape[-1] and 0 <= y < shape[-2])
    center = celestial.pixel_to_world((shape[-1] - 1) / 2, (shape[-2] - 1) / 2)
    axis = stokes_axis(header)
    axis_length = ""
    i_value = ""
    stokes_ok = True
    if axis is not None:
        numpy_axis = len(shape) - axis
        axis_length = str(shape[numpy_axis])
        values = [stokes_world_value(header, axis, index) for index in range(shape[numpy_axis])]
        i_values = [value for value in values if math.isclose(value, 1.0, abs_tol=1e-8)]
        i_value = f"{i_values[0]:.1f}" if len(i_values) == 1 else ""
        if row["transformation"].startswith("extract_"):
            stokes_ok = shape[numpy_axis] == 1 and len(i_values) == 1
    valid = (
        inside
        and str(header.get("TELESCOP", "")).strip().upper() == "ALMA"
        and str(header.get("BUNIT", "")).strip().casefold() == "jy/beam"
        and bool(str(header.get("DATE-OBS", "")))
        and stokes_ok
    )
    return {
        **row,
        "parent_sha256": sha256(parent),
        "output_sha256": sha256(output),
        "output_size_bytes": str(output.stat().st_size),
        "header_object": str(header.get("OBJECT", "")),
        "header_telescope": str(header.get("TELESCOP", "")),
        "header_bunit": str(header.get("BUNIT", "")),
        "header_date_obs": str(header.get("DATE-OBS", "")),
        "header_frequency_ghz": f"{header_frequency_ghz(header):.10g}",
        "header_shape": "x".join(str(value) for value in shape),
        "stokes_axis_length": axis_length,
        "stokes_i_world_value": i_value,
        "source_inside_image": "yes" if inside else "no",
        "source_to_image_center_arcsec": f"{source.separation(center).arcsec:.10g}",
        "validation_status": "validated" if valid else "review",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit-root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument(
        "--names-and-coordinates",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "namesandcoords.txt",
    )
    parser.add_argument("--workers", type=int, default=6)
    arguments = parser.parse_args()
    root = arguments.audit_root.resolve()
    review = read_csv(root / "alma_archive_mous_candidate_review.csv")
    products = read_csv(root / "alma_archive_fits_product_audit.csv")
    selected = select_rows(review, products)
    manifest_fields = tuple(key for key in FIELDS if key in selected[0])
    write_csv(root / "approved_alma_input_manifest.csv", selected, manifest_fields)

    reusable = {
        row["official_product_url"]: root / row["local_path"]
        for row in read_csv(root / "task28_downloaded_alma_provenance_audit.csv")
        if row["provenance_status"] == "unique_exact_size_match_and_header_valid"
    }
    parent_root = root / "approved_alma_parent_products"
    output_root = root / "approved_alma_products"
    parents: dict[str, Path] = {}
    with ThreadPoolExecutor(max_workers=arguments.workers) as executor:
        futures = {
            executor.submit(materialize_parent, row, parent_root, output_root, reusable): row
            for row in selected
        }
        for index, future in enumerate(as_completed(futures), start=1):
            row = futures[future]
            parents[row["output_relative_path"]] = future.result()
            print(f"materialized {index}/{len(selected)} {row['output_relative_path']}", flush=True)

    for row in selected:
        output = output_root / row["output_relative_path"]
        parent = parents[row["output_relative_path"]]
        if row["transformation"].startswith("extract_") and not output.is_file():
            extract_stokes_i(parent, output)

    coordinates = source_coordinates(arguments.names_and_coordinates.resolve())
    validated = [
        validate(
            row,
            parents[row["output_relative_path"]],
            output_root / row["output_relative_path"],
            coordinates,
        )
        for row in selected
    ]
    write_csv(root / "approved_alma_products_validated.csv", validated, FIELDS)
    review_rows = [row for row in validated if row["validation_status"] != "validated"]
    if review_rows:
        raise RuntimeError(
            "approved output validation failed: "
            + ", ".join(row["output_relative_path"] for row in review_rows)
        )
    derived = sum(row["transformation"].startswith("extract_") for row in validated)
    print(
        f"validated outputs={len(validated)} derived_stokes_i={derived} "
        f"bytes={sum(int(row['output_size_bytes']) for row in validated)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
