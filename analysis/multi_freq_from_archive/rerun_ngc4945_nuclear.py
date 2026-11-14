#!/usr/bin/env python3
"""Reproduce and apply the fixed-nucleus NGC 4945 continuum correction.

The historical continuum fitter selected Emig et al. (2020) source 17 as the
smallest-beam 349.4 GHz anchor and then used that position to seed every other
NGC 4945 fit.  Source 18, not source 17, is the compact source identified with
the AGN/maser nucleus.  This script deliberately performs a source-specific
correction instead of allowing another brightest-pixel search:

* transform the published source-18 ICRS coordinate through each FITS WCS;
* fix the Gaussian centroid, restoring-beam major/minor axes, and position
  angle (CASA estimates flag ``xyabp``), leaving only amplitude free;
* use the same 3-sigma/FWHM square fit-box convention as the released fitter;
* calculate classification S/N from the scaled MAD in the final 10-BMAJ,
  source-excluded local aperture used for generic publication fits, while
  retaining classic RMS, sigma, and mean as diagnostics; and
* preserve the rejected source-17 anchor, a free-centroid validation fit, and
  wider-region sensitivity checks in machine-readable audit sidecars.

CASA is needed only to regenerate the photometry.  ``--check`` validates the
checked-in records and canonical sorted table with the standard interpreter.
The archive root is supplied explicitly so no author-local absolute path is
written into a release product.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import math
import tempfile
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from continuum_measurements import (
    FINAL_NOISE_APERTURE_BMAJ_MULTIPLIER,
    FINAL_NOISE_MIN_INDEPENDENT_BEAMS,
    SCALED_MAD_NORMAL_FACTOR,
    THREE_SIGMA_FROM_FWHM,
    classify_ordinary_fit,
    effective_independent_noise_samples,
    final_noise_review_flags,
    gaussian_beam_area_arcsec2,
    independent_beam_areas,
    require_fixed_restoring_beam_shape,
    require_minimum_independent_beams,
    square_aperture_coverage,
)
from continuum_units import brightness_to_mjy_per_beam_factor


SOURCE = "NGC4945"
FIXED_RA_ICRS = "13:05:27.4798"
FIXED_DEC_ICRS = "-49:28:05.404"
COORDINATE_REFERENCE = "Emig2020_source18_doi:10.3847/1538-4357/abb67d"
COORDINATE_FRAME = "ICRS"
EXPECTED_ROWS = 45
EXPECTED_LEGACY_ROWS = 27
FIT_BOX_FACTOR = 3.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))
EXPECTED_INPUT_HASH_MANIFEST_SHA256 = (
    "545e4b39694fb3be57a5a28c35857731a5e090d2285c67a5792ad816c16f8f8a"
)
EXPECTED_COMPONENT_RECORD_IDS = frozenset(
    {
        "rejected-source17-archived-anchor",
        "accepted-source18-fixed-canonical",
        "source18-free-centroid-check",
        "source18-fixed-width-0p18",
        "source18-free-width-0p18",
        "source18-fixed-width-0p24",
        "source18-free-width-0p24",
    }
)

PACKAGE_DIR = Path(__file__).resolve().parent
DEFAULT_SORTED = PACKAGE_DIR / "fitsumfiles/sorted/fitsummary_fixtochosen_sorted.txt"
DEFAULT_PHOTOMETRY = PACKAGE_DIR / "ngc4945_nuclear_photometry.csv"
DEFAULT_COMPONENT_AUDIT = PACKAGE_DIR / "ngc4945_component_audit.csv"
DEFAULT_CORRECTIONS = PACKAGE_DIR / "validated_continuum_corrections.csv"

CORRECTION_SCHEMA_VERSION = "2"
LEGACY_CORRECTION_FIELDS = (
    "source",
    "telescope",
    "frequency_ghz",
    "date_label",
    "product",
    "old_flux_field",
    "new_flux_field",
    "old_snr",
    "new_snr",
    "correction_reasons",
    "validated_unrounded_flux_mjy",
    "validated_background_only_snr",
)
CORRECTION_FIELDS = (
    "schema_version",
    *LEGACY_CORRECTION_FIELDS[:-1],
    "validated_adopted_noise_snr",
)

PHOTOMETRY_FIELDS = (
    "source",
    "telescope",
    "frequency_ghz",
    "date_label",
    "product",
    "input_sha256",
    "coordinate_frame",
    "fixed_ra",
    "fixed_dec",
    "coordinate_reference",
    "component_identity",
    "centroid_policy",
    "beam_policy",
    "fit_box_policy",
    "pixel_x_zero_based",
    "pixel_y_zero_based",
    "beam_major_arcsec",
    "beam_minor_arcsec",
    "beam_pa_deg",
    "fit_region_width_arcsec",
    "old_flux_field",
    "new_flux_field",
    "old_snr",
    "new_snr",
    "unrounded_flux_mjy",
    "adopted_noise_estimator",
    "adopted_noise_mjy_per_beam",
    "noise_rms_mjy_per_beam",
    "brightness_unit",
    "noise_aperture_policy",
    "noise_aperture_width_arcsec",
    "noise_aperture_pre_exclusion_npts",
    "noise_aperture_selected_pixel_area_over_requested_area",
    "noise_mask_major_semiaxis_arcsec",
    "noise_mask_minor_semiaxis_arcsec",
    "noise_mask_pa_deg",
    "fit_beam_major_delta_arcsec",
    "fit_beam_minor_delta_arcsec",
    "fit_beam_pa_delta_deg",
    "fit_beam_pa_guard_applicable",
    "image_width_pixels",
    "image_height_pixels",
    "noise_aperture_half_width_x_pixels",
    "noise_aperture_half_width_y_pixels",
    "noise_aperture_left_clearance_pixels",
    "noise_aperture_right_clearance_pixels",
    "noise_aperture_bottom_clearance_pixels",
    "noise_aperture_top_clearance_pixels",
    "noise_aperture_minimum_edge_clearance_pixels",
    "noise_aperture_edge_status",
    "noise_npts",
    "pixel_width_arcsec",
    "pixel_height_arcsec",
    "pixel_area_arcsec2",
    "valid_unmasked_area_arcsec2",
    "gaussian_beam_area_arcsec2",
    "independent_beam_areas",
    "effective_independent_samples",
    "noise_mean_mjy_per_beam",
    "noise_sigma_mjy_per_beam",
    "noise_scaled_mad_mjy_per_beam",
    "noise_review_flags",
    "noise_review_disposition",
    "noise_guard_status",
    "row_disposition",
    "status",
    "rejected_seed_identity",
    "rejected_seed_ra_icrs",
    "rejected_seed_dec_icrs",
)

# The pre-local-noise release is accepted only as an exact historical lineage
# input.  Keeping its schema explicit prevents a partial/stale sidecar from
# silently supplying the preserved pre-component flux and S/N fields.
LEGACY_PHOTOMETRY_FIELDS = (
    "source",
    "telescope",
    "frequency_ghz",
    "date_label",
    "product",
    "input_sha256",
    "coordinate_frame",
    "fixed_ra",
    "fixed_dec",
    "coordinate_reference",
    "centroid_policy",
    "beam_policy",
    "fit_box_policy",
    "pixel_x_zero_based",
    "pixel_y_zero_based",
    "beam_major_arcsec",
    "beam_minor_arcsec",
    "beam_pa_deg",
    "fit_region_width_arcsec",
    "old_flux_field",
    "new_flux_field",
    "old_snr",
    "new_snr",
    "unrounded_flux_mjy",
    "background_rms_mjy_per_beam",
    "status",
    "rejected_seed_identity",
    "rejected_seed_ra_icrs",
    "rejected_seed_dec_icrs",
)

COMPONENT_FIELDS = (
    "record_id",
    "decision",
    "component_identity",
    "frequency_ghz",
    "product",
    "coordinate_frame",
    "ra",
    "dec",
    "coordinate_reference",
    "centroid_policy",
    "beam_policy",
    "fit_region_width_arcsec",
    "flux_mjy",
    "image_scatter_mjy_per_beam",
    "uncertainty_basis",
    "included_in_canonical_table",
    "caveat",
)


@dataclass(frozen=True)
class Measurement:
    flux_mjy: float
    rms_mjy_per_beam: float
    adopted_noise_mjy_per_beam: float
    snr: float
    pixel_x: float
    pixel_y: float
    beam_major_arcsec: float
    beam_minor_arcsec: float
    beam_pa_deg: float
    fit_region_width_arcsec: float
    fit_ra_icrs: str
    fit_dec_icrs: str
    brightness_unit: str
    noise_aperture_width_arcsec: float
    noise_aperture_pre_exclusion_npts: int
    noise_aperture_selected_pixel_area_over_requested_area: float
    noise_mask_major_semiaxis_arcsec: float
    noise_mask_minor_semiaxis_arcsec: float
    noise_mask_pa_deg: float
    fit_beam_major_delta_arcsec: float
    fit_beam_minor_delta_arcsec: float
    fit_beam_pa_delta_deg: float
    fit_beam_pa_guard_applicable: bool
    image_width_pixels: int
    image_height_pixels: int
    noise_aperture_half_width_x_pixels: float
    noise_aperture_half_width_y_pixels: float
    noise_aperture_left_clearance_pixels: float
    noise_aperture_right_clearance_pixels: float
    noise_aperture_bottom_clearance_pixels: float
    noise_aperture_top_clearance_pixels: float
    noise_aperture_minimum_edge_clearance_pixels: float
    noise_npts: int
    pixel_width_arcsec: float
    pixel_height_arcsec: float
    pixel_area_arcsec2: float
    valid_unmasked_area_arcsec2: float
    gaussian_beam_area_arcsec2: float
    independent_beam_areas: float
    effective_independent_samples: float
    noise_mean_mjy_per_beam: float
    noise_sigma_mjy_per_beam: float
    noise_scaled_mad_mjy_per_beam: float
    noise_review_flags: tuple[str, ...]


def format_two_decimal(value: float) -> str:
    """Match the canonical continuum table's ordinary two-decimal format."""

    return f"{value:.2f}"


def portable_product(value: object) -> str:
    product = str(value).strip().replace("\\", "/")
    marker = "/multi_freq_from_archive/"
    if marker in product:
        product = product.split(marker, 1)[1]
    return product


def read_sorted(path: Path) -> list[list[Any]]:
    rows = [
        ast.literal_eval(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    return rows


def ngc4945_rows(rows: Iterable[list[Any]]) -> list[list[Any]]:
    selected = [row for row in rows if str(row[0]) == SOURCE]
    if len(selected) != EXPECTED_ROWS:
        raise RuntimeError(
            f"expected {EXPECTED_ROWS} {SOURCE} rows, found {len(selected)}"
        )
    return selected


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _component_direction_icrs(component: dict[str, Any]) -> tuple[str, str]:
    from astropy import units as u
    from astropy.coordinates import FK5, SkyCoord
    from astropy.time import Time

    direction = component["shape"]["direction"]
    longitude = float(direction["m0"]["value"]) * u.rad
    latitude = float(direction["m1"]["value"]) * u.rad
    refer = str(direction.get("refer", "ICRS")).upper()
    if refer in {"J2000", "FK5"}:
        sky = SkyCoord(
            longitude,
            latitude,
            frame=FK5(equinox=Time("J2000")),
        ).icrs
    else:
        sky = SkyCoord(longitude, latitude, frame="icrs")
    return (
        sky.ra.to_string(unit=u.hourangle, sep=":", precision=6, pad=True),
        sky.dec.to_string(unit=u.deg, sep=":", precision=6, alwayssign=True, pad=True),
    )


def _fits_coordinate_and_beam(
    path: Path,
) -> tuple[Any, float, float, float, tuple[int, int], float, float, float]:
    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from astropy.wcs import FITSFixedWarning, WCS
    from astropy.wcs.utils import proj_plane_pixel_scales

    warnings.simplefilter("ignore", FITSFixedWarning)
    with fits.open(path, memmap=True) as hdul:
        header = hdul[0].header
        wcs = WCS(header).celestial
        pixel_scales_arcsec = tuple(
            float(value) * 3600.0 for value in proj_plane_pixel_scales(wcs)
        )
        if len(pixel_scales_arcsec) != 2 or any(
            not math.isfinite(value) or value <= 0.0
            for value in pixel_scales_arcsec
        ):
            raise RuntimeError(f"invalid celestial pixel scale for {path}")
        shape = (int(header["NAXIS2"]), int(header["NAXIS1"]))
        major = float(header["BMAJ"]) * 3600.0
        minor = float(header["BMIN"]) * 3600.0
        pa = float(header["BPA"])
        frequency_hz = math.nan
        for axis in range(1, int(header.get("NAXIS", 0)) + 1):
            if str(header.get(f"CTYPE{axis}", "")).upper().startswith("FREQ"):
                frequency_hz = float(header[f"CRVAL{axis}"])
                break

    coordinate = SkyCoord(
        f"{FIXED_RA_ICRS} {FIXED_DEC_ICRS}",
        unit=("hourangle", "deg"),
        frame="icrs",
    )
    x, y = (float(value) for value in wcs.world_to_pixel(coordinate))
    return (
        (x, y), major, minor, pa, shape, frequency_hz,
        pixel_scales_arcsec[0], pixel_scales_arcsec[1],
    )


def measure(
    path: Path,
    *,
    width_arcsec: float | None = None,
    fixed_centroid: bool = True,
    calculate_rms: bool = True,
) -> Measurement:
    """Fit one beam-shaped component at the source-18 nuclear position."""

    from casatasks import imfit, imhead, imstat

    (
        (x, y), major, minor, pa, shape, _, pixel_width, pixel_height,
    ) = _fits_coordinate_and_beam(path)
    image_height, image_width = shape
    width = major * FIT_BOX_FACTOR if width_arcsec is None else width_arcsec
    region = f"centerbox[[{x}pix,{y}pix],[{width}arcsec,{width}arcsec]]"
    fix = "xyabp" if fixed_centroid else "abp"

    with tempfile.TemporaryDirectory(prefix="ngc4945-nuclear-") as tmp:
        tmpdir = Path(tmp)
        estimates = tmpdir / "estimates.txt"
        estimates.write_text(
            f"0.001,{x},{y},{major}arcsec,{minor}arcsec,{pa}deg,{fix}\n",
            encoding="utf-8",
        )
        result = imfit(
            imagename=str(path),
            region=region,
            estimates=str(estimates),
        )
        if not bool(result["converged"][0]):
            raise RuntimeError(f"CASA imfit did not converge for {path}")

        deconvolved = result["deconvolved"]["component0"]
        flux_value = float(deconvolved["flux"]["value"][0])
        flux_unit = str(deconvolved["flux"]["unit"]).lower()
        if flux_unit == "jy":
            flux_mjy = flux_value * 1000.0
        elif flux_unit == "mjy":
            flux_mjy = flux_value
        else:
            raise RuntimeError(f"unsupported CASA component flux unit {flux_unit!r}")

        fitted = result["results"]["component0"]
        fit_x, fit_y = (float(value) for value in fitted["pixelcoords"])
        fit_ra, fit_dec = _component_direction_icrs(fitted)

        rms_mjy = math.nan
        noise_aperture_width = math.nan
        noise_aperture_pre_exclusion_npts = 0
        noise_aperture_selected_pixel_area_over_requested_area = math.nan
        noise_mask_major = math.nan
        noise_mask_minor = math.nan
        noise_mask_pa = math.nan
        noise_npts = 0
        pixel_area = pixel_width * pixel_height
        valid_unmasked_area = math.nan
        beam_area = gaussian_beam_area_arcsec2(major, minor)
        beam_areas = math.nan
        effective_samples = math.nan
        noise_mean = math.nan
        noise_sigma = math.nan
        noise_scaled_mad = math.nan
        noise_review_flags: tuple[str, ...] = ()
        brightness_unit_text = ""
        fit_major_delta = math.nan
        fit_minor_delta = math.nan
        fit_pa_delta = math.nan
        fit_pa_guard_applicable = False
        aperture_half_width_x = math.nan
        aperture_half_width_y = math.nan
        aperture_left_clearance = math.nan
        aperture_right_clearance = math.nan
        aperture_bottom_clearance = math.nan
        aperture_top_clearance = math.nan
        aperture_minimum_clearance = math.nan
        if calculate_rms:
            fit_major = float(fitted["shape"]["majoraxis"]["value"])
            fit_minor = float(fitted["shape"]["minoraxis"]["value"])
            fit_major_unit = str(fitted["shape"]["majoraxis"]["unit"]).lower()
            fit_minor_unit = str(fitted["shape"]["minoraxis"]["unit"]).lower()
            if fit_major_unit == "deg":
                fit_major *= 3600.0
            elif fit_major_unit != "arcsec":
                raise RuntimeError(
                    f"unsupported fitted major-axis unit {fit_major_unit!r}"
                )
            if fit_minor_unit == "deg":
                fit_minor *= 3600.0
            elif fit_minor_unit != "arcsec":
                raise RuntimeError(
                    f"unsupported fitted minor-axis unit {fit_minor_unit!r}"
                )
            fit_pa = float(fitted["shape"]["positionangle"]["value"])
            if str(fitted["shape"]["positionangle"]["unit"]).lower() != "deg":
                raise RuntimeError("unsupported fitted position-angle unit")
            fixed_shape = require_fixed_restoring_beam_shape(
                fit_major, fit_minor, fit_pa, major, minor, pa
            )
            fit_major_delta = fixed_shape.major_delta_arcsec
            fit_minor_delta = fixed_shape.minor_delta_arcsec
            fit_pa_delta = fixed_shape.pa_delta_deg
            fit_pa_guard_applicable = fixed_shape.pa_guard_applicable
            noise_aperture_width = major * FINAL_NOISE_APERTURE_BMAJ_MULTIPLIER
            coverage = square_aperture_coverage(
                fit_x,
                fit_y,
                noise_aperture_width,
                pixel_width,
                pixel_height,
                image_width,
                image_height,
            )
            aperture_half_width_x = coverage.half_width_x_pixels
            aperture_half_width_y = coverage.half_width_y_pixels
            aperture_left_clearance = coverage.left_clearance_pixels
            aperture_right_clearance = coverage.right_clearance_pixels
            aperture_bottom_clearance = coverage.bottom_clearance_pixels
            aperture_top_clearance = coverage.top_clearance_pixels
            aperture_minimum_clearance = coverage.minimum_edge_clearance_pixels
            noise_mask_major = fit_major * THREE_SIGMA_FROM_FWHM
            noise_mask_minor = fit_minor * THREE_SIGMA_FROM_FWHM
            noise_mask_pa = fit_pa
            noise_aperture_region = (
                f"centerbox[[{fit_x}pix,{fit_y}pix],"
                f"[{noise_aperture_width}arcsec,{noise_aperture_width}arcsec]]"
            )
            noise_region = tmpdir / "source_excluded_background.crtf"
            noise_region.write_text(
                "#CRTFv0\n"
                # CASAdocs: centerbox is full width; ellipse axes are semiaxes.
                f"{noise_aperture_region}\n"
                f"-ellipse[[{fit_x}pix,{fit_y}pix],"
                f"[{noise_mask_major}arcsec,"
                f"{noise_mask_minor}arcsec],{fit_pa}deg]\n",
                encoding="utf-8",
            )
            aperture_statistics = imstat(
                imagename=str(path), region=noise_aperture_region
            )
            statistics = imstat(imagename=str(path), region=str(noise_region))
            brightness_unit = imhead(
                imagename=str(path), mode="get", hdkey="bunit"
            )
            if isinstance(brightness_unit, dict):
                brightness_unit = brightness_unit.get(
                    "value", brightness_unit.get("unit", "")
                )
            brightness_unit_text = str(brightness_unit)
            conversion = brightness_to_mjy_per_beam_factor(brightness_unit)

            def stat_value(result: dict[str, Any], field: str) -> float:
                if field not in result:
                    raise RuntimeError(f"CASA imstat omitted {field} for {path}")
                try:
                    value = result[field][0]
                except (IndexError, KeyError, TypeError):
                    value = result[field]
                value = float(value)
                if not math.isfinite(value):
                    raise RuntimeError(f"CASA imstat {field} is non-finite for {path}")
                return value

            aperture_npts_value = stat_value(aperture_statistics, "npts")
            if aperture_npts_value < 0.0 or not aperture_npts_value.is_integer():
                raise RuntimeError(
                    f"invalid CASA aperture npts for {path}: {aperture_npts_value}"
                )
            noise_aperture_pre_exclusion_npts = int(aperture_npts_value)
            # CASA rasterizes the continuous square to selected pixels, so
            # this diagnostic area ratio can be slightly greater than one.
            noise_aperture_selected_pixel_area_over_requested_area = (
                noise_aperture_pre_exclusion_npts
                * pixel_area
                / (noise_aperture_width * noise_aperture_width)
            )
            noise_npts_value = stat_value(statistics, "npts")
            if noise_npts_value < 0.0 or not noise_npts_value.is_integer():
                raise RuntimeError(f"invalid CASA npts for {path}: {noise_npts_value}")
            noise_npts = int(noise_npts_value)
            rms_mjy = stat_value(statistics, "rms") * conversion
            noise_sigma = stat_value(statistics, "sigma") * conversion
            noise_mean = stat_value(statistics, "mean") * conversion
            noise_scaled_mad = (
                stat_value(statistics, "medabsdevmed")
                * conversion
                * SCALED_MAD_NORMAL_FACTOR
            )
            valid_unmasked_area = noise_npts * pixel_area
            beam_areas = independent_beam_areas(
                noise_npts, pixel_width, pixel_height, major, minor
            )
            effective_samples = effective_independent_noise_samples(
                noise_npts, pixel_width, pixel_height, major, minor
            )
            require_minimum_independent_beams(effective_samples)
            noise_review_flags = final_noise_review_flags(
                noise_mean,
                noise_sigma,
                rms_mjy,
                noise_scaled_mad,
                aperture_minimum_edge_clearance_pixels=(
                    coverage.minimum_edge_clearance_pixels
                ),
            )

    adopted_noise = noise_scaled_mad if calculate_rms else math.nan
    snr = flux_mjy / adopted_noise if calculate_rms else math.nan
    return Measurement(
        flux_mjy=flux_mjy,
        rms_mjy_per_beam=rms_mjy,
        adopted_noise_mjy_per_beam=adopted_noise,
        snr=snr,
        pixel_x=fit_x,
        pixel_y=fit_y,
        beam_major_arcsec=major,
        beam_minor_arcsec=minor,
        beam_pa_deg=pa,
        fit_region_width_arcsec=width,
        fit_ra_icrs=fit_ra,
        fit_dec_icrs=fit_dec,
        brightness_unit=brightness_unit_text,
        noise_aperture_width_arcsec=noise_aperture_width,
        noise_aperture_pre_exclusion_npts=noise_aperture_pre_exclusion_npts,
        noise_aperture_selected_pixel_area_over_requested_area=(
            noise_aperture_selected_pixel_area_over_requested_area
        ),
        noise_mask_major_semiaxis_arcsec=noise_mask_major,
        noise_mask_minor_semiaxis_arcsec=noise_mask_minor,
        noise_mask_pa_deg=noise_mask_pa,
        fit_beam_major_delta_arcsec=fit_major_delta,
        fit_beam_minor_delta_arcsec=fit_minor_delta,
        fit_beam_pa_delta_deg=fit_pa_delta,
        fit_beam_pa_guard_applicable=fit_pa_guard_applicable,
        image_width_pixels=image_width,
        image_height_pixels=image_height,
        noise_aperture_half_width_x_pixels=aperture_half_width_x,
        noise_aperture_half_width_y_pixels=aperture_half_width_y,
        noise_aperture_left_clearance_pixels=aperture_left_clearance,
        noise_aperture_right_clearance_pixels=aperture_right_clearance,
        noise_aperture_bottom_clearance_pixels=aperture_bottom_clearance,
        noise_aperture_top_clearance_pixels=aperture_top_clearance,
        noise_aperture_minimum_edge_clearance_pixels=aperture_minimum_clearance,
        noise_npts=noise_npts,
        pixel_width_arcsec=pixel_width,
        pixel_height_arcsec=pixel_height,
        pixel_area_arcsec2=pixel_area,
        valid_unmasked_area_arcsec2=valid_unmasked_area,
        gaussian_beam_area_arcsec2=beam_area,
        independent_beam_areas=beam_areas,
        effective_independent_samples=effective_samples,
        noise_mean_mjy_per_beam=noise_mean,
        noise_sigma_mjy_per_beam=noise_sigma,
        noise_scaled_mad_mjy_per_beam=noise_scaled_mad,
        noise_review_flags=noise_review_flags,
    )


def build_records(
    source_table: Path,
    archive_root: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows = ngc4945_rows(read_sorted(source_table))
    photometry: list[dict[str, object]] = []
    anchor_path: Path | None = None
    anchor_measurement: Measurement | None = None
    anchor_row: list[Any] | None = None

    for row in rows:
        product = portable_product(row[7])
        path = archive_root / product
        if not path.is_file():
            raise FileNotFoundError(f"missing archive product: {path}")
        measurement = measure(path)
        frequency = float(row[2])
        beam_geom = math.sqrt(
            measurement.beam_major_arcsec * measurement.beam_minor_arcsec
        )
        if not math.isclose(beam_geom, float(row[5]), rel_tol=0.0, abs_tol=5e-5):
            raise RuntimeError(
                f"beam mismatch for {product}: header {beam_geom}, table {row[5]}"
            )
        classification = classify_ordinary_fit(
            measurement.flux_mjy,
            measurement.adopted_noise_mjy_per_beam,
        )
        status = "detection" if classification.is_detection else "upper_limit"
        if status == "detection":
            flux_field = format_two_decimal(measurement.flux_mjy)
            snr_field = f"{classification.snr:.3f}"
        else:
            flux_field = format_two_decimal(classification.upper_limit_mjy) + "*"
            snr_field = "NA"

        photometry.append(
            {
                "source": SOURCE,
                "telescope": row[1],
                "frequency_ghz": row[2],
                "date_label": row[4],
                "product": product,
                "input_sha256": sha256(path),
                "coordinate_frame": COORDINATE_FRAME,
                "fixed_ra": FIXED_RA_ICRS,
                "fixed_dec": FIXED_DEC_ICRS,
                "coordinate_reference": COORDINATE_REFERENCE,
                "component_identity": "Emig2020_source18_AGN_maser_nucleus",
                "centroid_policy": "fixed_xy",
                "beam_policy": "fixed_restoring_beam_abp",
                "fit_box_policy": (
                    "legacy_square_full_width_3_gaussian_sigma_equals_1p274_bmaj"
                ),
                "pixel_x_zero_based": measurement.pixel_x,
                "pixel_y_zero_based": measurement.pixel_y,
                "beam_major_arcsec": measurement.beam_major_arcsec,
                "beam_minor_arcsec": measurement.beam_minor_arcsec,
                "beam_pa_deg": measurement.beam_pa_deg,
                "fit_region_width_arcsec": measurement.fit_region_width_arcsec,
                "old_flux_field": row[3],
                "new_flux_field": flux_field,
                "old_snr": row[6],
                "new_snr": snr_field,
                "unrounded_flux_mjy": measurement.flux_mjy,
                "adopted_noise_estimator": "scaled_mad",
                "adopted_noise_mjy_per_beam": (
                    measurement.adopted_noise_mjy_per_beam
                ),
                "noise_rms_mjy_per_beam": measurement.rms_mjy_per_beam,
                "brightness_unit": measurement.brightness_unit,
                "noise_aperture_policy": (
                    "fixed_component_centered_10bmaj_full_width_"
                    "minus_fitted_3sigma_semiaxis_ellipse"
                ),
                "noise_aperture_width_arcsec": measurement.noise_aperture_width_arcsec,
                "noise_aperture_pre_exclusion_npts": (
                    measurement.noise_aperture_pre_exclusion_npts
                ),
                "noise_aperture_selected_pixel_area_over_requested_area": (
                    measurement.noise_aperture_selected_pixel_area_over_requested_area
                ),
                "noise_mask_major_semiaxis_arcsec": (
                    measurement.noise_mask_major_semiaxis_arcsec
                ),
                "noise_mask_minor_semiaxis_arcsec": (
                    measurement.noise_mask_minor_semiaxis_arcsec
                ),
                "noise_mask_pa_deg": measurement.noise_mask_pa_deg,
                "fit_beam_major_delta_arcsec": (
                    measurement.fit_beam_major_delta_arcsec
                ),
                "fit_beam_minor_delta_arcsec": (
                    measurement.fit_beam_minor_delta_arcsec
                ),
                "fit_beam_pa_delta_deg": measurement.fit_beam_pa_delta_deg,
                "fit_beam_pa_guard_applicable": (
                    "yes" if measurement.fit_beam_pa_guard_applicable else "no"
                ),
                "image_width_pixels": measurement.image_width_pixels,
                "image_height_pixels": measurement.image_height_pixels,
                "noise_aperture_half_width_x_pixels": (
                    measurement.noise_aperture_half_width_x_pixels
                ),
                "noise_aperture_half_width_y_pixels": (
                    measurement.noise_aperture_half_width_y_pixels
                ),
                "noise_aperture_left_clearance_pixels": (
                    measurement.noise_aperture_left_clearance_pixels
                ),
                "noise_aperture_right_clearance_pixels": (
                    measurement.noise_aperture_right_clearance_pixels
                ),
                "noise_aperture_bottom_clearance_pixels": (
                    measurement.noise_aperture_bottom_clearance_pixels
                ),
                "noise_aperture_top_clearance_pixels": (
                    measurement.noise_aperture_top_clearance_pixels
                ),
                "noise_aperture_minimum_edge_clearance_pixels": (
                    measurement.noise_aperture_minimum_edge_clearance_pixels
                ),
                "noise_aperture_edge_status": (
                    "truncated_available_valid_intersection"
                    if measurement.noise_aperture_minimum_edge_clearance_pixels < -1e-9
                    else "fully_contained"
                ),
                "noise_npts": measurement.noise_npts,
                "pixel_width_arcsec": measurement.pixel_width_arcsec,
                "pixel_height_arcsec": measurement.pixel_height_arcsec,
                "pixel_area_arcsec2": measurement.pixel_area_arcsec2,
                "valid_unmasked_area_arcsec2": measurement.valid_unmasked_area_arcsec2,
                "gaussian_beam_area_arcsec2": measurement.gaussian_beam_area_arcsec2,
                "independent_beam_areas": measurement.independent_beam_areas,
                "effective_independent_samples": measurement.effective_independent_samples,
                "noise_mean_mjy_per_beam": measurement.noise_mean_mjy_per_beam,
                "noise_sigma_mjy_per_beam": measurement.noise_sigma_mjy_per_beam,
                "noise_scaled_mad_mjy_per_beam": (
                    measurement.noise_scaled_mad_mjy_per_beam
                ),
                "noise_review_flags": "+".join(measurement.noise_review_flags) or "none",
                "noise_review_disposition": (
                    "diagnostic_flags_recorded_not_classification_overrides"
                    if measurement.noise_review_flags
                    else "no_diagnostic_flags"
                ),
                "noise_guard_status": "passed",
                "row_disposition": "retained_source18_replacement",
                "status": status,
                "rejected_seed_identity": "Emig2020_source17_SSC_SNR_associated",
                "rejected_seed_ra_icrs": "13:05:27.4933",
                "rejected_seed_dec_icrs": "-49:28:05.1422",
            }
        )

        if math.isclose(frequency, 349.40, rel_tol=0.0, abs_tol=1e-6):
            anchor_path = path
            anchor_measurement = measurement
            anchor_row = row

    if anchor_path is None or anchor_measurement is None or anchor_row is None:
        raise RuntimeError("the 349.40 GHz NGC 4945 anchor row was not found")

    component_audit: list[dict[str, object]] = [
        {
            "record_id": "rejected-source17-archived-anchor",
            "decision": "rejected",
            "component_identity": "Emig2020_source17_SSC_SNR_associated",
            "frequency_ghz": anchor_row[2],
            "product": portable_product(anchor_row[7]),
            "coordinate_frame": "ICRS",
            "ra": "13:05:27.4933",
            "dec": "-49:28:05.1422",
            "coordinate_reference": "archived_CASA_analytic_centroid_and_Emig2020_source17",
            "centroid_policy": "historical_free_centroid",
            "beam_policy": "fixed_restoring_beam_abp",
            "fit_region_width_arcsec": 0.11434164558049752,
            "flux_mjy": 15.841014462851778,
            "image_scatter_mjy_per_beam": 15.84 / 61.054,
            "uncertainty_basis": "archived_table_flux_over_image_scatter_snr",
            "included_in_canonical_table": "no",
            "caveat": "wrong physical component; retained for audit only",
        }
    ]

    checks = (
        ("accepted-source18-fixed-canonical", True, None, "yes"),
        ("source18-free-centroid-check", False, None, "no"),
        ("source18-fixed-width-0p18", True, 0.18, "no"),
        ("source18-free-width-0p18", False, 0.18, "no"),
        ("source18-fixed-width-0p24", True, 0.24, "no"),
        ("source18-free-width-0p24", False, 0.24, "no"),
    )
    for record_id, fixed, width, canonical in checks:
        result = (
            anchor_measurement
            if fixed and width is None
            else measure(
                anchor_path,
                width_arcsec=width,
                fixed_centroid=fixed,
                calculate_rms=False,
            )
        )
        is_canonical = canonical == "yes"
        component_audit.append(
            {
                "record_id": record_id,
                "decision": "accepted" if is_canonical else "diagnostic",
                "component_identity": "Emig2020_source18_AGN_maser_nucleus",
                "frequency_ghz": anchor_row[2],
                "product": portable_product(anchor_row[7]),
                "coordinate_frame": "ICRS",
                "ra": result.fit_ra_icrs,
                "dec": result.fit_dec_icrs,
                "coordinate_reference": COORDINATE_REFERENCE,
                "centroid_policy": "fixed_xy" if fixed else "free_xy_validation_only",
                "beam_policy": "fixed_restoring_beam_abp",
                "fit_region_width_arcsec": result.fit_region_width_arcsec,
                "flux_mjy": result.flux_mjy,
                "image_scatter_mjy_per_beam": (
                    anchor_measurement.adopted_noise_mjy_per_beam
                ),
                "uncertainty_basis": (
                    "10bmaj_source_excluded_scaled_mad_once; "
                    "calibration_and_model_terms_excluded"
                ),
                "included_in_canonical_table": canonical,
                "caveat": (
                    "canonical released fit-box policy"
                    if is_canonical
                    else "region/centroid sensitivity diagnostic; not an additional independent error"
                ),
            }
        )
    return photometry, component_audit


def write_csv(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_exact_csv(
    path: Path,
    fields: tuple[str, ...],
    *,
    label: str,
) -> list[dict[str, str]]:
    """Read one audit CSV without accepting missing, extra, or reordered fields."""

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != fields:
            raise RuntimeError(
                f"{label} has fields {reader.fieldnames!r}; expected {fields!r}"
            )
        rows = list(reader)
    if any(any(row.get(field) is None for field in fields) for row in rows):
        raise RuntimeError(f"{label} contains an incomplete record")
    return rows


def _finite_float(
    record: dict[str, object],
    field: str,
    *,
    positive: bool = False,
) -> float:
    try:
        value = float(record[field])
    except (KeyError, TypeError, ValueError) as error:
        raise RuntimeError(f"invalid {field} in {record!r}") from error
    if not math.isfinite(value) or (positive and value <= 0.0):
        qualifier = "finite and positive" if positive else "finite"
        raise RuntimeError(f"{field} must be {qualifier} in {record!r}")
    return value


def input_hash_manifest_sha256(photometry: list[dict[str, object]]) -> str:
    """Hash the complete product-to-input-hash mapping in a stable order."""

    payload = "\n".join(
        f"{row['product']}\0{row['input_sha256']}"
        for row in sorted(photometry, key=lambda item: str(item["product"]))
    ) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def validate_photometry_records(photometry: list[dict[str, object]]) -> None:
    """Validate the fixed-nucleus sidecar independently of CASA and FITS."""

    if len(photometry) != EXPECTED_ROWS:
        raise RuntimeError(
            f"expected {EXPECTED_ROWS} photometry records, found {len(photometry)}"
        )
    if any(set(record) != set(PHOTOMETRY_FIELDS) for record in photometry):
        raise RuntimeError("photometry record fields do not match the released schema")

    products = [str(record["product"]) for record in photometry]
    if len(set(products)) != EXPECTED_ROWS:
        raise RuntimeError(
            f"photometry products are not exactly {EXPECTED_ROWS} unique records"
        )
    if any(portable_product(product) != product or product.startswith("/") for product in products):
        raise RuntimeError("photometry products must use release-portable paths")

    for record in photometry:
        product = str(record["product"])
        recorded_hash = str(record["input_sha256"])
        if len(recorded_hash) != 64 or any(
            character not in "0123456789abcdef" for character in recorded_hash
        ):
            raise RuntimeError(f"invalid input SHA-256 for {product}")

        expected_identity = {
            "source": SOURCE,
            "coordinate_frame": COORDINATE_FRAME,
            "fixed_ra": FIXED_RA_ICRS,
            "fixed_dec": FIXED_DEC_ICRS,
            "coordinate_reference": COORDINATE_REFERENCE,
            "component_identity": "Emig2020_source18_AGN_maser_nucleus",
            "centroid_policy": "fixed_xy",
            "beam_policy": "fixed_restoring_beam_abp",
            "fit_box_policy": (
                "legacy_square_full_width_3_gaussian_sigma_equals_1p274_bmaj"
            ),
            "noise_aperture_policy": (
                "fixed_component_centered_10bmaj_full_width_"
                "minus_fitted_3sigma_semiaxis_ellipse"
            ),
            "adopted_noise_estimator": "scaled_mad",
            "noise_guard_status": "passed",
            "row_disposition": "retained_source18_replacement",
            "rejected_seed_identity": "Emig2020_source17_SSC_SNR_associated",
            "rejected_seed_ra_icrs": "13:05:27.4933",
            "rejected_seed_dec_icrs": "-49:28:05.1422",
        }
        for field, expected in expected_identity.items():
            if str(record[field]) != expected:
                raise RuntimeError(f"unexpected {field} for {product}")
        if str(record["component_identity"]) == str(record["rejected_seed_identity"]):
            raise RuntimeError(
                f"retained source-18 and rejected source-17 identities collide for {product}"
            )

        flux = _finite_float(record, "unrounded_flux_mjy", positive=True)
        adopted_noise = _finite_float(
            record,
            "adopted_noise_mjy_per_beam",
            positive=True,
        )
        diagnostic_rms = _finite_float(
            record, "noise_rms_mjy_per_beam", positive=True
        )
        beam_major = _finite_float(record, "beam_major_arcsec", positive=True)
        beam_minor = _finite_float(record, "beam_minor_arcsec", positive=True)
        width = _finite_float(record, "fit_region_width_arcsec", positive=True)
        noise_width = _finite_float(record, "noise_aperture_width_arcsec", positive=True)
        aperture_pre_npts = _finite_float(
            record, "noise_aperture_pre_exclusion_npts", positive=True
        )
        if not aperture_pre_npts.is_integer():
            raise RuntimeError(
                f"noise_aperture_pre_exclusion_npts must be an integer for {product}"
            )
        aperture_selected_area_ratio = _finite_float(
            record,
            "noise_aperture_selected_pixel_area_over_requested_area",
            positive=True,
        )
        mask_major = _finite_float(
            record, "noise_mask_major_semiaxis_arcsec", positive=True
        )
        mask_minor = _finite_float(
            record, "noise_mask_minor_semiaxis_arcsec", positive=True
        )
        mask_pa = _finite_float(record, "noise_mask_pa_deg")
        recorded_major_delta = _finite_float(
            record, "fit_beam_major_delta_arcsec"
        )
        recorded_minor_delta = _finite_float(
            record, "fit_beam_minor_delta_arcsec"
        )
        recorded_pa_delta = _finite_float(record, "fit_beam_pa_delta_deg")
        if str(record["fit_beam_pa_guard_applicable"]) not in {"yes", "no"}:
            raise RuntimeError(f"invalid fit-beam PA guard flag for {product}")
        image_width = _finite_float(record, "image_width_pixels", positive=True)
        image_height = _finite_float(record, "image_height_pixels", positive=True)
        if not image_width.is_integer() or not image_height.is_integer():
            raise RuntimeError(f"image dimensions must be integers for {product}")
        recorded_half_x = _finite_float(
            record, "noise_aperture_half_width_x_pixels", positive=True
        )
        recorded_half_y = _finite_float(
            record, "noise_aperture_half_width_y_pixels", positive=True
        )
        recorded_left = _finite_float(
            record, "noise_aperture_left_clearance_pixels"
        )
        recorded_right = _finite_float(
            record, "noise_aperture_right_clearance_pixels"
        )
        recorded_bottom = _finite_float(
            record, "noise_aperture_bottom_clearance_pixels"
        )
        recorded_top = _finite_float(
            record, "noise_aperture_top_clearance_pixels"
        )
        recorded_minimum_clearance = _finite_float(
            record, "noise_aperture_minimum_edge_clearance_pixels"
        )
        noise_npts = _finite_float(record, "noise_npts", positive=True)
        if not noise_npts.is_integer():
            raise RuntimeError(f"noise_npts must be an integer for {product}")
        if noise_npts > aperture_pre_npts:
            raise RuntimeError(
                f"source-excluded noise_npts exceeds aperture npts for {product}"
            )
        pixel_width = _finite_float(record, "pixel_width_arcsec", positive=True)
        pixel_height = _finite_float(record, "pixel_height_arcsec", positive=True)
        pixel_area = _finite_float(record, "pixel_area_arcsec2", positive=True)
        valid_area = _finite_float(
            record, "valid_unmasked_area_arcsec2", positive=True
        )
        recorded_beam_area = _finite_float(
            record, "gaussian_beam_area_arcsec2", positive=True
        )
        recorded_beam_areas = _finite_float(
            record, "independent_beam_areas", positive=True
        )
        effective_samples = _finite_float(
            record, "effective_independent_samples", positive=True
        )
        noise_mean = _finite_float(record, "noise_mean_mjy_per_beam")
        noise_sigma = _finite_float(
            record, "noise_sigma_mjy_per_beam", positive=True
        )
        noise_scaled_mad = _finite_float(record, "noise_scaled_mad_mjy_per_beam")
        if noise_scaled_mad <= 0.0:
            raise RuntimeError(f"noise scaled MAD must be positive for {product}")
        if not math.isclose(
            adopted_noise,
            noise_scaled_mad,
            rel_tol=0.0,
            abs_tol=1e-15,
        ):
            raise RuntimeError(f"adopted scaled MAD is stale for {product}")
        _finite_float(record, "frequency_ghz", positive=True)
        center_x = _finite_float(record, "pixel_x_zero_based")
        center_y = _finite_float(record, "pixel_y_zero_based")
        beam_pa = _finite_float(record, "beam_pa_deg")
        brightness_to_mjy_per_beam_factor(str(record["brightness_unit"]))
        if not math.isclose(
            width,
            beam_major * FIT_BOX_FACTOR,
            rel_tol=1e-12,
            abs_tol=1e-12,
        ):
            raise RuntimeError(f"fit-box width is inconsistent with the beam for {product}")
        if beam_minor > beam_major:
            raise RuntimeError(f"beam minor axis exceeds the major axis for {product}")
        fixed_shape = require_fixed_restoring_beam_shape(
            mask_major / THREE_SIGMA_FROM_FWHM,
            mask_minor / THREE_SIGMA_FROM_FWHM,
            mask_pa,
            beam_major,
            beam_minor,
            beam_pa,
        )
        coverage = square_aperture_coverage(
            center_x,
            center_y,
            noise_width,
            pixel_width,
            pixel_height,
            image_width,
            image_height,
        )
        expected_noise_values = {
            "noise aperture width": (
                noise_width,
                beam_major * FINAL_NOISE_APERTURE_BMAJ_MULTIPLIER,
            ),
            "valid unmasked area": (valid_area, noise_npts * pixel_area),
            "Gaussian beam area": (
                recorded_beam_area,
                gaussian_beam_area_arcsec2(beam_major, beam_minor),
            ),
            "independent beam areas": (
                recorded_beam_areas,
                valid_area / recorded_beam_area,
            ),
            "effective independent samples": (
                effective_samples,
                min(noise_npts, recorded_beam_areas),
            ),
            "pixel area": (pixel_area, pixel_width * pixel_height),
            "selected-pixel/requested-aperture area ratio": (
                aperture_selected_area_ratio,
                aperture_pre_npts * pixel_area / (noise_width * noise_width),
            ),
            "fit-beam major delta": (
                recorded_major_delta,
                fixed_shape.major_delta_arcsec,
            ),
            "fit-beam minor delta": (
                recorded_minor_delta,
                fixed_shape.minor_delta_arcsec,
            ),
            "fit-beam PA delta": (recorded_pa_delta, fixed_shape.pa_delta_deg),
            "aperture half width x": (
                recorded_half_x,
                coverage.half_width_x_pixels,
            ),
            "aperture half width y": (
                recorded_half_y,
                coverage.half_width_y_pixels,
            ),
            "aperture left clearance": (
                recorded_left,
                coverage.left_clearance_pixels,
            ),
            "aperture right clearance": (
                recorded_right,
                coverage.right_clearance_pixels,
            ),
            "aperture bottom clearance": (
                recorded_bottom,
                coverage.bottom_clearance_pixels,
            ),
            "aperture top clearance": (
                recorded_top,
                coverage.top_clearance_pixels,
            ),
            "aperture minimum clearance": (
                recorded_minimum_clearance,
                coverage.minimum_edge_clearance_pixels,
            ),
        }
        for label, (actual, expected) in expected_noise_values.items():
            if not math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12):
                raise RuntimeError(f"stale {label} for {product}: {actual} vs {expected}")
        expected_pa_guard = "yes" if fixed_shape.pa_guard_applicable else "no"
        if str(record["fit_beam_pa_guard_applicable"]) != expected_pa_guard:
            raise RuntimeError(f"stale fit-beam PA guard flag for {product}")
        expected_edge_status = (
            "truncated_available_valid_intersection"
            if coverage.minimum_edge_clearance_pixels < -1e-9
            else "fully_contained"
        )
        if str(record["noise_aperture_edge_status"]) != expected_edge_status:
            raise RuntimeError(f"stale noise-aperture edge status for {product}")
        require_minimum_independent_beams(effective_samples)
        expected_review_flags = "+".join(
            final_noise_review_flags(
                noise_mean,
                noise_sigma,
                diagnostic_rms,
                noise_scaled_mad,
                aperture_minimum_edge_clearance_pixels=(
                    coverage.minimum_edge_clearance_pixels
                ),
            )
        ) or "none"
        if str(record["noise_review_flags"]) != expected_review_flags:
            raise RuntimeError(f"stale noise review flags for {product}")
        expected_review_disposition = (
            "diagnostic_flags_recorded_not_classification_overrides"
            if expected_review_flags != "none"
            else "no_diagnostic_flags"
        )
        if str(record["noise_review_disposition"]) != expected_review_disposition:
            raise RuntimeError(f"stale noise review disposition for {product}")

        classification = classify_ordinary_fit(flux, adopted_noise)
        expected_status = "detection" if classification.is_detection else "upper_limit"
        expected_flux = (
            format_two_decimal(flux)
            if expected_status == "detection"
            else format_two_decimal(classification.upper_limit_mjy) + "*"
        )
        expected_snr = (
            f"{classification.snr:.3f}" if expected_status == "detection" else "NA"
        )
        if str(record["status"]) != expected_status:
            raise RuntimeError(
                f"status is inconsistent with flux/adopted noise for {product}"
            )
        if str(record["new_flux_field"]) != expected_flux:
            raise RuntimeError(f"rounded flux is inconsistent for {product}")
        if str(record["new_snr"]) != expected_snr:
            raise RuntimeError(
                f"S/N is inconsistent with flux/adopted noise for {product}"
            )

    actual_manifest_hash = input_hash_manifest_sha256(photometry)
    if actual_manifest_hash != EXPECTED_INPUT_HASH_MANIFEST_SHA256:
        raise RuntimeError(
            "recorded input SHA-256 mapping is stale or has been tampered with: "
            f"{actual_manifest_hash}"
        )


def validate_component_records(
    components: list[dict[str, object]],
    photometry: list[dict[str, object]],
) -> None:
    """Require the complete seven-row component-identity audit."""

    if any(set(record) != set(COMPONENT_FIELDS) for record in components):
        raise RuntimeError("component-audit fields do not match the released schema")
    record_ids = [str(record["record_id"]) for record in components]
    if len(record_ids) != len(EXPECTED_COMPONENT_RECORD_IDS) or set(record_ids) != set(
        EXPECTED_COMPONENT_RECORD_IDS
    ):
        missing = sorted(EXPECTED_COMPONENT_RECORD_IDS - set(record_ids))
        extra = sorted(set(record_ids) - EXPECTED_COMPONENT_RECORD_IDS)
        raise RuntimeError(
            "component audit must contain exactly the seven expected records; "
            f"missing={missing}, extra={extra}, rows={len(record_ids)}"
        )

    by_id = {str(record["record_id"]): record for record in components}
    if len(by_id) != len(components):
        raise RuntimeError("component audit contains a duplicate record ID")
    expected_decisions = {
        "rejected-source17-archived-anchor": ("rejected", "no"),
        "accepted-source18-fixed-canonical": ("accepted", "yes"),
        "source18-free-centroid-check": ("diagnostic", "no"),
        "source18-fixed-width-0p18": ("diagnostic", "no"),
        "source18-free-width-0p18": ("diagnostic", "no"),
        "source18-fixed-width-0p24": ("diagnostic", "no"),
        "source18-free-width-0p24": ("diagnostic", "no"),
    }
    for record_id, (decision, included) in expected_decisions.items():
        record = by_id[record_id]
        if (
            str(record["decision"]) != decision
            or str(record["included_in_canonical_table"]) != included
        ):
            raise RuntimeError(f"unexpected decision for component record {record_id}")
        _finite_float(record, "frequency_ghz", positive=True)
        _finite_float(record, "fit_region_width_arcsec", positive=True)
        _finite_float(record, "flux_mjy", positive=True)
        _finite_float(
            record,
            "image_scatter_mjy_per_beam",
            positive=True,
        )

    anchor_rows = [
        record
        for record in photometry
        if math.isclose(
            _finite_float(record, "frequency_ghz", positive=True),
            349.40,
            rel_tol=0.0,
            abs_tol=1e-6,
        )
    ]
    if len(anchor_rows) != 1:
        raise RuntimeError("expected exactly one 349.40 GHz photometry anchor")
    anchor = anchor_rows[0]
    accepted = by_id["accepted-source18-fixed-canonical"]
    for field in ("product", "frequency_ghz"):
        if str(accepted[field]) != str(anchor[field]):
            raise RuntimeError(f"canonical component has stale {field}")
    if not math.isclose(
        _finite_float(accepted, "flux_mjy", positive=True),
        _finite_float(anchor, "unrounded_flux_mjy", positive=True),
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise RuntimeError("canonical component flux is stale")
    if not math.isclose(
        _finite_float(
            accepted,
            "image_scatter_mjy_per_beam",
            positive=True,
        ),
        _finite_float(
            anchor,
            "adopted_noise_mjy_per_beam",
            positive=True,
        ),
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise RuntimeError("canonical component adopted image scatter is stale")


def preserve_pre_component_fields(
    photometry: list[dict[str, object]],
    previous_path: Path,
    correction_manifest_path: Path,
) -> None:
    """Keep exact, ledger-authenticated pre-correction fields across reruns.

    The prior photometry sidecar authenticates the 27 previously released
    source-18 products.  The correction manifest is the authoritative
    pre-application ledger for their ``old_flux_field`` and ``old_snr``.  New
    Task 28 products retain the generic-fit fields produced in the current run.
    """

    if not previous_path.is_file():
        raise RuntimeError(
            f"previous NGC 4945 photometry lineage is missing: {previous_path}"
        )
    if not correction_manifest_path.is_file():
        raise RuntimeError(
            f"continuum correction lineage is missing: {correction_manifest_path}"
        )
    with previous_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = tuple(reader.fieldnames or ())
        if fields not in {LEGACY_PHOTOMETRY_FIELDS, PHOTOMETRY_FIELDS}:
            raise RuntimeError(
                "previous NGC 4945 photometry has an unrecognized lineage schema: "
                f"{fields!r}"
            )
        previous_rows = list(reader)
    if len(previous_rows) != EXPECTED_LEGACY_ROWS:
        raise RuntimeError(
            "previous NGC 4945 photometry lineage must contain exactly "
            f"{EXPECTED_LEGACY_ROWS} rows"
        )
    previous_products = [portable_product(str(row["product"])) for row in previous_rows]
    if len(set(previous_products)) != EXPECTED_LEGACY_ROWS:
        raise RuntimeError(
            "previous NGC 4945 photometry lineage contains duplicate products"
        )
    current_products = [portable_product(str(row["product"])) for row in photometry]
    if len(current_products) != EXPECTED_ROWS or len(set(current_products)) != EXPECTED_ROWS:
        raise RuntimeError(
            f"new NGC 4945 photometry must contain exactly {EXPECTED_ROWS} unique products "
            "before preserving lineage"
        )
    if not set(previous_products) < set(current_products):
        raise RuntimeError(
            "previous NGC 4945 product inventory is not a proper subset of the expanded inventory"
        )
    with correction_manifest_path.open(encoding="utf-8", newline="") as handle:
        correction_reader = csv.DictReader(handle)
        correction_fields = tuple(correction_reader.fieldnames or ())
        if correction_fields not in {LEGACY_CORRECTION_FIELDS, CORRECTION_FIELDS}:
            raise RuntimeError(
                "continuum correction lineage has an unrecognized schema: "
                f"{correction_fields!r}"
            )
        correction_rows = [
            row for row in correction_reader if str(row["source"]) == SOURCE
        ]
    if len(correction_rows) != EXPECTED_LEGACY_ROWS:
        raise RuntimeError(
            "continuum correction lineage must contain exactly "
            f"{EXPECTED_LEGACY_ROWS} previously released NGC 4945 rows"
        )
    correction_products = [
        portable_product(str(row["product"])) for row in correction_rows
    ]
    if len(set(correction_products)) != EXPECTED_LEGACY_ROWS:
        raise RuntimeError(
            "continuum correction lineage contains duplicate NGC 4945 products"
        )
    if set(correction_products) != set(previous_products):
        raise RuntimeError(
            "previous photometry and correction lineages have different NGC 4945 products"
        )
    correction_by_product = dict(zip(correction_products, correction_rows))
    for row in photometry:
        product = portable_product(str(row["product"]))
        old = correction_by_product.get(product)
        if old is None:
            # The current generic fit is the pre-source-18 measurement for a
            # newly admitted product, so its fields are already the correct
            # old-value lineage for the new correction-manifest row.
            continue
        old_flux_field = str(old["old_flux_field"]).strip()
        old_snr = str(old["old_snr"]).strip()
        if not old_flux_field:
            raise RuntimeError(f"blank historical old_flux_field for {product}")
        try:
            old_flux = float(old_flux_field.rstrip("*"))
        except ValueError as error:
            raise RuntimeError(
                f"invalid historical old_flux_field for {product}"
            ) from error
        if not math.isfinite(old_flux) or old_flux <= 0.0:
            raise RuntimeError(f"invalid historical old_flux_field for {product}")
        if old_snr != "NA":
            try:
                old_snr_value = float(old_snr)
            except ValueError as error:
                raise RuntimeError(f"invalid historical old_snr for {product}") from error
            if not math.isfinite(old_snr_value) or old_snr_value <= 0.0:
                raise RuntimeError(f"invalid historical old_snr for {product}")
        row["old_flux_field"] = old_flux_field
        row["old_snr"] = old_snr


def apply_sorted_table(path: Path, photometry: list[dict[str, object]]) -> None:
    replacements = {str(row["product"]): row for row in photometry}
    rows = read_sorted(path)
    applied: set[str] = set()
    for row in rows:
        product = portable_product(row[7])
        correction = replacements.get(product)
        if correction is None:
            continue
        row[3] = str(correction["new_flux_field"])
        row[6] = str(correction["new_snr"])
        applied.add(product)
    if applied != set(replacements):
        missing = sorted(set(replacements) - applied)
        raise RuntimeError(f"failed to apply NGC 4945 products: {missing}")
    # Keep the released table's human-readable grouping: one row per line and
    # two blank lines between sources (including after the final source).
    output: list[str] = []
    for index, row in enumerate(rows):
        output.append(f"{row!r}\n")
        next_source = rows[index + 1][0] if index + 1 < len(rows) else None
        if next_source != row[0]:
            output.append("\n\n")
    path.write_text("".join(output), encoding="utf-8")


def merge_correction_manifest(
    path: Path,
    photometry: list[dict[str, object]],
) -> None:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = tuple(reader.fieldnames or ())
        existing = list(reader)
    if not existing:
        raise RuntimeError("continuum correction manifest cannot be empty")
    if fields == LEGACY_CORRECTION_FIELDS:
        existing = [
            {
                "schema_version": CORRECTION_SCHEMA_VERSION,
                **{
                    field: row[field]
                    for field in LEGACY_CORRECTION_FIELDS[:-1]
                },
                "validated_adopted_noise_snr": row[
                    "validated_background_only_snr"
                ],
            }
            for row in existing
        ]
    elif fields == CORRECTION_FIELDS:
        if any(
            row["schema_version"] != CORRECTION_SCHEMA_VERSION
            for row in existing
        ):
            raise RuntimeError("continuum correction manifest has a stale schema version")
    else:
        raise RuntimeError(
            f"continuum correction manifest has unexpected fields {fields!r}"
        )
    by_product = {str(row["product"]): row for row in existing}
    if len(by_product) != len(existing):
        raise RuntimeError("continuum correction manifest contains duplicate products")

    for measurement in photometry:
        product = str(measurement["product"])
        old = by_product.get(product)
        if old is not None and (
            old["old_flux_field"] != str(measurement["old_flux_field"])
            or old["old_snr"] != str(measurement["old_snr"])
        ):
            raise RuntimeError(
                f"correction-manifest historical lineage disagrees with photometry "
                f"for {product}"
            )
        reasons = (
            "fixed_nuclear_component_source18",
            "single_local_scaled_mad",
            "source18_10bmaj_local_noise_aperture",
        )
        by_product[product] = {
            "schema_version": CORRECTION_SCHEMA_VERSION,
            "source": SOURCE,
            "telescope": measurement["telescope"],
            "frequency_ghz": measurement["frequency_ghz"],
            "date_label": measurement["date_label"],
            "product": product,
            "old_flux_field": measurement["old_flux_field"],
            "new_flux_field": measurement["new_flux_field"],
            "old_snr": measurement["old_snr"],
            "new_snr": measurement["new_snr"],
            "correction_reasons": "+".join(reasons),
            "validated_unrounded_flux_mjy": measurement["unrounded_flux_mjy"],
            "validated_adopted_noise_snr": measurement["new_snr"],
        }

    def sort_key(row: dict[str, object]) -> tuple[object, ...]:
        try:
            frequency: object = float(str(row["frequency_ghz"]))
        except ValueError:
            frequency = str(row["frequency_ghz"])
        return (
            str(row["source"]),
            str(row["telescope"]),
            frequency,
            str(row["date_label"]),
            str(row["product"]),
        )

    write_csv(path, CORRECTION_FIELDS, sorted(by_product.values(), key=sort_key))


def check_recorded(
    source_table: Path,
    photometry_path: Path,
    component_path: Path,
    correction_manifest_path: Path,
) -> None:
    photometry = read_exact_csv(
        photometry_path,
        PHOTOMETRY_FIELDS,
        label="NGC 4945 photometry sidecar",
    )
    validate_photometry_records(photometry)

    canonical = {
        portable_product(row[7]): row
        for row in ngc4945_rows(read_sorted(source_table))
    }
    recorded_products = {str(record["product"]) for record in photometry}
    if set(canonical) != recorded_products:
        raise RuntimeError(
            "canonical table and NGC 4945 photometry product inventories differ"
        )
    for record in photometry:
        row = canonical[record["product"]]
        if str(row[3]) != record["new_flux_field"] or str(row[6]) != record["new_snr"]:
            raise RuntimeError(f"canonical row is stale for {record['product']}")

    components = read_exact_csv(
        component_path,
        COMPONENT_FIELDS,
        label="NGC 4945 component audit",
    )
    validate_component_records(components, photometry)

    all_corrections = read_exact_csv(
        correction_manifest_path,
        CORRECTION_FIELDS,
        label="continuum correction manifest",
    )
    if any(
        row["schema_version"] != CORRECTION_SCHEMA_VERSION
        for row in all_corrections
    ):
        raise RuntimeError("continuum correction manifest has a stale schema version")
    correction_rows = [row for row in all_corrections if row["source"] == SOURCE]
    corrections = {row["product"]: row for row in correction_rows}
    if len(corrections) != len(correction_rows):
        raise RuntimeError("correction manifest contains a duplicate NGC 4945 product")
    if set(corrections) != {row["product"] for row in photometry}:
        raise RuntimeError("correction manifest does not contain exactly the NGC 4945 rerun")
    for record in photometry:
        correction = corrections[record["product"]]
        reasons = correction["correction_reasons"].split("+")
        expected_reasons = {
            "fixed_nuclear_component_source18",
            "single_local_scaled_mad",
            "source18_10bmaj_local_noise_aperture",
        }
        if set(reasons) != expected_reasons or len(reasons) != len(expected_reasons):
            raise RuntimeError(
                f"correction reasons are not the exact v2 set for {record['product']}"
            )
        if "fixed_nuclear_component_source18" not in reasons:
            raise RuntimeError(
                f"correction reason is missing for {record['product']}"
            )
        if "source18_10bmaj_local_noise_aperture" not in reasons:
            raise RuntimeError(
                f"local-noise correction reason is missing for {record['product']}"
            )
        if "single_local_scaled_mad" not in reasons:
            raise RuntimeError(
                f"local-scaled-MAD correction reason is missing for {record['product']}"
            )
        if "single_background_rms" in reasons or "single_local_image_scatter" in reasons:
            raise RuntimeError(f"stale correction reason is present for {record['product']}")
        if (
            correction["old_flux_field"] != record["old_flux_field"]
            or correction["old_snr"] != record["old_snr"]
            or correction["new_flux_field"] != record["new_flux_field"]
            or correction["new_snr"] != record["new_snr"]
            or correction["validated_adopted_noise_snr"] != record["new_snr"]
        ):
            raise RuntimeError(
                f"correction manifest is stale for {record['product']}"
            )
        if not math.isclose(
            _finite_float(correction, "validated_unrounded_flux_mjy", positive=True),
            _finite_float(record, "unrounded_flux_mjy", positive=True),
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise RuntimeError(
                f"correction manifest unrounded flux is stale for {record['product']}"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-root", type=Path)
    parser.add_argument("--source-table", type=Path, default=DEFAULT_SORTED)
    parser.add_argument("--output", type=Path, default=DEFAULT_PHOTOMETRY)
    parser.add_argument("--component-output", type=Path, default=DEFAULT_COMPONENT_AUDIT)
    parser.add_argument("--correction-manifest", type=Path, default=DEFAULT_CORRECTIONS)
    parser.add_argument(
        "--apply",
        action="store_true",
        help="apply new flux/SNR values to the sorted table and correction manifest",
    )
    parser.add_argument(
        "--apply-recorded",
        action="store_true",
        help=(
            "reapply the checked-in fixed-nucleus sidecar to a newly combined "
            "sorted table without rerunning CASA"
        ),
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="validate checked-in sidecars and canonical rows without CASA/FITS data",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.apply and args.apply_recorded:
        raise SystemExit("--apply and --apply-recorded are mutually exclusive")
    if args.check:
        check_recorded(
            args.source_table,
            args.output,
            args.component_output,
            args.correction_manifest,
        )
        print("NGC 4945 fixed-nucleus records and canonical rows are consistent")
        return
    if args.apply_recorded:
        photometry = read_exact_csv(
            args.output,
            PHOTOMETRY_FIELDS,
            label="NGC 4945 photometry sidecar",
        )
        components = read_exact_csv(
            args.component_output,
            COMPONENT_FIELDS,
            label="NGC 4945 component audit",
        )
        # Validate immutable inputs before mutating the combined sorted table.
        validate_photometry_records(photometry)
        validate_component_records(components, photometry)
        apply_sorted_table(args.source_table, photometry)
        merge_correction_manifest(args.correction_manifest, photometry)
        check_recorded(
            args.source_table,
            args.output,
            args.component_output,
            args.correction_manifest,
        )
        print(f"reapplied {len(photometry)} recorded fixed-nucleus measurements")
        return
    if args.archive_root is None:
        raise SystemExit(
            "--archive-root is required unless --check or --apply-recorded is used"
        )

    photometry, components = build_records(args.source_table, args.archive_root)
    preserve_pre_component_fields(
        photometry,
        args.output,
        args.correction_manifest,
    )
    validate_photometry_records(photometry)
    validate_component_records(components, photometry)
    write_csv(args.output, PHOTOMETRY_FIELDS, photometry)
    write_csv(args.component_output, COMPONENT_FIELDS, components)
    if args.apply:
        apply_sorted_table(args.source_table, photometry)
        merge_correction_manifest(args.correction_manifest, photometry)
    print(f"wrote {len(photometry)} fixed-nucleus measurements to {args.output}")
    print(f"wrote {len(components)} component-audit records to {args.component_output}")


if __name__ == "__main__":
    main()
