"""Fail-closed validation for final continuum noise-audit inventories."""

from __future__ import annotations

import ast
import csv
import math
from pathlib import Path
from typing import Any

from continuum_measurements import (
    FINAL_NOISE_APERTURE_BMAJ_MULTIPLIER,
    FINAL_NOISE_MIN_INDEPENDENT_BEAMS,
    THREE_SIGMA_FROM_FWHM,
    classify_ordinary_fit,
    effective_independent_noise_samples,
    final_noise_review_flags,
    format_two_decimal,
    gaussian_beam_area_arcsec2,
    independent_beam_areas,
    portable_release_product,
    require_fixed_restoring_beam_shape,
    square_aperture_coverage,
)
from continuum_units import brightness_to_mjy_per_beam_factor
from continuum_product_additions import (
    addition_products,
    load_product_addition_policy,
    publication_addition_products,
)
from continuum_boundary_rejections import (
    boundary_rejection_products,
    load_boundary_rejection_policy,
    validate_boundary_rejection_record_identity,
)


GENERIC_FINAL_NOISE_AUDIT_FIELDS = (
    "schema_version",
    "run_context",
    "source",
    "telescope",
    "frequency_ghz",
    "product",
    "forced_policy",
    "measurement_scope",
    "fit_selection_policy",
    "aperture_center_policy",
    "aperture_center_x_pixel",
    "aperture_center_y_pixel",
    "fit_center_reference_frame",
    "fit_center_ra_deg",
    "fit_center_dec_deg",
    "aperture_width_arcsec",
    "aperture_pre_exclusion_npts",
    "aperture_selected_pixel_area_over_requested_area",
    "mask_semantics",
    "mask_major_semiaxis_arcsec",
    "mask_minor_semiaxis_arcsec",
    "mask_pa_deg",
    "beam_major_arcsec",
    "beam_minor_arcsec",
    "beam_pa_deg",
    "fit_beam_major_delta_arcsec",
    "fit_beam_minor_delta_arcsec",
    "fit_beam_pa_delta_deg",
    "fit_beam_pa_guard_applicable",
    "image_width_pixels",
    "image_height_pixels",
    "aperture_half_width_x_pixels",
    "aperture_half_width_y_pixels",
    "aperture_left_clearance_pixels",
    "aperture_right_clearance_pixels",
    "aperture_bottom_clearance_pixels",
    "aperture_top_clearance_pixels",
    "aperture_minimum_edge_clearance_pixels",
    "aperture_edge_status",
    "npts",
    "pixel_width_arcsec",
    "pixel_height_arcsec",
    "pixel_area_arcsec2",
    "valid_unmasked_area_arcsec2",
    "gaussian_beam_area_arcsec2",
    "independent_beam_areas",
    "effective_independent_samples",
    "minimum_independent_beams",
    "brightness_unit",
    "rms_mjy_per_beam",
    "sigma_mjy_per_beam",
    "mean_mjy_per_beam",
    "scaled_mad_mjy_per_beam",
    "adopted_noise_estimator",
    "adopted_noise_mjy_per_beam",
    "mean_over_sigma",
    "rms_over_scaled_mad",
    "review_flags",
    "review_disposition",
    "flux_mjy",
    "snr",
    "classification",
    "ordinary_upper_limit_mjy",
    "reported_limit_mjy",
    "reported_flux_field",
    "row_disposition",
    "guard_status",
)

EXPECTED_PUBLICATION_ROWS = 406
EXPECTED_ORDINARY_ROWS = 394
EXPECTED_FORCED_DIRTY_SMA_ROWS = 12
EXPECTED_NGC4945_ROWS = 45
NGC4945 = "NGC4945"
SOURCE18_COMPONENT_IDENTITY = "Emig2020_source18_AGN_maser_nucleus"
SOURCE17_COMPONENT_IDENTITY = "Emig2020_source17_SSC_SNR_associated"
FORCED_FIT_POLICY_CONVERGED = "reviewed_forced_dirty_sma_fit_no_retry"
FORCED_FIT_POLICY_NONCONVERGED = (
    "reviewed_forced_dirty_sma_nonconverged_fit_no_retry"
)
_FORCED_SMA_PREFIX = (
    "analysis/sma_calibration/230316_128/images/nitercal_0/"
    "flag2_cal1/fields/science/"
)
EXPECTED_FORCED_FIT_POLICY_BY_PRODUCT = {
    f"{_FORCED_SMA_PREFIX}{basename}_test.image": (
        FORCED_FIT_POLICY_NONCONVERGED
        if basename in {"NGC6264", "UGC6093"}
        else FORCED_FIT_POLICY_CONVERGED
    )
    for basename in (
        "ESO558-G009",
        "MRK1029",
        "NGC1194",
        "NGC1320",
        "NGC2273",
        "NGC2960",
        "NGC4388",
        "NGC5765b",
        "NGC6264",
        "NGC6323",
        "UGC3789",
        "UGC6093",
    )
}


def _read_exact_csv(
    path: Path,
    fields: tuple[str, ...],
    *,
    label: str,
) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != fields:
            raise RuntimeError(
                f"{label} fields differ from the exact schema: "
                f"{tuple(reader.fieldnames or ())!r}"
            )
        rows = list(reader)
    if any(any(row.get(field) is None for field in fields) for row in rows):
        raise RuntimeError(f"{label} contains an incomplete row")
    return rows


def _read_publication_rows(path: Path) -> list[list[Any]]:
    rows: list[list[Any]] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            row = ast.literal_eval(line)
        except (SyntaxError, ValueError) as error:
            raise RuntimeError(
                f"invalid publication row at {path}:{line_number}"
            ) from error
        if not isinstance(row, list) or len(row) != 8:
            raise RuntimeError(
                f"publication row at {path}:{line_number} is not an eight-field list"
            )
        row[7] = portable_release_product(row[7])
        rows.append(row)
    products = [str(row[7]) for row in rows]
    if len(products) != len(set(products)):
        raise RuntimeError("publication table contains duplicate product identities")
    return rows


def _finite(record: dict[str, str], field: str, *, positive: bool = False) -> float:
    try:
        value = float(record[field])
    except (KeyError, TypeError, ValueError) as error:
        raise RuntimeError(f"invalid {field} for {record.get('product', '<unknown>')}") from error
    if not math.isfinite(value) or (positive and value <= 0.0):
        raise RuntimeError(f"non-finite or invalid {field} for {record.get('product')}")
    return value


def _require_close(
    product: str,
    label: str,
    actual: float,
    expected: float,
) -> None:
    if not math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12):
        raise RuntimeError(
            f"stale generic {label} for {product}: {actual} vs {expected}"
        )


def _validate_ordinary_noise_math(record: dict[str, str]) -> None:
    """Recompute every ordinary generic geometry/statistic audit field."""

    product = record["product"]
    if (
        record["aperture_center_policy"] != "fitted_component_pixelcoords"
        or record["mask_semantics"] != "crtf_ellipse_axes_are_semiaxes"
        or record["fit_selection_policy"]
        != "single_initial_fixed_beam_fit_no_widened_retry"
    ):
        raise RuntimeError(f"ordinary aperture/mask policy is stale for {product}")
    center_x = _finite(record, "aperture_center_x_pixel")
    center_y = _finite(record, "aperture_center_y_pixel")
    fit_center_frame = record["fit_center_reference_frame"]
    if fit_center_frame not in {"ICRS", "J2000"}:
        raise RuntimeError(
            f"ordinary fit-center reference frame is invalid for {product}"
        )
    fit_center_ra = _finite(record, "fit_center_ra_deg")
    fit_center_dec = _finite(record, "fit_center_dec_deg")
    if not 0.0 <= fit_center_ra < 360.0:
        raise RuntimeError(f"ordinary fit-center RA is out of range for {product}")
    if not -90.0 <= fit_center_dec <= 90.0:
        raise RuntimeError(f"ordinary fit-center Dec is out of range for {product}")
    beam_major = _finite(record, "beam_major_arcsec", positive=True)
    beam_minor = _finite(record, "beam_minor_arcsec", positive=True)
    beam_pa = _finite(record, "beam_pa_deg")
    aperture_width = _finite(record, "aperture_width_arcsec", positive=True)
    aperture_pre_npts = _finite(record, "aperture_pre_exclusion_npts", positive=True)
    npts = _finite(record, "npts", positive=True)
    if not aperture_pre_npts.is_integer() or not npts.is_integer():
        raise RuntimeError(f"ordinary aperture npts fields must be integers for {product}")
    if npts > aperture_pre_npts:
        raise RuntimeError(f"source-excluded npts exceeds aperture npts for {product}")
    pixel_width = _finite(record, "pixel_width_arcsec", positive=True)
    pixel_height = _finite(record, "pixel_height_arcsec", positive=True)
    pixel_area = _finite(record, "pixel_area_arcsec2", positive=True)
    image_width = _finite(record, "image_width_pixels", positive=True)
    image_height = _finite(record, "image_height_pixels", positive=True)
    if not image_width.is_integer() or not image_height.is_integer():
        raise RuntimeError(f"ordinary image dimensions must be integers for {product}")

    expected_aperture_width = beam_major * FINAL_NOISE_APERTURE_BMAJ_MULTIPLIER
    _require_close(product, "aperture width", aperture_width, expected_aperture_width)
    _require_close(product, "pixel area", pixel_area, pixel_width * pixel_height)
    coverage = square_aperture_coverage(
        center_x,
        center_y,
        aperture_width,
        pixel_width,
        pixel_height,
        image_width,
        image_height,
    )
    coverage_values = {
        "aperture half width x": (
            "aperture_half_width_x_pixels",
            coverage.half_width_x_pixels,
        ),
        "aperture half width y": (
            "aperture_half_width_y_pixels",
            coverage.half_width_y_pixels,
        ),
        "aperture left clearance": (
            "aperture_left_clearance_pixels",
            coverage.left_clearance_pixels,
        ),
        "aperture right clearance": (
            "aperture_right_clearance_pixels",
            coverage.right_clearance_pixels,
        ),
        "aperture bottom clearance": (
            "aperture_bottom_clearance_pixels",
            coverage.bottom_clearance_pixels,
        ),
        "aperture top clearance": (
            "aperture_top_clearance_pixels",
            coverage.top_clearance_pixels,
        ),
        "aperture minimum clearance": (
            "aperture_minimum_edge_clearance_pixels",
            coverage.minimum_edge_clearance_pixels,
        ),
    }
    for label, (field, expected) in coverage_values.items():
        _require_close(product, label, _finite(record, field), expected)
    expected_edge_status = (
        "truncated_available_valid_intersection"
        if coverage.minimum_edge_clearance_pixels < -1e-9
        else "fully_contained"
    )
    if record["aperture_edge_status"] != expected_edge_status:
        raise RuntimeError(f"stale generic aperture edge status for {product}")

    expected_valid_fraction = (
        aperture_pre_npts * pixel_area / (aperture_width * aperture_width)
    )
    _require_close(
        product,
        "selected-pixel/requested-aperture area ratio",
        _finite(
            record,
            "aperture_selected_pixel_area_over_requested_area",
            positive=True,
        ),
        expected_valid_fraction,
    )
    valid_area = npts * pixel_area
    beam_area = gaussian_beam_area_arcsec2(beam_major, beam_minor)
    beam_areas = independent_beam_areas(
        npts,
        pixel_width,
        pixel_height,
        beam_major,
        beam_minor,
    )
    effective_samples = effective_independent_noise_samples(
        npts,
        pixel_width,
        pixel_height,
        beam_major,
        beam_minor,
    )
    area_values = {
        "valid unmasked area": ("valid_unmasked_area_arcsec2", valid_area),
        "Gaussian beam area": ("gaussian_beam_area_arcsec2", beam_area),
        "independent beam areas": ("independent_beam_areas", beam_areas),
        "effective independent samples": (
            "effective_independent_samples",
            effective_samples,
        ),
        "minimum independent beams": (
            "minimum_independent_beams",
            FINAL_NOISE_MIN_INDEPENDENT_BEAMS,
        ),
    }
    for label, (field, expected) in area_values.items():
        _require_close(product, label, _finite(record, field), expected)
    if effective_samples < FINAL_NOISE_MIN_INDEPENDENT_BEAMS:
        raise RuntimeError(f"ordinary final-noise N_eff is below 50 for {product}")

    mask_major = _finite(record, "mask_major_semiaxis_arcsec", positive=True)
    mask_minor = _finite(record, "mask_minor_semiaxis_arcsec", positive=True)
    mask_pa = _finite(record, "mask_pa_deg")
    fixed_shape = require_fixed_restoring_beam_shape(
        mask_major / THREE_SIGMA_FROM_FWHM,
        mask_minor / THREE_SIGMA_FROM_FWHM,
        mask_pa,
        beam_major,
        beam_minor,
        beam_pa,
    )
    shape_values = {
        "fit-beam major delta": (
            "fit_beam_major_delta_arcsec",
            fixed_shape.major_delta_arcsec,
        ),
        "fit-beam minor delta": (
            "fit_beam_minor_delta_arcsec",
            fixed_shape.minor_delta_arcsec,
        ),
        "fit-beam PA delta": (
            "fit_beam_pa_delta_deg",
            fixed_shape.pa_delta_deg,
        ),
    }
    for label, (field, expected) in shape_values.items():
        _require_close(product, label, _finite(record, field), expected)
    expected_pa_guard = "yes" if fixed_shape.pa_guard_applicable else "no"
    if record["fit_beam_pa_guard_applicable"] != expected_pa_guard:
        raise RuntimeError(f"stale generic fit-beam PA applicability for {product}")

    brightness_to_mjy_per_beam_factor(record["brightness_unit"])
    rms = _finite(record, "rms_mjy_per_beam", positive=True)
    sigma = _finite(record, "sigma_mjy_per_beam", positive=True)
    mean = _finite(record, "mean_mjy_per_beam")
    scaled_mad = _finite(record, "scaled_mad_mjy_per_beam", positive=True)
    if record["adopted_noise_estimator"] != "scaled_mad":
        raise RuntimeError(f"ordinary adopted-noise estimator is stale for {product}")
    _require_close(
        product,
        "adopted scaled-MAD noise",
        _finite(record, "adopted_noise_mjy_per_beam", positive=True),
        scaled_mad,
    )
    _require_close(
        product,
        "mean/sigma ratio",
        _finite(record, "mean_over_sigma"),
        abs(mean) / sigma,
    )
    _require_close(
        product,
        "RMS/scaled-MAD ratio",
        _finite(record, "rms_over_scaled_mad", positive=True),
        rms / scaled_mad,
    )
    expected_flags = "+".join(
        final_noise_review_flags(
            mean,
            sigma,
            rms,
            scaled_mad,
            aperture_minimum_edge_clearance_pixels=(
                coverage.minimum_edge_clearance_pixels
            ),
        )
    ) or "none"
    if record["review_flags"] != expected_flags:
        raise RuntimeError(f"stale generic review flags for {product}")
    expected_review_disposition = (
        "diagnostic_flags_recorded_not_classification_overrides"
        if expected_flags != "none"
        else "no_diagnostic_flags"
    )
    if record["review_disposition"] != expected_review_disposition:
        raise RuntimeError(f"stale generic review disposition for {product}")


def _validate_forced_noise_math(record: dict[str, str]) -> float:
    """Recompute the retained dirty-SMA one-beam RMS diagnostics."""

    product = record["product"]
    expected_text = {
        "aperture_center_policy": "catalog_or_anchor_coordinate",
        "aperture_width_arcsec": "one_restoring_beam_rotbox",
        "mask_semantics": "none_source_not_excluded_reviewed_forced_policy",
        "minimum_independent_beams": "not_applicable",
    }
    for field, expected in expected_text.items():
        if record[field] != expected:
            raise RuntimeError(f"forced dirty-SMA {field} is stale for {product}")
    expected_blank = (
        "aperture_center_x_pixel",
        "aperture_center_y_pixel",
        "fit_center_reference_frame",
        "fit_center_ra_deg",
        "fit_center_dec_deg",
        "aperture_pre_exclusion_npts",
        "aperture_selected_pixel_area_over_requested_area",
        "mask_major_semiaxis_arcsec",
        "mask_minor_semiaxis_arcsec",
        "mask_pa_deg",
        "fit_beam_major_delta_arcsec",
        "fit_beam_minor_delta_arcsec",
        "fit_beam_pa_delta_deg",
        "fit_beam_pa_guard_applicable",
        "image_width_pixels",
        "image_height_pixels",
        "aperture_half_width_x_pixels",
        "aperture_half_width_y_pixels",
        "aperture_left_clearance_pixels",
        "aperture_right_clearance_pixels",
        "aperture_bottom_clearance_pixels",
        "aperture_top_clearance_pixels",
        "aperture_minimum_edge_clearance_pixels",
        "aperture_edge_status",
        "mean_over_sigma",
        "rms_over_scaled_mad",
        "snr",
        "ordinary_upper_limit_mjy",
    )
    stale_nonblank = [field for field in expected_blank if record[field] != ""]
    if stale_nonblank:
        raise RuntimeError(
            f"forced dirty-SMA not-applicable fields are nonblank for {product}: "
            f"{stale_nonblank}"
        )

    beam_major = _finite(record, "beam_major_arcsec", positive=True)
    beam_minor = _finite(record, "beam_minor_arcsec", positive=True)
    _finite(record, "beam_pa_deg")
    npts = _finite(record, "npts", positive=True)
    if not npts.is_integer():
        raise RuntimeError(f"forced dirty-SMA npts is not an integer for {product}")
    pixel_width = _finite(record, "pixel_width_arcsec", positive=True)
    pixel_height = _finite(record, "pixel_height_arcsec", positive=True)
    pixel_area = _finite(record, "pixel_area_arcsec2", positive=True)
    valid_area = _finite(record, "valid_unmasked_area_arcsec2", positive=True)
    beam_area = _finite(record, "gaussian_beam_area_arcsec2", positive=True)
    beam_areas = _finite(record, "independent_beam_areas", positive=True)
    effective_samples = _finite(
        record, "effective_independent_samples", positive=True
    )
    expected_values = {
        "pixel area": (pixel_area, pixel_width * pixel_height),
        "valid unmasked area": (valid_area, npts * pixel_area),
        "Gaussian beam area": (
            beam_area,
            gaussian_beam_area_arcsec2(beam_major, beam_minor),
        ),
        "independent beam areas": (
            beam_areas,
            independent_beam_areas(
                npts, pixel_width, pixel_height, beam_major, beam_minor
            ),
        ),
        "effective independent samples": (
            effective_samples,
            effective_independent_noise_samples(
                npts, pixel_width, pixel_height, beam_major, beam_minor
            ),
        ),
    }
    for label, (actual, expected) in expected_values.items():
        _require_close(product, f"forced {label}", actual, expected)

    brightness_to_mjy_per_beam_factor(record["brightness_unit"])
    diagnostic_rms = _finite(record, "rms_mjy_per_beam", positive=True)
    _finite(record, "sigma_mjy_per_beam", positive=True)
    _finite(record, "mean_mjy_per_beam")
    scaled_mad = _finite(record, "scaled_mad_mjy_per_beam")
    if scaled_mad < 0.0:
        raise RuntimeError(f"forced dirty-SMA scaled MAD is negative for {product}")
    return diagnostic_rms


def _validate_generic_record(
    record: dict[str, str],
    table_row: list[Any] | None,
    *,
    expected_disposition: str | None = None,
) -> str:
    product = str(record["product"])
    if record["schema_version"] != "1" or not record["run_context"]:
        raise RuntimeError(f"invalid generic audit schema/run context for {product}")
    if table_row is not None:
        if (
            record["source"] != str(table_row[0])
            or record["telescope"] != str(table_row[1])
            or not math.isclose(
                _finite(record, "frequency_ghz", positive=True),
                float(table_row[2]),
                rel_tol=0.0,
                abs_tol=5e-3,
            )
        ):
            raise RuntimeError(
                f"generic audit identity fields disagree with table for {product}"
            )
    replaces_generic_component = record["source"] == NGC4945
    if (
        table_row is not None
        and not replaces_generic_component
        and record["reported_flux_field"] != str(table_row[3])
    ):
        raise RuntimeError(f"generic audit flux field disagrees with table for {product}")

    if expected_disposition is None:
        expected_disposition = (
            "rejected_source17_requires_source18_join"
            if record["source"] == NGC4945
            else "generic_table_retained"
        )
    if record["row_disposition"] != expected_disposition:
        raise RuntimeError(
            f"nonterminal or incorrect generic disposition for {product}: "
            f"{record['row_disposition']!r}"
        )

    if record["forced_policy"] == "no":
        flux = _finite(record, "flux_mjy")
        if record["measurement_scope"] != "publication_final_10bmaj_source_excluded":
            raise RuntimeError(f"ordinary audit has the wrong scope for {product}")
        if record["guard_status"] != "passed":
            raise RuntimeError(f"ordinary final-noise guard did not pass for {product}")
        _validate_ordinary_noise_math(record)
        adopted_noise = _finite(record, "adopted_noise_mjy_per_beam", positive=True)
        classification = classify_ordinary_fit(flux, adopted_noise)
        expected_class = "detection" if classification.is_detection else "ordinary_upper_limit"
        if record["classification"] != expected_class:
            raise RuntimeError(f"ordinary classification is stale for {product}")
        if not math.isclose(
            _finite(record, "snr"),
            classification.snr,
            rel_tol=1e-12,
            abs_tol=1e-12,
        ):
            raise RuntimeError(f"ordinary S/N is stale for {product}")
        if not math.isclose(
            _finite(record, "ordinary_upper_limit_mjy", positive=True),
            classification.upper_limit_mjy,
            rel_tol=1e-12,
            abs_tol=1e-12,
        ):
            raise RuntimeError(f"ordinary upper-limit invariant is stale for {product}")
        expected_table_snr = f"{classification.snr:.3f}" if classification.is_detection else "NA"
        if (
            table_row is not None
            and not replaces_generic_component
            and str(table_row[6]) != expected_table_snr
        ):
            raise RuntimeError(f"ordinary table S/N is stale for {product}")
        expected_field = (
            format_two_decimal(flux)
            if classification.is_detection
            else f"{format_two_decimal(classification.upper_limit_mjy)}*"
        )
        if record["reported_flux_field"] != expected_field:
            raise RuntimeError(f"ordinary reported flux/limit is stale for {product}")
        if classification.is_detection:
            if record["reported_limit_mjy"] != "":
                raise RuntimeError(f"detection carries a reported limit for {product}")
        elif not math.isclose(
            _finite(record, "reported_limit_mjy", positive=True),
            classification.upper_limit_mjy,
            rel_tol=1e-12,
            abs_tol=1e-12,
        ):
            raise RuntimeError(f"ordinary reported limit is stale for {product}")
        return "ordinary"

    if record["forced_policy"] != "yes_dirty_sma_nitercal_0":
        raise RuntimeError(f"unknown forced policy for {product}")
    forced_fit_policy = record["fit_selection_policy"]
    if forced_fit_policy == FORCED_FIT_POLICY_CONVERGED:
        _finite(record, "flux_mjy")
    elif forced_fit_policy == FORCED_FIT_POLICY_NONCONVERGED:
        if record["flux_mjy"] != "":
            raise RuntimeError(
                f"nonconverged forced dirty-SMA fit carries a fitted flux for {product}"
            )
    else:
        raise RuntimeError(
            f"forced dirty-SMA fit-selection policy is invalid for {product}"
        )
    diagnostic_rms = _validate_forced_noise_math(record)
    if (
        record["telescope"] != "SMA"
        or "/nitercal_0/" not in f"/{product}"
        or record["measurement_scope"] != "reviewed_forced_dirty_sma_one_beam_region"
        or record["classification"] != "forced_upper_limit"
        or record["guard_status"] != "not_applicable_reviewed_forced_policy"
        or record["review_flags"] != "reviewed_forced_policy"
        or record["review_disposition"] != "preapproved_forced_policy"
        or record["adopted_noise_estimator"]
        != "classic_rms_reviewed_forced_policy"
        or record["snr"] not in {"", "NA"}
        or (table_row is not None and str(table_row[6]) != "NA")
    ):
        raise RuntimeError(f"forced dirty-SMA policy fields are invalid for {product}")
    _require_close(
        product,
        "forced adopted RMS",
        _finite(record, "adopted_noise_mjy_per_beam", positive=True),
        diagnostic_rms,
    )
    limit = _finite(record, "reported_limit_mjy", positive=True)
    if not math.isclose(
        limit, 3.0 * diagnostic_rms, rel_tol=1e-12, abs_tol=1e-12
    ):
        raise RuntimeError(f"forced dirty-SMA limit is not exactly 3 RMS for {product}")
    if record["reported_flux_field"] != f"{format_two_decimal(limit)}*":
        raise RuntimeError(f"forced dirty-SMA display limit is stale for {product}")
    if (
        table_row is not None
        and not replaces_generic_component
        and record["reported_flux_field"] != str(table_row[3])
    ):
        raise RuntimeError(f"forced dirty-SMA table limit is stale for {product}")
    return "forced"


def validate_final_noise_inventory(
    publication_table_path: Path,
    generic_audit_path: Path,
    source18_photometry_path: Path,
    addition_policy_path: Path | None = None,
    boundary_rejection_policy_path: Path | None = None,
    analysis_root: Path | None = None,
) -> dict[str, int]:
    """Validate retained rows plus the exact reviewed boundary diagnostics."""

    addition_rows = (
        []
        if addition_policy_path is None
        else load_product_addition_policy(addition_policy_path)
    )
    all_additions = addition_products(addition_rows)
    publication_additions = publication_addition_products(addition_rows)
    boundary_rows = (
        []
        if boundary_rejection_policy_path is None
        else load_boundary_rejection_policy(
            boundary_rejection_policy_path,
            analysis_root=analysis_root,
        )
    )
    boundary_by_product = {row["product"]: row for row in boundary_rows}
    boundary_products = boundary_rejection_products(boundary_rows)
    expected_publication_rows = EXPECTED_PUBLICATION_ROWS + len(
        publication_additions
    )
    publication_rows = _read_publication_rows(publication_table_path)
    if len(publication_rows) != expected_publication_rows:
        raise RuntimeError(
            f"expected {expected_publication_rows} publication rows; "
            f"found {len(publication_rows)}"
        )
    table_by_product = {str(row[7]): row for row in publication_rows}
    if set(table_by_product) & all_additions != publication_additions:
        raise RuntimeError(
            "publication addition inventory differs from the exact disposition policy"
        )

    generic = _read_exact_csv(
        generic_audit_path,
        GENERIC_FINAL_NOISE_AUDIT_FIELDS,
        label="generic final-noise audit",
    )
    generic_products = [portable_release_product(row["product"]) for row in generic]
    if generic_products != sorted(generic_products):
        raise RuntimeError("generic final-noise audit is not deterministically product-sorted")
    if len(generic_products) != len(set(generic_products)):
        raise RuntimeError("generic final-noise audit contains duplicate products")
    if set(table_by_product) & boundary_products:
        raise RuntimeError(
            "boundary-rejection products must be absent from the publication table"
        )
    expected_generic_products = set(table_by_product) | boundary_products
    if set(generic_products) != expected_generic_products:
        missing = sorted(expected_generic_products - set(generic_products))
        extra = sorted(set(generic_products) - expected_generic_products)
        raise RuntimeError(
            f"generic audit/table product inventories differ: missing={missing}, extra={extra}"
        )

    ordinary = 0
    ordinary_publication = 0
    forced = 0
    boundary_rejected = 0
    generic_ngc4945: set[str] = set()
    for record, product in zip(generic, generic_products):
        record["product"] = product
        if product in boundary_by_product:
            validate_boundary_rejection_record_identity(
                record, boundary_by_product[product]
            )
            kind = _validate_generic_record(
                record,
                None,
                expected_disposition="rejected_on_boundary",
            )
            if kind != "ordinary":
                raise RuntimeError(
                    f"boundary-rejection diagnostic is not ordinary for {product}"
                )
            boundary_rejected += 1
        else:
            kind = _validate_generic_record(record, table_by_product[product])
            ordinary_publication += kind == "ordinary"
        ordinary += kind == "ordinary"
        forced += kind == "forced"
        if record["source"] == NGC4945:
            generic_ngc4945.add(product)
    expected_ordinary = EXPECTED_ORDINARY_ROWS + len(publication_additions)
    if (
        ordinary_publication != expected_ordinary
        or ordinary != expected_ordinary + len(boundary_rows)
        or boundary_rejected != len(boundary_rows)
        or forced != EXPECTED_FORCED_DIRTY_SMA_ROWS
    ):
        raise RuntimeError(
            f"expected {expected_ordinary} publication ordinary, "
            f"{len(boundary_rows)} boundary-rejected ordinary, and "
            f"{EXPECTED_FORCED_DIRTY_SMA_ROWS} forced rows; found "
            f"{ordinary_publication}/{boundary_rejected}/{forced}"
        )
    expected_forced = {
        product
        for product, row in table_by_product.items()
        if str(row[1]) == "SMA" and "/nitercal_0/" in f"/{product}"
    }
    actual_forced = {
        row["product"] for row in generic if row["forced_policy"] != "no"
    }
    if actual_forced != expected_forced or len(expected_forced) != EXPECTED_FORCED_DIRTY_SMA_ROWS:
        raise RuntimeError("forced dirty-SMA product inventory is not the exact reviewed set")
    actual_forced_fit_policies = {
        row["product"]: row["fit_selection_policy"]
        for row in generic
        if row["forced_policy"] != "no"
    }
    if actual_forced_fit_policies != EXPECTED_FORCED_FIT_POLICY_BY_PRODUCT:
        raise RuntimeError(
            "forced dirty-SMA convergence inventory differs from the exact "
            "reviewed 10-converged/2-nonconverged set"
        )
    if len(generic_ngc4945) != EXPECTED_NGC4945_ROWS:
        raise RuntimeError(
            f"generic audit does not contain exactly {EXPECTED_NGC4945_ROWS} "
            "rejected NGC 4945 rows"
        )

    # Importing the source-specific module is CASA-independent.  Reuse its
    # exact schema and mathematical validator instead of weakening that audit.
    from rerun_ngc4945_nuclear import (  # pylint: disable=import-outside-toplevel
        PHOTOMETRY_FIELDS,
        validate_photometry_records,
    )

    source18 = _read_exact_csv(
        source18_photometry_path,
        PHOTOMETRY_FIELDS,
        label="NGC 4945 source-18 photometry audit",
    )
    validate_photometry_records(source18)
    source18_products = [portable_release_product(row["product"]) for row in source18]
    if len(source18_products) != EXPECTED_NGC4945_ROWS or len(set(source18_products)) != EXPECTED_NGC4945_ROWS:
        raise RuntimeError(
            f"source-18 audit does not contain {EXPECTED_NGC4945_ROWS} unique products"
        )
    if set(source18_products) != generic_ngc4945:
        raise RuntimeError("source-18 and rejected generic source-17 products differ")
    table_ngc4945_order = [
        str(row[7]) for row in publication_rows if str(row[0]) == NGC4945
    ]
    if source18_products != table_ngc4945_order:
        raise RuntimeError("source-18 audit order/identity does not match the publication table")
    for record in source18:
        product = portable_release_product(record["product"])
        table_row = table_by_product[product]
        if (
            record["component_identity"] != SOURCE18_COMPONENT_IDENTITY
            or record["rejected_seed_identity"] != SOURCE17_COMPONENT_IDENTITY
            or record["component_identity"] == record["rejected_seed_identity"]
            or record["row_disposition"] != "retained_source18_replacement"
            or record["noise_guard_status"] != "passed"
        ):
            raise RuntimeError(
                f"source-18 retained component identity/disposition is invalid for "
                f"{record['product']}"
            )
        if (
            record["source"] != str(table_row[0])
            or record["telescope"] != str(table_row[1])
            or not math.isclose(
                _finite(record, "frequency_ghz", positive=True),
                float(table_row[2]),
                rel_tol=0.0,
                abs_tol=5e-3,
            )
            or record["date_label"] != str(table_row[4])
            or record["new_flux_field"] != str(table_row[3])
            or record["new_snr"] != str(table_row[6])
        ):
            raise RuntimeError(
                f"retained source-18 audit disagrees with publication row for {product}"
            )

    retained_generic = set(table_by_product) - generic_ngc4945
    retained_union = retained_generic | set(source18_products)
    if retained_union != set(table_by_product) or len(retained_union) != len(table_by_product):
        raise RuntimeError("each publication product must have exactly one retained diagnostic")
    return {
        "publication_rows": len(publication_rows),
        "admitted_input_additions": len(addition_rows),
        "publication_additions": len(publication_additions),
        "input_only_additions": len(addition_rows) - len(publication_additions),
        "generic_audit_rows": len(generic),
        "generic_publication_product_rows": len(generic) - boundary_rejected,
        "ordinary_generic_rows": ordinary,
        "ordinary_publication_product_rows": ordinary_publication,
        "forced_dirty_sma_rows": forced,
        "boundary_rejection_policy_rows": len(boundary_rows),
        "boundary_rejected_noise_audit_rows": boundary_rejected,
        "generic_ngc4945_rejected_rows": len(generic_ngc4945),
        "source18_retained_rows": len(source18),
        "retained_union_rows": len(retained_union),
    }
