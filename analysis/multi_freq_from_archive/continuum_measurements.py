"""Small, CASA-independent helpers for continuum measurements."""

from __future__ import annotations

import math
from datetime import date
from numbers import Integral
from typing import Any, NamedTuple


FINAL_NOISE_APERTURE_BMAJ_MULTIPLIER = 10.0
FINAL_NOISE_MIN_INDEPENDENT_BEAMS = 50.0
FINAL_NOISE_MEAN_SIGMA_REVIEW_THRESHOLD = 0.2
FINAL_NOISE_RMS_XMAD_REVIEW_THRESHOLD = 0.2
FIXED_BEAM_AXIS_REL_TOLERANCE = 1e-6
FIXED_BEAM_AXIS_ABS_TOLERANCE_ARCSEC = 1e-6
# CASA serializes a genuinely fixed ``abp`` position angle with observed
# round-trip differences up to 4.23e-4 deg in the reviewed image inventory.
# Keep PA's angular tolerance distinct from the arcsecond axis tolerance.
FIXED_BEAM_PA_ABS_TOLERANCE_DEG = 1e-3
GAUSSIAN_SIGMA_FROM_FWHM = 1.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))
THREE_SIGMA_FROM_FWHM = 3.0 * GAUSSIAN_SIGMA_FROM_FWHM
SCALED_MAD_NORMAL_FACTOR = 1.482602218505602


class OrdinaryFitClassification(NamedTuple):
    """Classification derived from one adopted source-excluded scatter."""

    is_detection: bool
    snr: float
    upper_limit_mjy: float


class FixedBeamShapeCheck(NamedTuple):
    """Numerical agreement between a fixed CASA fit and restoring beam."""

    major_delta_arcsec: float
    minor_delta_arcsec: float
    pa_delta_deg: float
    pa_guard_applicable: bool


class SquareApertureCoverage(NamedTuple):
    """Pixel-plane clearance of an axis-aligned square aperture."""

    half_width_x_pixels: float
    half_width_y_pixels: float
    left_clearance_pixels: float
    right_clearance_pixels: float
    bottom_clearance_pixels: float
    top_clearance_pixels: float
    minimum_edge_clearance_pixels: float


class ForcedFitDiagnostic(NamedTuple):
    """Optional component flux from the reviewed forced-limit fit."""

    converged: bool
    flux_value: float | None
    flux_unit: str | None


def forced_fit_diagnostic(imfit_result: Any) -> ForcedFitDiagnostic:
    """Extract a diagnostic flux only when a forced dirty-SMA fit converged.

    The forced publication limit is derived independently from the one-beam
    RMS.  A nonconverged Gaussian therefore has no scientifically meaningful
    fitted-flux diagnostic and must remain blank rather than inheriting an
    initial estimate, a zero, or an all-zero model value.
    """

    try:
        converged_values = imfit_result["converged"]
        if len(converged_values) != 1:
            raise ValueError("expected exactly one convergence flag")
        convergence_flag = converged_values[0]
        flag_type = type(convergence_flag)
        if flag_type is bool:
            converged = convergence_flag
        elif (
            flag_type.__name__ in {"bool", "bool_"}
            and flag_type.__module__.startswith("numpy")
        ):
            converged = bool(convergence_flag)
        else:
            raise ValueError("convergence flag is not boolean")
        results = imfit_result["results"]
        raw_nelements = results["nelements"]
        if isinstance(raw_nelements, bool) or not isinstance(
            raw_nelements, Integral
        ):
            raise ValueError("component count is not an integer")
        nelements = int(raw_nelements)
    except (KeyError, IndexError, TypeError, ValueError) as error:
        raise ValueError("Forced dirty-SMA imfit result lacks convergence status") from error
    if not converged:
        if nelements != 0 or "component0" in results or "deconvolved" in imfit_result:
            raise ValueError(
                "Nonconverged forced dirty-SMA imfit result contains a component"
            )
        return ForcedFitDiagnostic(False, None, None)

    if nelements != 1:
        raise ValueError(
            "Converged forced dirty-SMA imfit result is not single-component"
        )

    try:
        flux = results["component0"]["flux"]
        value = float(flux["value"][0])
        raw_unit = flux["unit"]
        if not isinstance(raw_unit, str):
            raise ValueError("component0 flux unit is not text")
        unit = raw_unit.strip()
    except (KeyError, IndexError, TypeError, ValueError) as error:
        raise ValueError(
            "Converged forced dirty-SMA imfit result lacks component0 flux"
        ) from error
    if not math.isfinite(value) or not unit:
        raise ValueError(
            "Converged forced dirty-SMA imfit component0 flux is invalid"
        )
    if "deconvolved" in imfit_result:
        deconvolved = imfit_result["deconvolved"]
        if deconvolved is None:
            raise ValueError("Forced dirty-SMA deconvolved record is null")
        try:
            deconvolved_flux = deconvolved["component0"]["flux"]
            deconvolved_value = float(deconvolved_flux["value"][0])
            raw_deconvolved_unit = deconvolved_flux["unit"]
            if not isinstance(raw_deconvolved_unit, str):
                raise ValueError("deconvolved flux unit is not text")
            deconvolved_unit = raw_deconvolved_unit.strip()
        except (KeyError, IndexError, TypeError, ValueError) as error:
            raise ValueError(
                "Forced dirty-SMA deconvolved component flux is malformed"
            ) from error
        if (
            not math.isfinite(deconvolved_value)
            or deconvolved_unit != unit
            or not math.isclose(
                deconvolved_value,
                value,
                rel_tol=1e-12,
                abs_tol=1e-15,
            )
        ):
            raise ValueError(
                "Forced dirty-SMA fitted and deconvolved flux diagnostics disagree"
            )
    return ForcedFitDiagnostic(True, value, unit)


def format_two_decimal(value: Any) -> str:
    """Return ordinary two-decimal formatting without deleting zero digits."""

    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Cannot format a non-finite continuum value: {value!r}")
    return f"{number:.2f}"


def background_rms_snr(flux_mjy: Any, background_rms_mjy_per_beam: Any) -> float:
    """Compute fit S/N from the source-excluded image RMS exactly once.

    The CASA formal component-fit error is residual-derived and is not known
    to be independent of this separately measured image RMS, so it must not be
    added to the RMS in quadrature.
    """

    flux = float(flux_mjy)
    rms = float(background_rms_mjy_per_beam)
    if not math.isfinite(flux):
        raise ValueError(f"Continuum flux must be finite: {flux_mjy!r}")
    if not math.isfinite(rms) or rms <= 0:
        raise ValueError(
            f"Background RMS must be finite and positive: {background_rms_mjy_per_beam!r}"
        )
    return flux / rms


def adopted_noise_snr(flux_mjy: Any, adopted_noise_mjy_per_beam: Any) -> float:
    """Compute S/N from the explicitly selected final noise estimator."""

    flux = float(flux_mjy)
    noise = float(adopted_noise_mjy_per_beam)
    if not math.isfinite(flux):
        raise ValueError(f"Continuum flux must be finite: {flux_mjy!r}")
    if not math.isfinite(noise) or noise <= 0.0:
        raise ValueError(
            f"Adopted noise must be finite and positive: {adopted_noise_mjy_per_beam!r}"
        )
    return flux / noise


def classify_ordinary_fit(
    flux_mjy: Any,
    adopted_noise_mjy_per_beam: Any,
) -> OrdinaryFitClassification:
    """Classify an ordinary fit from one adopted noise measurement.

    Final publication callers pass the local source-excluded scaled MAD;
    discovery callers retain their reviewed legacy RMS.  The returned upper
    limit is the larger of three times the adopted scatter and the fitted
    flux, so it cannot fall below the fitted signal at the threshold edge.
    Forced dirty-SMA limits remain a separate RMS policy.
    """

    flux = float(flux_mjy)
    noise = float(adopted_noise_mjy_per_beam)
    snr = adopted_noise_snr(flux, noise)
    upper_limit_mjy = max(3.0 * noise, flux)
    return OrdinaryFitClassification(
        is_detection=snr >= 3.0,
        snr=snr,
        upper_limit_mjy=upper_limit_mjy,
    )


def angular_value_to_arcsec(value: Any, unit: Any) -> float:
    """Convert one angular value to arcseconds with a closed unit set."""

    number = float(value)
    normalized = str(unit).strip().casefold()
    factors = {
        "rad": 180.0 / math.pi * 3600.0,
        "deg": 3600.0,
        "arcsec": 1.0,
    }
    if not math.isfinite(number):
        raise ValueError(f"Angular value must be finite: {value!r}")
    if normalized not in factors:
        raise ValueError(
            f"Unsupported angular unit {unit!r}; expected rad, deg, or arcsec"
        )
    return number * factors[normalized]


def normalize_right_ascension_degrees(value: Any) -> float:
    """Normalize a finite native longitude to the conventional [0, 360) range.

    CASA may serialize an otherwise valid celestial longitude as a negative
    angle (for example, -149.23 deg instead of 210.77 deg).  These are the
    same sky coordinate, so normalization belongs at the representation
    boundary; non-finite values still fail closed.
    """

    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Right ascension must be finite: {value!r}")
    return number % 360.0


def gaussian_beam_area_arcsec2(
    beam_major_arcsec: Any,
    beam_minor_arcsec: Any,
) -> float:
    """Return the Gaussian restoring-beam area, pi BMAJ BMIN/(4 ln 2)."""

    major = float(beam_major_arcsec)
    minor = float(beam_minor_arcsec)
    if not math.isfinite(major) or major <= 0.0:
        raise ValueError(f"Beam major axis must be finite and positive: {major!r}")
    if not math.isfinite(minor) or minor <= 0.0:
        raise ValueError(f"Beam minor axis must be finite and positive: {minor!r}")
    if minor > major:
        raise ValueError(
            f"Beam minor axis cannot exceed major axis: major={major}, minor={minor}"
        )
    return math.pi * major * minor / (4.0 * math.log(2.0))


def independent_beam_areas(
    valid_pixel_count: Any,
    pixel_width_arcsec: Any,
    pixel_height_arcsec: Any,
    beam_major_arcsec: Any,
    beam_minor_arcsec: Any,
) -> float:
    """Convert valid unmasked pixel area to Gaussian independent-beam areas."""

    count = float(valid_pixel_count)
    width = abs(float(pixel_width_arcsec))
    height = abs(float(pixel_height_arcsec))
    if not math.isfinite(count) or count < 0.0 or not count.is_integer():
        raise ValueError(f"Valid pixel count must be a finite integer: {valid_pixel_count!r}")
    if not math.isfinite(width) or width <= 0.0:
        raise ValueError(f"Pixel width must be finite and positive: {pixel_width_arcsec!r}")
    if not math.isfinite(height) or height <= 0.0:
        raise ValueError(f"Pixel height must be finite and positive: {pixel_height_arcsec!r}")
    return count * width * height / gaussian_beam_area_arcsec2(
        beam_major_arcsec,
        beam_minor_arcsec,
    )


def effective_independent_noise_samples(
    valid_pixel_count: Any,
    pixel_width_arcsec: Any,
    pixel_height_arcsec: Any,
    beam_major_arcsec: Any,
    beam_minor_arcsec: Any,
) -> float:
    """Bound beam-area samples by the number of actually measured pixels.

    The usual image-noise estimate is the valid area divided by the Gaussian
    beam area.  An undersampled image cannot contain more independent samples
    than valid pixels, however, so the fail-closed guard uses the smaller of
    those two quantities and records both in the audit.
    """

    count = float(valid_pixel_count)
    beam_areas = independent_beam_areas(
        count,
        pixel_width_arcsec,
        pixel_height_arcsec,
        beam_major_arcsec,
        beam_minor_arcsec,
    )
    return min(count, beam_areas)


def square_aperture_coverage(
    center_x_zero_based: Any,
    center_y_zero_based: Any,
    full_width_arcsec: Any,
    pixel_width_arcsec: Any,
    pixel_height_arcsec: Any,
    image_width_pixels: Any,
    image_height_pixels: Any,
) -> SquareApertureCoverage:
    """Measure how a requested square intersects the image footprint.

    CASA/FITS component coordinates are zero based.  Pixel cells span from
    ``-0.5`` through ``N - 0.5`` along each image axis, so those cell edges
    define whether the requested continuous 10-BMAJ square is geometrically
    complete.  Negative clearance records a requested aperture truncated by
    the finite image; callers retain the available valid intersection only
    when its independently measured area still passes the N_eff hard gate.
    This is deliberately separate from valid-pixel fraction because blanked
    pixels may reduce that fraction even for a fully contained aperture.
    """

    center_x = float(center_x_zero_based)
    center_y = float(center_y_zero_based)
    width = float(full_width_arcsec)
    pixel_width = abs(float(pixel_width_arcsec))
    pixel_height = abs(float(pixel_height_arcsec))
    image_width = float(image_width_pixels)
    image_height = float(image_height_pixels)
    finite_values = (
        center_x,
        center_y,
        width,
        pixel_width,
        pixel_height,
        image_width,
        image_height,
    )
    if any(not math.isfinite(value) for value in finite_values):
        raise ValueError(f"Aperture geometry must be finite: {finite_values!r}")
    if width <= 0.0 or pixel_width <= 0.0 or pixel_height <= 0.0:
        raise ValueError("Aperture width and pixel scales must be positive")
    if (
        image_width <= 0.0
        or image_height <= 0.0
        or not image_width.is_integer()
        or not image_height.is_integer()
    ):
        raise ValueError("Image dimensions must be positive integers")

    half_x = width / (2.0 * pixel_width)
    half_y = width / (2.0 * pixel_height)
    left = center_x - half_x + 0.5
    right = image_width - 0.5 - center_x - half_x
    bottom = center_y - half_y + 0.5
    top = image_height - 0.5 - center_y - half_y
    minimum = min(left, right, bottom, top)
    coverage = SquareApertureCoverage(
        half_width_x_pixels=half_x,
        half_width_y_pixels=half_y,
        left_clearance_pixels=left,
        right_clearance_pixels=right,
        bottom_clearance_pixels=bottom,
        top_clearance_pixels=top,
        minimum_edge_clearance_pixels=minimum,
    )
    return coverage


def require_minimum_independent_beams(
    measured_independent_beams: Any,
    minimum_independent_beams: Any = FINAL_NOISE_MIN_INDEPENDENT_BEAMS,
) -> float:
    """Fail closed when a final noise aperture has too little usable area."""

    measured = float(measured_independent_beams)
    minimum = float(minimum_independent_beams)
    if not math.isfinite(measured) or measured < 0.0:
        raise ValueError(
            f"Independent-beam area must be finite and nonnegative: {measured!r}"
        )
    if not math.isfinite(minimum) or minimum <= 0.0:
        raise ValueError(
            f"Minimum independent-beam area must be finite and positive: {minimum!r}"
        )
    if measured < minimum:
        raise ValueError(
            f"Final source-excluded aperture has {measured:.6f} independent beams; "
            f"at least {minimum:.6f} are required"
        )
    return measured


def position_angle_delta_deg(first_deg: Any, second_deg: Any) -> float:
    """Return the smallest difference between ellipse PAs (period 180 deg)."""

    first = float(first_deg)
    second = float(second_deg)
    if not math.isfinite(first) or not math.isfinite(second):
        raise ValueError(f"Position angles must be finite: {first_deg!r}, {second_deg!r}")
    return abs((first - second + 90.0) % 180.0 - 90.0)


def require_fixed_restoring_beam_shape(
    fitted_major_arcsec: Any,
    fitted_minor_arcsec: Any,
    fitted_pa_deg: Any,
    beam_major_arcsec: Any,
    beam_minor_arcsec: Any,
    beam_pa_deg: Any,
    *,
    axis_relative_tolerance: Any = FIXED_BEAM_AXIS_REL_TOLERANCE,
    axis_absolute_tolerance_arcsec: Any = FIXED_BEAM_AXIS_ABS_TOLERANCE_ARCSEC,
    pa_absolute_tolerance_deg: Any = FIXED_BEAM_PA_ABS_TOLERANCE_DEG,
) -> FixedBeamShapeCheck:
    """Require the convolved final fit to retain the restoring-beam shape.

    Only this condition makes integrated component mJy directly comparable to
    a local image RMS in mJy/beam.  The tolerances are deliberately frozen and
    tight enough to catch a lost CASA ``abp`` constraint while allowing normal
    floating-point roundoff in a genuinely fixed fit.
    """

    fitted_major = float(fitted_major_arcsec)
    fitted_minor = float(fitted_minor_arcsec)
    beam_major = float(beam_major_arcsec)
    beam_minor = float(beam_minor_arcsec)
    relative = float(axis_relative_tolerance)
    absolute = float(axis_absolute_tolerance_arcsec)
    pa_tolerance = float(pa_absolute_tolerance_deg)
    values = (fitted_major, fitted_minor, beam_major, beam_minor)
    if any(not math.isfinite(value) or value <= 0.0 for value in values):
        raise ValueError(f"Fit and beam axes must be finite and positive: {values!r}")
    if any(
        not math.isfinite(value) or value < 0.0
        for value in (relative, absolute, pa_tolerance)
    ):
        raise ValueError("Fixed-beam tolerances must be finite and nonnegative")

    pa_delta = position_angle_delta_deg(fitted_pa_deg, beam_pa_deg)
    # Position angle is undefined for a circular restoring beam.  The same
    # frozen axis tolerances used by the shape guard decide whether the beam is
    # numerically circular; PA is still recorded but is not a hard condition
    # in that degenerate case.
    pa_guard_applicable = not math.isclose(
        beam_major,
        beam_minor,
        rel_tol=relative,
        abs_tol=absolute,
    )
    check = FixedBeamShapeCheck(
        major_delta_arcsec=abs(fitted_major - beam_major),
        minor_delta_arcsec=abs(fitted_minor - beam_minor),
        pa_delta_deg=pa_delta,
        pa_guard_applicable=pa_guard_applicable,
    )
    major_matches = math.isclose(
        fitted_major,
        beam_major,
        rel_tol=relative,
        abs_tol=absolute,
    )
    minor_matches = math.isclose(
        fitted_minor,
        beam_minor,
        rel_tol=relative,
        abs_tol=absolute,
    )
    if (
        not major_matches
        or not minor_matches
        or (pa_guard_applicable and check.pa_delta_deg > pa_tolerance)
    ):
        raise ValueError(
            "Final convolved component is not fixed to the restoring beam: "
            f"major={fitted_major} vs {beam_major} arcsec, "
            f"minor={fitted_minor} vs {beam_minor} arcsec, "
            f"PA delta={check.pa_delta_deg} deg"
        )
    return check


def final_noise_review_flags(
    mean_mjy_per_beam: Any,
    sigma_mjy_per_beam: Any,
    rms_mjy_per_beam: Any,
    scaled_mad_mjy_per_beam: Any,
    *,
    mean_sigma_threshold: Any = FINAL_NOISE_MEAN_SIGMA_REVIEW_THRESHOLD,
    rms_xmad_threshold: Any = FINAL_NOISE_RMS_XMAD_REVIEW_THRESHOLD,
    aperture_minimum_edge_clearance_pixels: Any | None = None,
) -> tuple[str, ...]:
    """Return deterministic diagnostics alongside the adopted scaled MAD."""

    mean = float(mean_mjy_per_beam)
    sigma = float(sigma_mjy_per_beam)
    rms = float(rms_mjy_per_beam)
    scaled_mad = float(scaled_mad_mjy_per_beam)
    mean_threshold = float(mean_sigma_threshold)
    mad_threshold = float(rms_xmad_threshold)
    if any(not math.isfinite(value) for value in (mean, sigma, rms, scaled_mad)):
        raise ValueError("Final noise statistics must all be finite")
    if sigma <= 0.0 or rms <= 0.0 or scaled_mad <= 0.0:
        raise ValueError(
            "Final sigma, RMS, and scaled MAD must be positive"
        )
    if mean_threshold < 0.0 or mad_threshold < 0.0:
        raise ValueError("Final-noise review thresholds must be nonnegative")
    edge_clearance = None
    if aperture_minimum_edge_clearance_pixels is not None:
        edge_clearance = float(aperture_minimum_edge_clearance_pixels)
        if not math.isfinite(edge_clearance):
            raise ValueError("Aperture edge clearance must be finite")

    flags = []
    if abs(mean) / sigma > mean_threshold:
        flags.append("nonzero_background_mean")
    if abs(rms / scaled_mad - 1.0) > mad_threshold:
        flags.append("rms_scaled_mad_disagreement")
    # A 1e-9-pixel allowance treats only floating-point boundary equality as
    # fully contained.  Truncation is reviewable, not an automatic RMS
    # substitution or failure; N_eff remains the hard available-area guard.
    if edge_clearance is not None and edge_clearance < -1e-9:
        flags.append("aperture_truncated_by_image")
    return tuple(flags)


def square_bmaj_region_width_arcsec(
    beam_major_arcsec: Any,
    multiplier: Any,
    *,
    minimum_arcsec: Any | None = None,
    maximum_arcsec: Any | None = None,
    override_arcsec: Any | None = None,
) -> float:
    """Return the full side width of the intentional square search region.

    The continuum search is isotropic on the sky and conservatively uses the
    restoring-beam major-axis FWHM for both square sides.  It is not a rotated
    beam ellipse; the minor axis therefore has no role in this construction.
    CASA ``centerbox`` and the square ``rotbox`` fallback both interpret this
    value as the full side width.
    """

    major = float(beam_major_arcsec)
    scale = float(multiplier)
    if not math.isfinite(major) or major <= 0:
        raise ValueError(f"Beam major axis must be finite and positive: {beam_major_arcsec!r}")
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError(f"Region multiplier must be finite and positive: {multiplier!r}")
    width = major * scale
    if minimum_arcsec is not None:
        minimum = float(minimum_arcsec)
        if not math.isfinite(minimum) or minimum <= 0:
            raise ValueError(f"Minimum region width must be positive: {minimum_arcsec!r}")
        width = max(width, minimum)
    if maximum_arcsec is not None:
        maximum = float(maximum_arcsec)
        if not math.isfinite(maximum) or maximum <= 0:
            raise ValueError(f"Maximum region width must be positive: {maximum_arcsec!r}")
        width = min(width, maximum)
    if override_arcsec is not None:
        override = float(override_arcsec)
        if not math.isfinite(override) or override <= 0:
            raise ValueError(f"Override region width must be positive: {override_arcsec!r}")
        width = override
    return width


def calendar_date_sort_key(value: Any) -> tuple[int, int, int]:
    """Return a strict calendar-date key while ignoring the ``_N`` row suffix.

    The suffix distinguishes multiple products from one calendar date; it is
    not part of the date itself.  Missing dates sort before real dates, while
    malformed dates fail closed instead of silently changing row order.
    """

    raw = str(value).strip().strip(",")
    if not raw or raw.upper() == "NA":
        return (0, 0, 0)

    calendar_text, separator, suffix = raw.partition("_")
    if separator and (not suffix.isdigit() or int(suffix) < 1):
        raise ValueError(f"Invalid continuum date suffix: {value!r}")
    pieces = calendar_text.split("/")
    if len(pieces) != 3:
        raise ValueError(f"Invalid continuum calendar date: {value!r}")
    try:
        parsed = date(*(int(piece) for piece in pieces))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid continuum calendar date: {value!r}") from exc
    return (parsed.year, parsed.month, parsed.day)


def portable_release_product(value: Any) -> str:
    """Return a release-relative continuum product identifier.

    Historical SMA summary rows recorded the author's absolute checkout path.
    Preserve the repository-relative scientific provenance below ``analysis``
    while rejecting unknown absolute paths instead of leaking them silently.
    ALMA/VLA rows already stored relative to ``multi_freq_from_archive`` remain
    unchanged.
    """

    product = str(value).replace("\\n", "").replace("\n", "").strip()
    product = product.replace("\\", "/")
    for marker, prefix in (
        ("/analysis/sma_calibration/", "analysis/sma_calibration/"),
        ("/analysis/multi_freq_from_archive/", ""),
    ):
        if marker in product:
            return prefix + product.split(marker, 1)[1]
    windows_absolute = len(product) >= 3 and product[1:3] == ":/"
    if product.startswith("/") or windows_absolute:
        raise ValueError(
            "Cannot make an unknown absolute continuum product path portable: "
            f"{product!r}"
        )
    return product
