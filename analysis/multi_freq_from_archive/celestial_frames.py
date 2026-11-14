"""Pure helpers for routing coordinates through an image's celestial frame.

CASA reports FITS ``RADESYS=FK5`` images as ``J2000`` direction coordinate
systems.  The continuum fitter keeps coherent ICRS and J2000 coordinate pairs
and uses these helpers to select the pair matching the image metadata rather
than guessing from the telescope name.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any


ICRS = "ICRS"
J2000 = "J2000"


def colon_sexagesimal_declination(raw_declination: Any) -> str:
    """Normalize a sexagesimal declination without losing fractional seconds.

    CASA reports declinations with dot-separated degree and minute fields,
    while analysisUtils returns colon-separated fields and callers may also
    supply CASA ``d/m/s`` notation.  Only the first two separators belong to
    the sexagesimal fields; a decimal point in the seconds field must remain a
    decimal point.
    """

    value = str(raw_declination).strip()
    if not value:
        raise ValueError("Declination must be non-empty")

    if "d" in value:
        if value.count("d") != 1 or value.count("m") != 1:
            raise ValueError(f"Cannot format declination {raw_declination!r}")
        value = value.rstrip("s")
        value = value.replace("d", ":", 1).replace("m", ":", 1)

    if ":" in value:
        parts = value.split(":")
    else:
        # ``maxsplit=2`` is essential: it preserves the fractional-seconds
        # decimal point in values such as ``-65.20.20.970``.
        parts = value.split(".", 2)

    if len(parts) != 3 or not all(part.strip() for part in parts):
        raise ValueError(f"Cannot format declination {raw_declination!r}")
    try:
        float(parts[0])
        float(parts[1])
        float(parts[2])
    except ValueError as error:
        raise ValueError(f"Cannot format declination {raw_declination!r}") from error
    return ":".join(part.strip() for part in parts)


def casa_sexagesimal_declination(raw_declination: Any) -> str:
    """Return a CASA ``d/m/s`` declination with fractional seconds intact."""

    degrees, minutes, seconds = colon_sexagesimal_declination(raw_declination).split(
        ":"
    )
    return f"{degrees}d{minutes}m{seconds}s"


def coordinate_candidate_attempt_limit(getfiles: Any, *, choosefile: bool) -> int:
    """Return a finite upper bound for coordinate-image selection attempts."""

    if choosefile:
        return 1
    return sum(
        1
        for telescope_group in list(getfiles)[1:]
        for telescope in telescope_group
        for candidate in telescope
        if candidate is not None and len(candidate) > 0
    )


def discard_failed_coordinate_candidate_or_raise(
    getfiles: Any,
    best_file: Any,
    *,
    choosefile: bool,
    error: Exception,
) -> None:
    """Ensure a coordinate-selection failure cannot retry the same file forever.

    Automatic selection discards the failed candidate before the next bounded
    attempt.  A specifically chosen image cannot fall back to another file, so
    its failure is raised immediately.
    """

    if choosefile:
        raise RuntimeError(
            f"Chosen coordinate image failed and cannot be retried: {best_file}"
        ) from error

    for telescope_group in list(getfiles)[1:]:
        for telescope in telescope_group:
            for index, candidate in enumerate(telescope):
                if candidate == best_file:
                    telescope[index] = None
                    return
    raise RuntimeError(
        f"Failed coordinate image was not present in the candidate set: {best_file}"
    ) from error


def normalize_celestial_reference_code(reference_code: Any) -> str:
    """Normalize CASA/FITS equatorial reference codes to ICRS or J2000.

    ``reference_code`` may be the string returned by a FITS header or the
    one-element sequence returned by
    ``coordsys.referencecode(type='direction')``.  FK5 is represented by
    CASA as J2000 for these images, so both spellings route to the same
    coordinate pair.  Unknown, missing, or contradictory metadata fails
    closed instead of silently selecting a frame.
    """

    if isinstance(reference_code, str):
        raw_values = [reference_code]
    elif isinstance(reference_code, Iterable):
        raw_values = list(reference_code)
    else:
        raw_values = [reference_code]

    normalized: set[str] = set()
    for raw in raw_values:
        value = str(raw).strip().upper().replace("_", "")
        if not value:
            continue
        if value == ICRS:
            normalized.add(ICRS)
        elif value in {J2000, "FK5", "FK5J2000", "EQUATORIALJ2000"}:
            normalized.add(J2000)
        else:
            raise ValueError(f"Unsupported celestial reference code: {raw!r}")

    if len(normalized) != 1:
        raise ValueError(
            "Celestial reference metadata must identify exactly one of "
            f"ICRS or J2000/FK5; received {reference_code!r}"
        )
    return normalized.pop()


def coordinate_pair_for_reference_code(
    reference_code: Any,
    *,
    icrs: Sequence[str],
    j2000: Sequence[str],
) -> tuple[str, str]:
    """Return the coherent coordinate pair matching image frame metadata."""

    if len(icrs) != 2 or len(j2000) != 2:
        raise ValueError("ICRS and J2000 coordinates must each contain RA and Dec")
    icrs_pair = (str(icrs[0]).strip(), str(icrs[1]).strip())
    j2000_pair = (str(j2000[0]).strip(), str(j2000[1]).strip())
    if not all(icrs_pair) or not all(j2000_pair):
        raise ValueError("ICRS and J2000 coordinate values must be non-empty")

    frame = normalize_celestial_reference_code(reference_code)
    return icrs_pair if frame == ICRS else j2000_pair


def casa_direction_for_reference_code(
    reference_code: Any,
    *,
    icrs: Sequence[str],
    j2000: Sequence[str],
) -> str:
    """Build a CASA direction string from coherent, metadata-routed values."""

    frame = normalize_celestial_reference_code(reference_code)
    ra, dec = coordinate_pair_for_reference_code(
        frame,
        icrs=icrs,
        j2000=j2000,
    )
    return f"{frame} {ra} {dec}"
