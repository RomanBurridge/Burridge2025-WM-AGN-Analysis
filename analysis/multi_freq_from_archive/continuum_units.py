"""Strict unit conversion for continuum fit and image-noise quantities."""

from __future__ import annotations

import math
from typing import Any


def _normalized_unit(unit: Any) -> str:
    text = str(unit).strip().replace(" ", "").replace("µ", "u").casefold()
    if not text:
        raise ValueError("A continuum quantity has no reported unit")
    return text


def brightness_to_mjy_per_beam_factor(unit: Any) -> float:
    """Return the vectorizable multiplicative factor to mJy/beam."""

    factors = {
        "jy/beam": 1000.0,
        "jybeam-1": 1000.0,
        "mjy/beam": 1.0,
        "mjybeam-1": 1.0,
        "ujy/beam": 0.001,
        "ujybeam-1": 0.001,
    }
    normalized = _normalized_unit(unit)
    if normalized not in factors:
        raise ValueError(
            f"Unsupported image-brightness unit {unit!r}; expected a per-beam "
            "Jy, mJy, or uJy unit"
        )
    return factors[normalized]


def flux_density_to_mjy(value: Any, unit: Any) -> float:
    """Convert a CASA component flux density to mJy or fail explicitly."""

    factors = {"jy": 1000.0, "mjy": 1.0, "ujy": 0.001}
    normalized = _normalized_unit(unit)
    if normalized not in factors:
        raise ValueError(
            f"Unsupported component-flux unit {unit!r}; expected Jy, mJy, or uJy"
        )
    converted = float(value) * factors[normalized]
    if not math.isfinite(converted) and not math.isnan(converted):
        raise ValueError(f"Non-finite component flux: {value!r} {unit}")
    return converted


def brightness_to_mjy_per_beam(value: Any, unit: Any) -> float:
    """Convert an image value to mJy/beam or fail on incompatible units."""

    converted = float(value) * brightness_to_mjy_per_beam_factor(unit)
    if not math.isfinite(converted):
        raise ValueError(f"Non-finite image brightness: {value!r} {unit}")
    return converted
