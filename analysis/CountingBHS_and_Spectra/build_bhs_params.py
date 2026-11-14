#!/usr/bin/env python3
"""Build the audited 21-source black-hole shadow parameter table.

The source CSV records published mass normalizations, adopted distances, and
the evidence-based uncertainty treatment for every source.  Paired-chain
angular gravitational-radius quantiles live in a separate author-supplied
input table.  External distance factors cancel from ``M/D`` for every
``shared_distance`` row; they are retained only in the distance and rescaled
mass records.

Non-probabilistic sensitivity ranges and conditional literature fits that are
not adopted are validated from separate sidecars and attached only after the
primary mass, shadow-size, and uncertainty fields have been constructed.
Checks that cannot be attached because a source has no robust primary are
validated from a distinct audit-only sidecar and never enter generated values.

All calculations use the IAU 2015 nominal solar mass parameter, the exact
speed of light, the exact astronomical unit, and the parsec definition.
"""

from __future__ import annotations

import argparse
import csv
import io
import math
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "analysis" / "CountingBHS_and_Spectra"
DEFAULT_INPUT = DATA_DIR / "bhs_source_inputs.csv"
DEFAULT_QUANTILES = DATA_DIR / "pesce_theta_g_quantiles.csv"
DEFAULT_MODEL_CHECKS = DATA_DIR / "bhs_model_sensitivity_checks.csv"
DEFAULT_AUDIT_ONLY_MODEL_CHECKS = (
    DATA_DIR / "bhs_audit_only_model_sensitivity_checks.csv"
)
DEFAULT_NONADOPTED_FITS = DATA_DIR / "bhs_nonadopted_ratio_fits.csv"
DEFAULT_OUTPUT = DATA_DIR / "BHS_params.csv"

# IAU 2015 Resolution B3 nominal solar mass parameter (m^3 s^-2).
NOMINAL_SOLAR_MASS_PARAMETER = 1.3271244e20
# SI definition, exact (m s^-1).
SPEED_OF_LIGHT = 299_792_458.0
# IAU astronomical-unit definition, exact (m).
ASTRONOMICAL_UNIT_M = 149_597_870_700.0
# Parsecs are defined by one astronomical unit subtending one arcsecond.
PARSEC_M = ASTRONOMICAL_UNIT_M * 648_000.0 / math.pi
MPC_M = 1_000_000.0 * PARSEC_M
RAD_TO_UAS = 180.0 / math.pi * 3_600.0 * 1_000_000.0
SHADOW_TO_GRAVITATIONAL_RADIUS = 6.0 * math.sqrt(3.0)

TREATMENT_COUNTS = {
    "paired_posterior": 6,
    "shared_distance": 13,
    "no_robust_mass": 2,
}

INPUT_FIELDS = [
    "source_order",
    "source_key",
    "source_label",
    "treatment",
    "reference_mass_msun",
    "reference_mass_unc_minus_msun",
    "reference_mass_unc_plus_msun",
    "reference_distance_mpc",
    "adopted_distance_mpc",
    "distance_measurement_unc_minus_mpc",
    "distance_measurement_unc_plus_mpc",
    "distance_flow_unc_mpc",
    "distance_shared_systematic_unc_mpc",
    "size_stat_fraction_minus_override",
    "size_stat_fraction_plus_override",
    "size_systematic_fraction_minus",
    "size_systematic_fraction_plus",
    "systematic_status",
    "mass_reference_doi",
    "mass_uncertainty_reference_doi",
    "distance_reference_doi",
    "uncertainty_note",
    "mass_display_interval_status",
    "distance_stat_component_status",
    "reference_mass_systematic_unc_minus_msun",
    "reference_mass_systematic_unc_plus_msun",
    "mass_systematic_status",
]

# The first twelve columns are retained for compatibility with existing
# continuum/contamination consumers.  The explicit component columns are the
# canonical uncertainty record.
OUTPUT_FIELDS = [
    "source_key",
    "source_label",
    "size_uas",
    "size_unc_independent_uas",
    "size_unc_correlated_uas",
    "size_unc_adopted_uas",
    "adopted_uncertainty_model",
    "log_mass",
    "log_mass_unc",
    "distance_mpc",
    "distance_unc_mpc",
    "mass_frac_percent",
    "mass_msun",
    "mass_unc_minus_msun",
    "mass_unc_plus_msun",
    "log_mass_unc_minus",
    "log_mass_unc_plus",
    "mass_systematic_unc_minus_msun",
    "mass_systematic_unc_plus_msun",
    "log_mass_systematic_unc_minus",
    "log_mass_systematic_unc_plus",
    "mass_systematic_status",
    "distance_unc_minus_mpc",
    "distance_unc_plus_mpc",
    "parameter_reference_doi",
    "mass_uncertainty_reference_doi",
    "uncertainty_status",
    "size_posterior_q16_uas",
    "size_posterior_q50_uas",
    "size_posterior_q84_uas",
    "size_unc_minus_uas",
    "size_unc_plus_uas",
    "treatment",
    "mass_reference_distance_mpc",
    "mass_fit_unc_minus_msun",
    "mass_fit_unc_plus_msun",
    "mass_shared_distance_unc_minus_msun",
    "mass_shared_distance_unc_plus_msun",
    "mass_total_unc_minus_msun",
    "mass_total_unc_plus_msun",
    "mass_display_unc_minus_msun",
    "mass_display_unc_plus_msun",
    "mass_display_interval_status",
    "distance_stat_unc_minus_mpc",
    "distance_stat_unc_plus_mpc",
    "distance_stat_component_status",
    "distance_systematic_unc_minus_mpc",
    "distance_systematic_unc_plus_mpc",
    "size_stat_unc_minus_uas",
    "size_stat_unc_plus_uas",
    "size_stat_component_status",
    "size_systematic_unc_minus_uas",
    "size_systematic_unc_plus_uas",
    "size_total_unc_minus_uas",
    "size_total_unc_plus_uas",
    "size_systematic_status",
    "size_total_status",
    "legacy_uncertainty_basis",
    "distance_reference_doi",
    "ratio_provenance",
    "source_order",
    "size_model_check_key",
    "size_model_check_type",
    "size_model_check_data_relation",
    "size_model_check_range_basis",
    "size_model_check_accepted_configuration_count",
    "size_model_check_excluded_configuration_ids",
    "size_model_check_lower_uas",
    "size_model_check_upper_uas",
    "size_model_check_probability_status",
    "size_model_check_adoption_status",
    "size_model_check_reference_key",
    "size_model_check_reference_doi",
    "nonadopted_fit_key",
    "nonadopted_fit_ratio_msun_per_mpc",
    "nonadopted_fit_ratio_unc_minus_msun_per_mpc",
    "nonadopted_fit_ratio_unc_plus_msun_per_mpc",
    "nonadopted_fit_size_uas",
    "nonadopted_fit_size_unc_minus_uas",
    "nonadopted_fit_size_unc_plus_uas",
    "nonadopted_fit_interval_status",
    "nonadopted_fit_adoption_status",
    "nonadopted_fit_basis",
    "nonadopted_fit_reference_key",
    "nonadopted_fit_reference_doi",
]

MASS_DISPLAY_INTERVAL_STATUSES = {
    "decomposed_components",
    "separate_decomposed_and_geometry_systematic_components",
    "separate_decomposed_and_position_systematic_components",
    "separate_statistical_and_systematic_components",
    "source_reported_total",
    "not_available",
}
MASS_SYSTEMATIC_STATUSES = {
    "source_reported_absolute_position_systematic",
    "source_reported_midline_inclination_systematic",
    "source_reported_model_selection_systematic",
    "not_available",
}
DISTANCE_STAT_COMPONENT_STATUSES = {
    "statistical",
    "reported_legacy_multimethod",
    "reported_legacy_tully_fisher",
    "not_available",
}
FIT_PLUS_ABSOLUTE_POSITION_SOURCES = {"NGC4388", "ESO558-G009"}
NGC4388_UNCERTAINTY_NOTE = (
    "kuo2011_table4_mass_8p4e6_preferred_over_section3p3p5_table6_mass_"
    "8p5e6_discrepancy_"
    "fit_plus_absolute_position_component_tully_fisher_cepheid_calibration_"
    "systematic_rms_five_spot_no_systemic_sparse_disk_geometry_unquantified"
)
SIZE_STAT_COMPONENT_STATUSES = {
    "paired_chain_quantile_interval",
    "reported_undecomposed_uncertainty",
    "fit_plus_absolute_position_mixed_component",
    "fixed_distance_formal_fit_interval",
    "angular_radius_only_component",
    "no_formal_uncertainty",
    "not_available",
}


def optional_float(value: str | None) -> float | None:
    """Parse one optional numeric CSV field."""

    if value is None or value.strip() in {"", "NA"}:
        return None
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"non-finite numeric input: {value!r}")
    return parsed


def mean_optional(lower: float | None, upper: float | None) -> float | None:
    if lower is None and upper is None:
        return None
    if lower is None or upper is None:
        raise ValueError("asymmetric uncertainty requires both bounds")
    return (lower + upper) / 2.0


def positive_quadrature(*values: float | None) -> float | None:
    retained = [value for value in values if value is not None and value > 0.0]
    if not retained:
        return None
    return math.sqrt(sum(value * value for value in retained))


def shadow_size_uas(mass_msun: float, distance_mpc: float) -> float:
    """Return the diameter of a Schwarzschild photon capture shadow in uas."""

    if mass_msun <= 0.0 or distance_mpc <= 0.0:
        raise ValueError("mass and distance must be positive")
    theta_radians = (
        SHADOW_TO_GRAVITATIONAL_RADIUS
        * NOMINAL_SOLAR_MASS_PARAMETER
        * mass_msun
        / (distance_mpc * MPC_M * SPEED_OF_LIGHT**2)
    )
    return theta_radians * RAD_TO_UAS


def theta_g_nas_to_shadow_uas(theta_g_nas: float) -> float:
    if theta_g_nas <= 0.0:
        raise ValueError("angular gravitational radius must be positive")
    return theta_g_nas * SHADOW_TO_GRAVITATIONAL_RADIUS / 1_000.0


def _validate_header(actual: list[str] | None, expected: list[str], path: Path) -> None:
    if actual != expected:
        raise ValueError(
            f"unexpected columns in {path}: expected {expected!r}, found {actual!r}"
        )


def load_source_inputs(path: Path = DEFAULT_INPUT) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        _validate_header(reader.fieldnames, INPUT_FIELDS, path)
        rows = list(reader)

    if any(None in row for row in rows):
        raise ValueError(f"one or more rows in {path} contain extra columns")
    if len(rows) != 21:
        raise ValueError(f"expected 21 source inputs, found {len(rows)}")

    orders = [int(row["source_order"]) for row in rows]
    if orders != list(range(1, 22)):
        raise ValueError("source_order must be the consecutive sequence 1..21")

    keys = [row["source_key"] for row in rows]
    if len(set(keys)) != len(keys):
        raise ValueError("source_key values must be unique")

    counts = Counter(row["treatment"] for row in rows)
    if dict(counts) != TREATMENT_COUNTS:
        raise ValueError(
            f"treatment counts must be {TREATMENT_COUNTS!r}, found {dict(counts)!r}"
        )
    if "independent" in counts:
        raise ValueError("the audited sample contains no independent-error cases")

    valid_statuses = {
        "paired_quantile_file",
        "none_adopted",
        "reported_undecomposed_uncertainty",
        "numeric_systematic",
        "model_sensitivity_envelope",
        "unquantified_model",
        "unquantified_astrometric",
        "pending_nakai",
        "no_formal_uncertainty",
        "not_applicable",
    }
    for row in rows:
        key = row["source_key"]
        distance = optional_float(row["adopted_distance_mpc"])
        if distance is None or distance <= 0.0:
            raise ValueError(f"{key}: adopted distance must be positive")
        if row["systematic_status"] not in valid_statuses:
            raise ValueError(f"{key}: invalid systematic_status")
        mass_display_status = row["mass_display_interval_status"]
        if mass_display_status not in MASS_DISPLAY_INTERVAL_STATUSES:
            raise ValueError(f"{key}: invalid mass_display_interval_status")
        mass_systematic_status = row["mass_systematic_status"]
        if mass_systematic_status not in MASS_SYSTEMATIC_STATUSES:
            raise ValueError(f"{key}: invalid mass_systematic_status")
        mass_systematic_minus = optional_float(
            row["reference_mass_systematic_unc_minus_msun"]
        )
        mass_systematic_plus = optional_float(
            row["reference_mass_systematic_unc_plus_msun"]
        )
        if (mass_systematic_minus is None) != (mass_systematic_plus is None):
            raise ValueError(f"{key}: asymmetric mass systematic requires both bounds")
        if mass_systematic_minus is not None and (
            mass_systematic_minus <= 0.0 or mass_systematic_plus <= 0.0
        ):
            raise ValueError(f"{key}: mass systematic bounds must be positive")
        if mass_systematic_status == "not_available":
            if mass_systematic_minus is not None:
                raise ValueError(
                    f"{key}: unavailable mass systematic cannot carry numerical bounds"
                )
        elif mass_systematic_minus is None:
            raise ValueError(
                f"{key}: reported mass systematic requires numerical bounds"
            )
        distance_component_status = row["distance_stat_component_status"]
        if distance_component_status not in DISTANCE_STAT_COMPONENT_STATUSES:
            raise ValueError(f"{key}: invalid distance_stat_component_status")
        distance_record = distance_components(row)
        has_distance_stat = (
            distance_record["stat_minus"] is not None
            and distance_record["stat_plus"] is not None
        )
        if distance_component_status == "not_available":
            if has_distance_stat:
                raise ValueError(
                    f"{key}: not-available distance status cannot carry a numerical component"
                )
        elif not has_distance_stat:
            raise ValueError(
                f"{key}: distance component status requires numerical bounds"
            )
        if not row["mass_reference_doi"] or not row["distance_reference_doi"]:
            raise ValueError(f"{key}: mass and distance references are required")

        mass_uncertainty_reference = row["mass_uncertainty_reference_doi"]
        if not mass_uncertainty_reference:
            raise ValueError(
                f"{key}: mass uncertainty reference must be a DOI or literal NA"
            )
        if (
            mass_uncertainty_reference != "NA"
            and not mass_uncertainty_reference.startswith("10.")
        ):
            raise ValueError(
                f"{key}: mass uncertainty reference must be a DOI or literal NA"
            )

        mass = optional_float(row["reference_mass_msun"])
        reference_distance = optional_float(row["reference_distance_mpc"])
        if row["treatment"] == "no_robust_mass":
            if mass is not None or reference_distance is not None:
                raise ValueError(f"{key}: no-robust-mass row must not supply a mass")
            if mass_systematic_status != "not_available":
                raise ValueError(
                    f"{key}: no-robust-mass row must not supply a mass systematic"
                )
            if mass_uncertainty_reference != "NA":
                raise ValueError(
                    f"{key}: no-robust-mass row must use NA uncertainty provenance"
                )
            if mass_display_status != "not_available":
                raise ValueError(
                    f"{key}: no-robust-mass row must not display a mass interval"
                )
            continue
        if mass is None or mass <= 0.0 or reference_distance is None or reference_distance <= 0.0:
            raise ValueError(f"{key}: positive reference mass and distance are required")
        if row["treatment"] == "paired_posterior" and not math.isclose(
            distance, reference_distance, rel_tol=0.0, abs_tol=1e-12
        ):
            raise ValueError(f"{key}: paired marginal distance must equal adopted distance")

        input_uncertainty_minus = optional_float(row["reference_mass_unc_minus_msun"])
        input_uncertainty_plus = optional_float(row["reference_mass_unc_plus_msun"])
        override_minus = optional_float(row["size_stat_fraction_minus_override"])
        override_plus = optional_float(row["size_stat_fraction_plus_override"])
        has_mass_uncertainty = (
            input_uncertainty_minus is not None
            or input_uncertainty_plus is not None
            or override_minus is not None
            or override_plus is not None
            or mass_systematic_minus is not None
            or mass_systematic_plus is not None
        )
        if has_mass_uncertainty and mass_uncertainty_reference == "NA":
            raise ValueError(f"{key}: numerical mass uncertainty requires a DOI")
        if not has_mass_uncertainty and mass_uncertainty_reference != "NA":
            raise ValueError(
                f"{key}: absent mass uncertainty must use literal NA provenance"
            )
        if mass_display_status == "not_available":
            if has_mass_uncertainty:
                raise ValueError(
                    f"{key}: not-available mass display cannot carry uncertainty inputs"
                )
        elif not has_mass_uncertainty:
            raise ValueError(
                f"{key}: mass display status requires numerical uncertainty inputs"
            )
        if mass_display_status == "source_reported_total" and (
            input_uncertainty_minus is None or input_uncertainty_plus is None
        ):
            raise ValueError(
                f"{key}: source-reported mass total requires both published bounds"
            )

        for suffix in ("minus", "plus"):
            override = optional_float(row[f"size_stat_fraction_{suffix}_override"])
            systematic = optional_float(row[f"size_systematic_fraction_{suffix}"])
            if override is not None and override < 0.0:
                raise ValueError(f"{key}: statistical fraction cannot be negative")
            if systematic is not None and systematic < 0.0:
                raise ValueError(f"{key}: systematic fraction cannot be negative")

    reported_sources = {
        row["source_key"]
        for row in rows
        if row["systematic_status"] == "reported_undecomposed_uncertainty"
    }
    if reported_sources != {"NGC1194", "NGC2273"}:
        raise ValueError(
            "reported undecomposed uncertainties must apply exactly to "
            "NGC1194 and NGC2273"
        )

    unquantified_astrometric_sources = {
        row["source_key"]
        for row in rows
        if row["systematic_status"] == "unquantified_astrometric"
    }
    if unquantified_astrometric_sources != {"NGC2960"}:
        raise ValueError(
            "unquantified astrometric completeness must apply exactly to NGC2960"
        )

    source_reported_mass_totals = {
        row["source_key"]
        for row in rows
        if row["mass_display_interval_status"] == "source_reported_total"
    }
    if source_reported_mass_totals != {"Circinus"}:
        raise ValueError(
            "source-reported mass display totals must apply exactly to Circinus"
        )
    separate_mass_components = {
        row["source_key"]
        for row in rows
        if row["mass_display_interval_status"]
        == "separate_statistical_and_systematic_components"
    }
    if separate_mass_components != {"CGCG074-064"}:
        raise ValueError(
            "separate statistical and systematic mass components must apply "
            "exactly to CGCG074-064"
        )
    separate_position_components = {
        row["source_key"]
        for row in rows
        if row["mass_display_interval_status"]
        == "separate_decomposed_and_position_systematic_components"
    }
    if separate_position_components != {"UGC6093"}:
        raise ValueError(
            "separate decomposed and position-systematic mass components must "
            "apply exactly to UGC6093"
        )
    separate_geometry_components = {
        row["source_key"]
        for row in rows
        if row["mass_display_interval_status"]
        == "separate_decomposed_and_geometry_systematic_components"
    }
    if separate_geometry_components != {"ESO558-G009"}:
        raise ValueError(
            "separate decomposed and geometry-systematic mass components must "
            "apply exactly to ESO558-G009"
        )
    mass_systematic_sources = {
        row["source_key"]
        for row in rows
        if row["mass_systematic_status"]
        == "source_reported_model_selection_systematic"
    }
    if mass_systematic_sources != {"CGCG074-064"}:
        raise ValueError(
            "source-reported mass model-selection systematic must apply exactly "
            "to CGCG074-064"
        )
    position_mass_systematic_sources = {
        row["source_key"]
        for row in rows
        if row["mass_systematic_status"]
        == "source_reported_absolute_position_systematic"
    }
    if position_mass_systematic_sources != {"UGC6093"}:
        raise ValueError(
            "source-reported absolute-position mass systematic must apply "
            "exactly to UGC6093"
        )
    geometry_mass_systematic_sources = {
        row["source_key"]
        for row in rows
        if row["mass_systematic_status"]
        == "source_reported_midline_inclination_systematic"
    }
    if geometry_mass_systematic_sources != {"ESO558-G009"}:
        raise ValueError(
            "source-reported midline/inclination mass systematic must apply "
            "exactly to ESO558-G009"
        )
    reported_legacy_multimethod_distances = {
        row["source_key"]
        for row in rows
        if row["distance_stat_component_status"]
        == "reported_legacy_multimethod"
    }
    if reported_legacy_multimethod_distances != {"Circinus"}:
        raise ValueError(
            "reported legacy multi-indicator distance status must apply exactly to Circinus"
        )
    reported_legacy_tully_fisher_distances = {
        row["source_key"]
        for row in rows
        if row["distance_stat_component_status"]
        == "reported_legacy_tully_fisher"
    }
    if reported_legacy_tully_fisher_distances != {"NGC3079"}:
        raise ValueError(
            "reported legacy Tully--Fisher distance status must apply exactly "
            "to NGC3079"
        )

    circinus = next(row for row in rows if row["source_key"] == "Circinus")
    reviewed_circinus_numeric = {
        "reference_mass_msun": 1_700_000.0,
        "reference_mass_unc_minus_msun": 300_000.0,
        "reference_mass_unc_plus_msun": 300_000.0,
        "reference_distance_mpc": 4.2,
        "adopted_distance_mpc": 4.2,
        "distance_measurement_unc_minus_mpc": 0.8,
        "distance_measurement_unc_plus_mpc": 0.8,
        "size_stat_fraction_minus_override": 0.018867924528301886,
        "size_stat_fraction_plus_override": 0.018867924528301886,
    }
    for field, expected in reviewed_circinus_numeric.items():
        actual = optional_float(circinus[field])
        if actual != expected:
            raise ValueError(f"Circinus: {field} does not match the reviewed input")
    reviewed_circinus_text = {
        "treatment": "shared_distance",
        "systematic_status": "unquantified_model",
        "mass_reference_doi": "10.1086/374862",
        "mass_uncertainty_reference_doi": "10.1086/374862",
        "distance_reference_doi": "10.1086/374862",
        "uncertainty_note": "angular_radius_error_only_lower_limit_geometry",
        "mass_display_interval_status": "source_reported_total",
        "distance_stat_component_status": "reported_legacy_multimethod",
    }
    for field, expected in reviewed_circinus_text.items():
        if circinus[field] != expected:
            raise ValueError(f"Circinus: {field} does not match the reviewed input")

    ngc4388 = next(row for row in rows if row["source_key"] == "NGC4388")
    reviewed_ngc4388_numeric = {
        "reference_mass_msun": 8_400_000.0,
        "reference_mass_unc_minus_msun": 200_000.0,
        "reference_mass_unc_plus_msun": 200_000.0,
        "reference_distance_mpc": 19.0,
        "adopted_distance_mpc": 19.0,
        "distance_flow_unc_mpc": 0.0,
        "distance_shared_systematic_unc_mpc": 2.09,
    }
    for field, expected in reviewed_ngc4388_numeric.items():
        actual = optional_float(ngc4388[field])
        if actual != expected:
            raise ValueError(f"NGC4388: {field} does not match the reviewed input")
    for field in (
        "distance_measurement_unc_minus_mpc",
        "distance_measurement_unc_plus_mpc",
    ):
        if optional_float(ngc4388[field]) is not None:
            raise ValueError(
                f"NGC4388: {field} must be blank because the 11% distance term "
                "is a shared calibration systematic"
            )
    reviewed_ngc4388_text = {
        "treatment": "shared_distance",
        "systematic_status": "unquantified_model",
        "mass_reference_doi": "10.1088/0004-637X/727/1/20",
        "mass_uncertainty_reference_doi": "10.1088/0004-637X/727/1/20",
        "distance_reference_doi": "10.1088/0004-637X/727/1/20",
        "uncertainty_note": NGC4388_UNCERTAINTY_NOTE,
        "mass_display_interval_status": "decomposed_components",
        "distance_stat_component_status": "not_available",
    }
    for field, expected in reviewed_ngc4388_text.items():
        if ngc4388[field] != expected:
            raise ValueError(f"NGC4388: {field} does not match the reviewed input")

    ugc6093 = next(row for row in rows if row["source_key"] == "UGC6093")
    reviewed_ugc6093_numeric = {
        "reference_mass_msun": 25_800_000.0,
        "reference_mass_unc_minus_msun": 900_000.0,
        "reference_mass_unc_plus_msun": 1_100_000.0,
        "reference_mass_systematic_unc_minus_msun": 1_290_000.0,
        "reference_mass_systematic_unc_plus_msun": 1_290_000.0,
        "reference_distance_mpc": 151.7,
        "adopted_distance_mpc": 143.7235453315291,
        "distance_measurement_unc_minus_mpc": 0.16441136671177264,
        "distance_measurement_unc_plus_mpc": 0.14614343707713126,
        "distance_flow_unc_mpc": 3.3829499323410013,
        "distance_shared_systematic_unc_mpc": 5.8345146954612623,
        "size_systematic_fraction_minus": 0.05,
        "size_systematic_fraction_plus": 0.05,
    }
    for field, expected in reviewed_ugc6093_numeric.items():
        actual = optional_float(ugc6093[field])
        if actual != expected:
            raise ValueError(f"UGC6093: {field} does not match the reviewed input")
    for field in (
        "size_stat_fraction_minus_override",
        "size_stat_fraction_plus_override",
    ):
        if optional_float(ugc6093[field]) is not None:
            raise ValueError(f"UGC6093: {field} must remain blank")
    reviewed_ugc6093_text = {
        "source_order": "13",
        "source_label": "UGC 6093",
        "treatment": "shared_distance",
        "systematic_status": "numeric_systematic",
        "mass_reference_doi": "10.3847/1538-4357/aaa95c",
        "mass_uncertainty_reference_doi": "10.3847/1538-4357/aaa95c",
        "distance_reference_doi": "10.3847/1538-4357/aaa95c",
        "uncertainty_note": "absolute_position_systematic",
        "mass_display_interval_status": (
            "separate_decomposed_and_position_systematic_components"
        ),
        "distance_stat_component_status": "statistical",
        "mass_systematic_status": (
            "source_reported_absolute_position_systematic"
        ),
    }
    for field, expected in reviewed_ugc6093_text.items():
        if ugc6093[field] != expected:
            raise ValueError(f"UGC6093: {field} does not match the reviewed input")

    eso558 = next(row for row in rows if row["source_key"] == "ESO558-G009")
    reviewed_eso558_numeric = {
        "reference_mass_msun": 17_000_000.0,
        "reference_mass_systematic_unc_minus_msun": 1_020_566.5093466471,
        "reference_mass_systematic_unc_plus_msun": 1_020_566.5093466471,
        "reference_distance_mpc": 107.6,
        "adopted_distance_mpc": 106.28958051420837,
        "distance_measurement_unc_minus_mpc": 0.3653585926928281,
        "distance_measurement_unc_plus_mpc": 0.3653585926928281,
        "distance_flow_unc_mpc": 3.3829499323410013,
        "distance_shared_systematic_unc_mpc": 4.314867950509135,
        "size_stat_fraction_minus_override": 0.023,
        "size_stat_fraction_plus_override": 0.023,
        "size_systematic_fraction_minus": 0.06003332407921454,
        "size_systematic_fraction_plus": 0.06003332407921454,
    }
    for field, expected in reviewed_eso558_numeric.items():
        actual = optional_float(eso558[field])
        if actual != expected:
            raise ValueError(f"ESO558-G009: {field} does not match the reviewed input")
    for field in ("reference_mass_unc_minus_msun", "reference_mass_unc_plus_msun"):
        if optional_float(eso558[field]) is not None:
            raise ValueError(
                f"ESO558-G009: {field} must remain blank because the reviewed "
                "2.3 percent fit-plus-position fraction is explicit"
            )
    reviewed_eso558_text = {
        "source_order": "14",
        "source_label": "ESO 558-G009",
        "treatment": "shared_distance",
        "systematic_status": "numeric_systematic",
        "mass_reference_doi": "10.3847/1538-4357/834/1/52",
        "mass_uncertainty_reference_doi": "10.3847/1538-4357/834/1/52",
        "distance_reference_doi": "10.3847/1538-4357/834/1/52",
        "uncertainty_note": (
            "fit_plus_absolute_position_separate_from_midline_and_"
            "inclination_systematics"
        ),
        "mass_display_interval_status": (
            "separate_decomposed_and_geometry_systematic_components"
        ),
        "distance_stat_component_status": "statistical",
        "mass_systematic_status": (
            "source_reported_midline_inclination_systematic"
        ),
    }
    for field, expected in reviewed_eso558_text.items():
        if eso558[field] != expected:
            raise ValueError(f"ESO558-G009: {field} does not match the reviewed input")

    ngc2960 = next(row for row in rows if row["source_key"] == "NGC2960")
    reviewed_ngc2960_numeric = {
        "reference_mass_msun": 12_370_000.0,
        "reference_mass_unc_minus_msun": 60_000.0,
        "reference_mass_unc_plus_msun": 70_000.0,
        "reference_distance_mpc": 72.2,
        "adopted_distance_mpc": 71.32070365358592,
        "distance_measurement_unc_minus_mpc": 15.0 / 73.9,
        "distance_measurement_unc_plus_mpc": 15.0 / 73.9,
        "distance_flow_unc_mpc": 250.0 / 73.9,
        "distance_shared_systematic_unc_mpc": 2.8952924351929332,
    }
    for field, expected in reviewed_ngc2960_numeric.items():
        actual = optional_float(ngc2960[field])
        if actual != expected:
            raise ValueError(f"NGC2960: {field} does not match the reviewed input")
    for field in (
        "size_stat_fraction_minus_override",
        "size_stat_fraction_plus_override",
        "size_systematic_fraction_minus",
        "size_systematic_fraction_plus",
        "reference_mass_systematic_unc_minus_msun",
        "reference_mass_systematic_unc_plus_msun",
    ):
        if optional_float(ngc2960[field]) is not None:
            raise ValueError(f"NGC2960: {field} must remain blank")
    reviewed_ngc2960_text = {
        "source_order": "15",
        "source_label": "NGC 2960",
        "treatment": "shared_distance",
        "systematic_status": "unquantified_astrometric",
        "mass_reference_doi": "10.3847/1538-4357/aabff1",
        "mass_uncertainty_reference_doi": "10.3847/1538-4357/aabff1",
        "distance_reference_doi": "10.1088/0004-637X/727/1/20",
        "uncertainty_note": (
            "later_three_dimensional_disk_fit_formal_interval_only_absolute_"
            "coordinate_error_floor_completeness_unreported_legacy_three_to_"
            "five_percent_not_transplanted_disk_model_checks_not_adopted"
        ),
        "mass_display_interval_status": "decomposed_components",
        "distance_stat_component_status": "statistical",
        "mass_systematic_status": "not_available",
    }
    for field, expected in reviewed_ngc2960_text.items():
        if ngc2960[field] != expected:
            raise ValueError(f"NGC2960: {field} does not match the reviewed input")

    ic2560 = next(row for row in rows if row["source_key"] == "IC2560")
    reviewed_ic2560_numeric = {
        "reference_mass_msun": 4_140_000.0,
        "reference_mass_unc_minus_msun": 170_000.0,
        "reference_mass_unc_plus_msun": 170_000.0,
        "reference_distance_mpc": 44.5,
        "adopted_distance_mpc": 44.5,
        "distance_measurement_unc_minus_mpc": 6.1,
        "distance_measurement_unc_plus_mpc": 6.1,
        "distance_flow_unc_mpc": 0.0,
        "distance_shared_systematic_unc_mpc": 0.0,
    }
    for field, expected in reviewed_ic2560_numeric.items():
        actual = optional_float(ic2560[field])
        if actual != expected:
            raise ValueError(f"IC2560: {field} does not match the reviewed input")
    for field in (
        "size_stat_fraction_minus_override",
        "size_stat_fraction_plus_override",
        "size_systematic_fraction_minus",
        "size_systematic_fraction_plus",
        "reference_mass_systematic_unc_minus_msun",
        "reference_mass_systematic_unc_plus_msun",
    ):
        if optional_float(ic2560[field]) is not None:
            raise ValueError(f"IC2560: {field} must remain blank")
    reviewed_ic2560_text = {
        "source_order": "16",
        "source_label": "IC 2560",
        "treatment": "shared_distance",
        "systematic_status": "unquantified_model",
        "mass_reference_doi": "10.1017/S1743921323002971",
        "mass_uncertainty_reference_doi": "10.1017/S1743921323002971",
        "distance_reference_doi": "10.1017/S1743921323002971",
        "uncertainty_note": (
            "formal_fit_error_known_no_separately_quantified_model_"
            "systematic_available"
        ),
        "mass_display_interval_status": "decomposed_components",
        "distance_stat_component_status": "statistical",
        "mass_systematic_status": "not_available",
    }
    for field, expected in reviewed_ic2560_text.items():
        if ic2560[field] != expected:
            raise ValueError(f"IC2560: {field} does not match the reviewed input")

    ngc3079 = next(row for row in rows if row["source_key"] == "NGC3079")
    reviewed_ngc3079_numeric = {
        "reference_mass_msun": 2_000_000.0,
        "reference_distance_mpc": 17.3,
        "adopted_distance_mpc": 17.3,
        "distance_measurement_unc_minus_mpc": 1.5,
        "distance_measurement_unc_plus_mpc": 1.5,
        "distance_flow_unc_mpc": 0.0,
        "distance_shared_systematic_unc_mpc": 0.0,
    }
    for field, expected in reviewed_ngc3079_numeric.items():
        actual = optional_float(ngc3079[field])
        if actual != expected:
            raise ValueError(f"NGC3079: {field} does not match the reviewed input")
    for field in (
        "reference_mass_unc_minus_msun",
        "reference_mass_unc_plus_msun",
        "size_stat_fraction_minus_override",
        "size_stat_fraction_plus_override",
        "size_systematic_fraction_minus",
        "size_systematic_fraction_plus",
        "reference_mass_systematic_unc_minus_msun",
        "reference_mass_systematic_unc_plus_msun",
    ):
        if optional_float(ngc3079[field]) is not None:
            raise ValueError(f"NGC3079: {field} must remain blank")
    reviewed_ngc3079_text = {
        "source_order": "17",
        "source_label": "NGC 3079",
        "treatment": "shared_distance",
        "systematic_status": "unquantified_model",
        "mass_reference_doi": "10.1086/426101",
        "mass_uncertainty_reference_doi": "NA",
        "distance_reference_doi": "10.1086/426101",
        "uncertainty_note": (
            "approximate_enclosed_mass_black_hole_proxy_lower_limit_factor_of_"
            "order_unity_flat_rotation_disk_mass_significant_no_published_"
            "formal_interval"
        ),
        "mass_display_interval_status": "not_available",
        "distance_stat_component_status": "reported_legacy_tully_fisher",
        "mass_systematic_status": "not_available",
    }
    for field, expected in reviewed_ngc3079_text.items():
        if ngc3079[field] != expected:
            raise ValueError(f"NGC3079: {field} does not match the reviewed input")

    ngc5765b = next(row for row in rows if row["source_key"] == "NGC5765b")
    reviewed_ngc5765b_numeric = {
        "reference_mass_msun": 41_500_000.0,
        "reference_mass_unc_minus_msun": 1_900_000.0,
        "reference_mass_unc_plus_msun": 2_000_000.0,
        "reference_distance_mpc": 112.2,
        "adopted_distance_mpc": 112.2,
        "distance_measurement_unc_minus_mpc": 5.1,
        "distance_measurement_unc_plus_mpc": 5.4,
        "distance_flow_unc_mpc": 0.0,
        "distance_shared_systematic_unc_mpc": 0.0,
    }
    for field, expected in reviewed_ngc5765b_numeric.items():
        actual = optional_float(ngc5765b[field])
        if actual != expected:
            raise ValueError(f"NGC5765b: {field} does not match the reviewed input")
    for field in (
        "size_stat_fraction_minus_override",
        "size_stat_fraction_plus_override",
        "size_systematic_fraction_minus",
        "size_systematic_fraction_plus",
    ):
        if optional_float(ngc5765b[field]) is not None:
            raise ValueError(
                f"NGC5765b: {field} must be blank for the paired-chain treatment"
            )
    reviewed_ngc5765b_text = {
        "source_order": "7",
        "source_label": "NGC 5765b",
        "treatment": "paired_posterior",
        "systematic_status": "paired_quantile_file",
        "mass_reference_doi": "10.3847/2041-8213/ab75f0",
        "mass_uncertainty_reference_doi": "10.3847/2041-8213/ab75f0",
        "distance_reference_doi": "10.3847/2041-8213/ab75f0",
        "uncertainty_note": "paired_theta_g_quantiles_plus_one_percent_systematic",
        "mass_display_interval_status": "decomposed_components",
        "distance_stat_component_status": "statistical",
    }
    for field, expected in reviewed_ngc5765b_text.items():
        if ngc5765b[field] != expected:
            raise ValueError(f"NGC5765b: {field} does not match the reviewed input")

    ngc6323 = next(row for row in rows if row["source_key"] == "NGC6323")
    reviewed_ngc6323_numeric = {
        "reference_mass_msun": 10_200_000.0,
        "reference_mass_unc_minus_msun": 2_200_000.0,
        "reference_mass_unc_plus_msun": 3_200_000.0,
        "reference_distance_mpc": 109.4,
        "adopted_distance_mpc": 109.4,
        "distance_measurement_unc_minus_mpc": 23.4,
        "distance_measurement_unc_plus_mpc": 34.2,
        "distance_flow_unc_mpc": 0.0,
        "distance_shared_systematic_unc_mpc": 0.0,
    }
    for field, expected in reviewed_ngc6323_numeric.items():
        actual = optional_float(ngc6323[field])
        if actual != expected:
            raise ValueError(f"NGC6323: {field} does not match the reviewed input")
    for field in (
        "reference_mass_systematic_unc_minus_msun",
        "reference_mass_systematic_unc_plus_msun",
        "size_stat_fraction_minus_override",
        "size_stat_fraction_plus_override",
        "size_systematic_fraction_minus",
        "size_systematic_fraction_plus",
    ):
        if optional_float(ngc6323[field]) is not None:
            raise ValueError(
                f"NGC6323: {field} must be blank for the paired-chain treatment"
            )
    reviewed_ngc6323_text = {
        "source_order": "18",
        "source_label": "NGC 6323",
        "treatment": "paired_posterior",
        "systematic_status": "paired_quantile_file",
        "mass_reference_doi": "10.3847/2041-8213/ab75f0",
        "mass_uncertainty_reference_doi": "10.3847/2041-8213/ab75f0",
        "distance_reference_doi": "10.3847/2041-8213/ab75f0",
        "uncertainty_note": "paired_theta_g_quantiles_plus_one_percent_systematic",
        "mass_display_interval_status": "decomposed_components",
        "distance_stat_component_status": "statistical",
        "mass_systematic_status": "not_available",
    }
    for field, expected in reviewed_ngc6323_text.items():
        if ngc6323[field] != expected:
            raise ValueError(f"NGC6323: {field} does not match the reviewed input")

    j0437 = next(row for row in rows if row["source_key"] == "J0437")
    reviewed_j0437_numeric = {
        "reference_mass_msun": 2_860_000.0,
        "reference_mass_unc_minus_msun": 200_000.0,
        "reference_mass_unc_plus_msun": 200_000.0,
        "reference_distance_mpc": 70.0,
        "adopted_distance_mpc": 64.50473612990527,
        "distance_measurement_unc_minus_mpc": 40.0 / 73.9,
        "distance_measurement_unc_plus_mpc": 40.0 / 73.9,
        "distance_flow_unc_mpc": 250.0 / 73.9,
        "distance_shared_systematic_unc_mpc": 2.6185955127160461,
    }
    for field, expected in reviewed_j0437_numeric.items():
        actual = optional_float(j0437[field])
        if actual != expected:
            raise ValueError(f"J0437: {field} does not match the reviewed input")
    for field in (
        "size_stat_fraction_minus_override",
        "size_stat_fraction_plus_override",
        "size_systematic_fraction_minus",
        "size_systematic_fraction_plus",
        "reference_mass_systematic_unc_minus_msun",
        "reference_mass_systematic_unc_plus_msun",
    ):
        if optional_float(j0437[field]) is not None:
            raise ValueError(f"J0437: {field} must remain blank")
    reviewed_j0437_text = {
        "source_order": "19",
        "source_label": "J0437+2456",
        "treatment": "shared_distance",
        "systematic_status": "unquantified_model",
        "mass_reference_doi": "10.3847/1538-4357/abde3d",
        "mass_uncertainty_reference_doi": "10.3847/1538-4357/abde3d",
        "distance_reference_doi": "10.3847/1538-4357/834/1/52",
        "uncertainty_note": (
            "updated_nine_parameter_high_velocity_only_thin_disk_fit_formal_"
            "interval_shared_external_distance_inclination_unconstrained_and_"
            "model_completeness_unquantified"
        ),
        "mass_display_interval_status": "decomposed_components",
        "distance_stat_component_status": "statistical",
        "mass_systematic_status": "not_available",
    }
    for field, expected in reviewed_j0437_text.items():
        if j0437[field] != expected:
            raise ValueError(f"J0437: {field} does not match the reviewed input")

    mrk1029 = next(row for row in rows if row["source_key"] == "Mrk1029")
    reviewed_mrk1029_numeric = {
        "adopted_distance_mpc": 119.32882273342354,
        "distance_measurement_unc_minus_mpc": 32.0 / 73.9,
        "distance_measurement_unc_plus_mpc": 32.0 / 73.9,
        "distance_flow_unc_mpc": 250.0 / 73.9,
        "distance_shared_systematic_unc_mpc": 4.8442011935083977,
    }
    for field, expected in reviewed_mrk1029_numeric.items():
        actual = optional_float(mrk1029[field])
        if actual != expected:
            raise ValueError(f"Mrk1029: {field} does not match the reviewed input")
    for field in (
        "reference_mass_msun",
        "reference_mass_unc_minus_msun",
        "reference_mass_unc_plus_msun",
        "reference_distance_mpc",
        "size_stat_fraction_minus_override",
        "size_stat_fraction_plus_override",
        "size_systematic_fraction_minus",
        "size_systematic_fraction_plus",
        "reference_mass_systematic_unc_minus_msun",
        "reference_mass_systematic_unc_plus_msun",
    ):
        if optional_float(mrk1029[field]) is not None:
            raise ValueError(f"Mrk1029: {field} must remain blank")
    reviewed_mrk1029_text = {
        "source_order": "20",
        "source_label": "Mrk 1029",
        "treatment": "no_robust_mass",
        "systematic_status": "not_applicable",
        "mass_reference_doi": "10.3847/1538-4357/834/1/52",
        "mass_uncertainty_reference_doi": "NA",
        "distance_reference_doi": "10.3847/1538-4357/834/1/52",
        "uncertainty_note": (
            "no_robust_primary_mass_huchra1999_9076_plusminus32_optical_"
            "velocity_component_selected_amid_two_velocity_ambiguity"
        ),
        "mass_display_interval_status": "not_available",
        "distance_stat_component_status": "statistical",
        "mass_systematic_status": "not_available",
    }
    for field, expected in reviewed_mrk1029_text.items():
        if mrk1029[field] != expected:
            raise ValueError(f"Mrk1029: {field} does not match the reviewed input")

    ngc1320 = next(row for row in rows if row["source_key"] == "NGC1320")
    reviewed_ngc1320_numeric = {
        "adopted_distance_mpc": 34.2 * 73.0 / 73.9,
        "distance_measurement_unc_minus_mpc": 16.0 / 73.9,
        "distance_measurement_unc_plus_mpc": 16.0 / 73.9,
        "distance_flow_unc_mpc": 250.0 / 73.9,
        "distance_shared_systematic_unc_mpc": (34.2 * 73.0 / 73.9) * 3.0 / 73.9,
    }
    for field, expected in reviewed_ngc1320_numeric.items():
        actual = optional_float(ngc1320[field])
        if actual != expected:
            raise ValueError(f"NGC1320: {field} does not match the reviewed input")
    for field in (
        "reference_mass_msun",
        "reference_mass_unc_minus_msun",
        "reference_mass_unc_plus_msun",
        "reference_distance_mpc",
        "size_stat_fraction_minus_override",
        "size_stat_fraction_plus_override",
        "size_systematic_fraction_minus",
        "size_systematic_fraction_plus",
        "reference_mass_systematic_unc_minus_msun",
        "reference_mass_systematic_unc_plus_msun",
    ):
        if optional_float(ngc1320[field]) is not None:
            raise ValueError(f"NGC1320: {field} must remain blank")
    reviewed_ngc1320_text = {
        "source_order": "21",
        "source_label": "NGC 1320",
        "treatment": "no_robust_mass",
        "systematic_status": "not_applicable",
        "mass_reference_doi": "10.3847/1538-4357/834/1/52",
        "mass_uncertainty_reference_doi": "NA",
        "distance_reference_doi": "10.3847/1538-4357/834/1/52",
        "uncertainty_note": (
            "no_robust_primary_mass_gao2017_2649_plusminus16_optical_"
            "velocity_component"
        ),
        "mass_display_interval_status": "not_available",
        "distance_stat_component_status": "statistical",
        "mass_systematic_status": "not_available",
    }
    for field, expected in reviewed_ngc1320_text.items():
        if ngc1320[field] != expected:
            raise ValueError(f"NGC1320: {field} does not match the reviewed input")

    cgcg074064 = next(row for row in rows if row["source_key"] == "CGCG074-064")
    reviewed_cgcg074064_numeric = {
        "reference_mass_msun": 24_200_000.0,
        "reference_mass_unc_minus_msun": 2_000_000.0,
        "reference_mass_unc_plus_msun": 2_200_000.0,
        "reference_mass_systematic_unc_minus_msun": 500_000.0,
        "reference_mass_systematic_unc_plus_msun": 500_000.0,
        "reference_distance_mpc": 87.6,
        "adopted_distance_mpc": 87.6,
        "distance_measurement_unc_minus_mpc": 7.2,
        "distance_measurement_unc_plus_mpc": 7.9,
        "distance_flow_unc_mpc": 0.0,
        "distance_shared_systematic_unc_mpc": 1.5,
    }
    for field, expected in reviewed_cgcg074064_numeric.items():
        actual = optional_float(cgcg074064[field])
        if actual != expected:
            raise ValueError(
                f"CGCG074-064: {field} does not match the reviewed input"
            )
    for field in (
        "size_stat_fraction_minus_override",
        "size_stat_fraction_plus_override",
        "size_systematic_fraction_minus",
        "size_systematic_fraction_plus",
    ):
        if optional_float(cgcg074064[field]) is not None:
            raise ValueError(
                f"CGCG074-064: {field} must be blank for the paired-chain treatment"
            )
    reviewed_cgcg074064_text = {
        "source_order": "10",
        "source_label": "CGCG 074-064",
        "treatment": "paired_posterior",
        "systematic_status": "paired_quantile_file",
        "mass_systematic_status": "source_reported_model_selection_systematic",
        "mass_reference_doi": "10.3847/1538-4357/ab6bcd",
        "mass_uncertainty_reference_doi": "10.3847/1538-4357/ab6bcd",
        "distance_reference_doi": "10.3847/1538-4357/ab6bcd",
        "uncertainty_note": "paired_theta_g_quantiles_plus_one_percent_systematic",
        "mass_display_interval_status": (
            "separate_statistical_and_systematic_components"
        ),
        "distance_stat_component_status": "statistical",
    }
    for field, expected in reviewed_cgcg074064_text.items():
        if cgcg074064[field] != expected:
            raise ValueError(
                f"CGCG074-064: {field} does not match the reviewed input"
            )

    return rows


QUANTILE_FIELDS = [
    "source_key",
    "q16_theta_g_nas",
    "q50_theta_g_nas",
    "q84_theta_g_nas",
    "systematic_fraction",
    "provenance",
]
REVIEWED_PAIRED_QUANTILES = {
    "NGC5765b": {
        "q16_theta_g_nas": 3.6373,
        "q50_theta_g_nas": 3.6493,
        "q84_theta_g_nas": 3.6623,
        "systematic_fraction": 0.01,
        "provenance": "D_Pesce_private_communication_2026-08-17",
    },
    "CGCG074-064": {
        "q16_theta_g_nas": 2.7101,
        "q50_theta_g_nas": 2.7277,
        "q84_theta_g_nas": 2.7474,
        "systematic_fraction": 0.01,
        "provenance": "D_Pesce_private_communication_2026-08-17",
    },
    "NGC6323": {
        "q16_theta_g_nas": 0.9135,
        "q50_theta_g_nas": 0.9189,
        "q84_theta_g_nas": 0.9248,
        "systematic_fraction": 0.01,
        "provenance": "D_Pesce_private_communication_2026-08-17",
    },
}

MODEL_CHECK_FIELDS = [
    "source_key",
    "check_key",
    "check_type",
    "data_relation",
    "range_basis",
    "endpoint_derivation",
    "sensitivity_direction",
    "accepted_configuration_count",
    "excluded_configuration_ids",
    "lower_configuration_id",
    "lower_mass_msun",
    "lower_distance_mpc",
    "upper_configuration_id",
    "upper_mass_msun",
    "upper_distance_mpc",
    "lower_primary_ratio_factor",
    "upper_primary_ratio_factor",
    "assumed_inner_radius_pc",
    "assumed_turbulent_velocity_kms",
    "assumed_surface_density_power_law_index",
    "assumed_rotation_profile",
    "mass_bias_eta",
    "assumption_profile",
    "probability_status",
    "adoption_status",
    "reference_key",
    "reference_doi",
    "provenance_note",
]

MODEL_CHECK_TYPES = {
    "data_selection_component_tracking",
    "disk_mass_model_choice",
    "turbulent_support_mass_bias",
}
MODEL_CHECK_DATA_RELATIONS = {
    "subset_reanalysis_of_argon2007_vlbi_epochs_used_by_primary",
    "published_analytic_model_applied_to_circinus_disk_parameters",
    "published_analytic_model_applied_to_ngc4388_disk_parameters",
    "published_analytic_model_applied_to_ngc3079_disk_parameters",
    "published_analytic_model_applied_to_ngc4945_disk_parameters",
    "predecessor_fit_data_alternative_accretion_disk_models",
    "reanalysis_of_published_gao2017_vlbi_data_data_selection_unreported",
    "same_fit_data_alternative_accretion_disk_models",
}
MODEL_CHECK_RANGE_BASES = {
    "accepted_configuration_medians",
    "one_sided_fractional_mass_underestimate",
    "published_physically_allowed_model_envelope",
    "reported_symmetric_mass_error_span",
}
MODEL_CHECK_ENDPOINT_DERIVATIONS = {
    "published_mass_distance_endpoints",
    "adopted_primary_ratio_factors",
    "reported_mbh_plusminus_error_at_fixed_reported_distance",
}
MODEL_CHECK_SENSITIVITY_DIRECTIONS = {
    "downward_only",
    "two_sided_range",
    "upward_only",
}
MODEL_CHECK_ROTATION_PROFILES = {"flat_rotation", "keplerian_rotation"}
MODEL_CHECK_ASSUMPTION_PROFILES = {
    "constant_turbulent_velocity_and_turbulent_pressure_dominates_thermal",
    "edge_on_coplanar_black_hole_plus_mestel_disk",
}
MODEL_CHECK_PROBABILITY_STATUSES = {"non_probabilistic"}
MODEL_CHECK_ADOPTION_STATUSES = {"reported_separately_not_adopted"}
MODEL_CHECK_REFERENCE_DOIS = {
    "Kuo2018Disk": "10.3847/1538-4357/aabff1",
    "vanDerBoom2026": "10.1093/mnras/stag1327",
    "Kuo2020": "10.1093/mnras/staa2260",
    "Deepshikha2024": "10.1017/S1743921323002491",
}

AUDIT_ONLY_MODEL_CHECK_FIELDS = [
    "source_key",
    "check_key",
    "check_type",
    "data_relation",
    "range_basis",
    "endpoint_derivation",
    "sensitivity_direction",
    "lower_reference_ratio_factor",
    "upper_reference_ratio_factor",
    "assumed_inner_radius_pc",
    "assumed_turbulent_velocity_kms",
    "assumed_surface_density_power_law_index",
    "assumed_rotation_profile",
    "mass_bias_eta",
    "assumption_profile",
    "probability_status",
    "adoption_status",
    "attachment_status",
    "reference_key",
    "reference_doi",
    "provenance_note",
]

# This result is scientifically relevant but cannot be attached to the
# generated NGC1320 row: Gao et al. did not establish a robust primary mass to
# which the one-percent conditional correction could be applied.  The exact
# record is nevertheless required during every deterministic build.
REVIEWED_AUDIT_ONLY_MODEL_CHECKS = {
    "kuo_2020_ngc1320_turbulent_support_bias": {
        "source_key": "NGC1320",
        "check_type": "turbulent_support_mass_bias",
        "data_relation": (
            "published_analytic_model_applied_to_ngc1320_indicative_disk_"
            "parameters"
        ),
        "range_basis": "one_sided_fractional_mass_underestimate",
        "endpoint_derivation": "conditional_reference_ratio_factors",
        "sensitivity_direction": "upward_only",
        "lower_reference_ratio_factor": 1.0,
        "upper_reference_ratio_factor": 1.01,
        "assumed_inner_radius_pc": 0.07,
        "assumed_turbulent_velocity_kms": 40.0,
        "assumed_surface_density_power_law_index": -1.0,
        "assumed_rotation_profile": "keplerian_rotation",
        "mass_bias_eta": 2.5,
        "assumption_profile": (
            "constant_turbulent_velocity_and_turbulent_pressure_dominates_"
            "thermal"
        ),
        "probability_status": "non_probabilistic",
        "adoption_status": "reported_separately_not_adopted",
        "attachment_status": "audit_only_unattached_no_robust_primary",
        "reference_key": "Kuo2020",
        "reference_doi": "10.1093/mnras/staa2260",
        "provenance_note": (
            "appendix_a_ngc1320_rin_0p07pc_vturb_40kms_sigma_r_minus1_"
            "keplerian_eta_2p5_delta_m_turb_over_m_enc_0p01_conditional_"
            "on_gao2017_indicative_enclosed_mass_no_adopted_primary"
        ),
    }
}

REVIEWED_MODEL_CHECKS = {
    "NGC4258": {
        "check_key": "van_der_boom_kushnir_2026_component_tracking",
        "check_type": "data_selection_component_tracking",
        "data_relation": "subset_reanalysis_of_argon2007_vlbi_epochs_used_by_primary",
        "range_basis": "accepted_configuration_medians",
        "endpoint_derivation": "published_mass_distance_endpoints",
        "sensitivity_direction": "two_sided_range",
        "accepted_configuration_count": 15,
        "excluded_configuration_ids": (2, 17),
        "lower_configuration_id": 1,
        "lower_mass_msun": 42_340_000.0,
        "lower_distance_mpc": 8.091,
        "upper_configuration_id": 8,
        "upper_mass_msun": 38_430_000.0,
        "upper_distance_mpc": 7.182,
        "lower_primary_ratio_factor": None,
        "upper_primary_ratio_factor": None,
        "assumed_inner_radius_pc": None,
        "assumed_turbulent_velocity_kms": None,
        "assumed_surface_density_power_law_index": None,
        "assumed_rotation_profile": None,
        "mass_bias_eta": None,
        "assumption_profile": None,
        "probability_status": "non_probabilistic",
        "adoption_status": "reported_separately_not_adopted",
        "reference_key": "vanDerBoom2026",
        "reference_doi": "10.1093/mnras/stag1327",
        "provenance_note": "derived_from_published_table_2_marginal_medians_using_iau_2015_shadow_conversion",
    },
    "NGC3079": {
        "check_key": "kuo_2020_ngc3079_turbulent_support_bias",
        "check_type": "turbulent_support_mass_bias",
        "data_relation": "published_analytic_model_applied_to_ngc3079_disk_parameters",
        "range_basis": "one_sided_fractional_mass_underestimate",
        "endpoint_derivation": "adopted_primary_ratio_factors",
        "sensitivity_direction": "upward_only",
        "accepted_configuration_count": None,
        "excluded_configuration_ids": None,
        "lower_configuration_id": None,
        "lower_mass_msun": None,
        "lower_distance_mpc": None,
        "upper_configuration_id": None,
        "upper_mass_msun": None,
        "upper_distance_mpc": None,
        "lower_primary_ratio_factor": 1.0,
        "upper_primary_ratio_factor": 1.34,
        "assumed_inner_radius_pc": 0.40,
        "assumed_turbulent_velocity_kms": 60.0,
        "assumed_surface_density_power_law_index": -1.0,
        "assumed_rotation_profile": "flat_rotation",
        "mass_bias_eta": 2.0,
        "assumption_profile": (
            "constant_turbulent_velocity_and_turbulent_pressure_dominates_thermal"
        ),
        "probability_status": "non_probabilistic",
        "adoption_status": "reported_separately_not_adopted",
        "reference_key": "Kuo2020",
        "reference_doi": "10.1093/mnras/staa2260",
        "provenance_note": "appendix_a_delta_m_turb_over_m_enc_0p34_under_stated_assumptions",
    },
    "NGC4388": {
        "check_key": "kuo_2020_ngc4388_turbulent_support_bias",
        "check_type": "turbulent_support_mass_bias",
        "data_relation": "published_analytic_model_applied_to_ngc4388_disk_parameters",
        "range_basis": "one_sided_fractional_mass_underestimate",
        "endpoint_derivation": "adopted_primary_ratio_factors",
        "sensitivity_direction": "upward_only",
        "accepted_configuration_count": None,
        "excluded_configuration_ids": None,
        "lower_configuration_id": None,
        "lower_mass_msun": None,
        "lower_distance_mpc": None,
        "upper_configuration_id": None,
        "upper_mass_msun": None,
        "upper_distance_mpc": None,
        "lower_primary_ratio_factor": 1.0,
        "upper_primary_ratio_factor": 1.001,
        "assumed_inner_radius_pc": 0.24,
        "assumed_turbulent_velocity_kms": 10.0,
        "assumed_surface_density_power_law_index": -1.0,
        "assumed_rotation_profile": "keplerian_rotation",
        "mass_bias_eta": 2.5,
        "assumption_profile": (
            "constant_turbulent_velocity_and_turbulent_pressure_dominates_thermal"
        ),
        "probability_status": "non_probabilistic",
        "adoption_status": "reported_separately_not_adopted",
        "reference_key": "Kuo2020",
        "reference_doi": "10.1093/mnras/staa2260",
        "provenance_note": (
            "appendix_a_delta_m_turb_over_m_enc_0p001_under_stated_assumptions"
        ),
    },
    "NGC4945": {
        "check_key": "kuo_2020_ngc4945_turbulent_support_bias",
        "check_type": "turbulent_support_mass_bias",
        "data_relation": "published_analytic_model_applied_to_ngc4945_disk_parameters",
        "range_basis": "one_sided_fractional_mass_underestimate",
        "endpoint_derivation": "adopted_primary_ratio_factors",
        "sensitivity_direction": "upward_only",
        "accepted_configuration_count": None,
        "excluded_configuration_ids": None,
        "lower_configuration_id": None,
        "lower_mass_msun": None,
        "lower_distance_mpc": None,
        "upper_configuration_id": None,
        "upper_mass_msun": None,
        "upper_distance_mpc": None,
        "lower_primary_ratio_factor": 1.0,
        "upper_primary_ratio_factor": 1.24,
        "assumed_inner_radius_pc": 0.16,
        "assumed_turbulent_velocity_kms": 60.0,
        "assumed_surface_density_power_law_index": -1.0,
        "assumed_rotation_profile": "keplerian_rotation",
        "mass_bias_eta": 2.5,
        "assumption_profile": (
            "constant_turbulent_velocity_and_turbulent_pressure_dominates_thermal"
        ),
        "probability_status": "non_probabilistic",
        "adoption_status": "reported_separately_not_adopted",
        "reference_key": "Kuo2020",
        "reference_doi": "10.1093/mnras/staa2260",
        "provenance_note": "appendix_a_delta_m_turb_over_m_enc_0p24_under_stated_assumptions",
    },
    "Circinus": {
        "check_key": "kuo_2020_circinus_turbulent_support_bias",
        "check_type": "turbulent_support_mass_bias",
        "data_relation": "published_analytic_model_applied_to_circinus_disk_parameters",
        "range_basis": "one_sided_fractional_mass_underestimate",
        "endpoint_derivation": "adopted_primary_ratio_factors",
        "sensitivity_direction": "upward_only",
        "accepted_configuration_count": None,
        "excluded_configuration_ids": None,
        "lower_configuration_id": None,
        "lower_mass_msun": None,
        "lower_distance_mpc": None,
        "upper_configuration_id": None,
        "upper_mass_msun": None,
        "upper_distance_mpc": None,
        "lower_primary_ratio_factor": 1.0,
        "upper_primary_ratio_factor": 1.03,
        "assumed_inner_radius_pc": 0.11,
        "assumed_turbulent_velocity_kms": 30.0,
        "assumed_surface_density_power_law_index": -1.0,
        "assumed_rotation_profile": "keplerian_rotation",
        "mass_bias_eta": 2.5,
        "assumption_profile": (
            "constant_turbulent_velocity_and_turbulent_pressure_dominates_thermal"
        ),
        "probability_status": "non_probabilistic",
        "adoption_status": "reported_separately_not_adopted",
        "reference_key": "Kuo2020",
        "reference_doi": "10.1093/mnras/staa2260",
        "provenance_note": "appendix_a_delta_m_turb_over_m_enc_0p03_under_stated_assumptions",
    },
    "NGC2960": {
        "check_key": "kuo_2018_ngc2960_physically_allowed_disk_mass_models",
        "check_type": "disk_mass_model_choice",
        "data_relation": "same_fit_data_alternative_accretion_disk_models",
        "range_basis": "published_physically_allowed_model_envelope",
        "endpoint_derivation": "published_mass_distance_endpoints",
        "sensitivity_direction": "downward_only",
        "accepted_configuration_count": 3,
        "excluded_configuration_ids": (2, 4),
        "lower_configuration_id": 3,
        "lower_mass_msun": 12_320_000.0,
        "lower_distance_mpc": 72.2,
        "upper_configuration_id": 1,
        "upper_mass_msun": 12_370_000.0,
        "upper_distance_mpc": 72.2,
        "lower_primary_ratio_factor": None,
        "upper_primary_ratio_factor": None,
        "assumed_inner_radius_pc": None,
        "assumed_turbulent_velocity_kms": None,
        "assumed_surface_density_power_law_index": None,
        "assumed_rotation_profile": None,
        "mass_bias_eta": None,
        "assumption_profile": None,
        "probability_status": "non_probabilistic",
        "adoption_status": "reported_separately_not_adopted",
        "reference_key": "Kuo2018Disk",
        "reference_doi": "10.3847/1538-4357/aabff1",
        "provenance_note": (
            "table_2_configuration_ids_1_fiducial_2_linear_no_density_excluded_"
            "3_linear_physical_lower_4_herrnstein_no_density_excluded_5_"
            "herrnstein_physical_intermediate_upper_endpoint_fiducial"
        ),
    },
    "NGC6323": {
        "check_key": "kuo_2018_ngc6323_physically_allowed_disk_mass_models",
        "check_type": "disk_mass_model_choice",
        "data_relation": "predecessor_fit_data_alternative_accretion_disk_models",
        "range_basis": "published_physically_allowed_model_envelope",
        "endpoint_derivation": "published_mass_distance_endpoints",
        "sensitivity_direction": "downward_only",
        "accepted_configuration_count": 3,
        "excluded_configuration_ids": (2, 4),
        "lower_configuration_id": 3,
        "lower_mass_msun": 9_840_000.0,
        "lower_distance_mpc": 106.7,
        "upper_configuration_id": 1,
        "upper_mass_msun": 9_920_000.0,
        "upper_distance_mpc": 106.7,
        "lower_primary_ratio_factor": None,
        "upper_primary_ratio_factor": None,
        "assumed_inner_radius_pc": None,
        "assumed_turbulent_velocity_kms": None,
        "assumed_surface_density_power_law_index": None,
        "assumed_rotation_profile": None,
        "mass_bias_eta": None,
        "assumption_profile": None,
        "probability_status": "non_probabilistic",
        "adoption_status": "reported_separately_not_adopted",
        "reference_key": "Kuo2018Disk",
        "reference_doi": "10.3847/1538-4357/aabff1",
        "provenance_note": (
            "table_2_configuration_ids_1_fiducial_2_linear_no_density_excluded_"
            "3_linear_physical_lower_4_herrnstein_no_density_excluded_5_"
            "herrnstein_physical_intermediate_upper_endpoint_fiducial"
        ),
    },
    "J0437": {
        "check_key": "deepshikha_2024_j0437_mestel_disk_mass_model",
        "check_type": "disk_mass_model_choice",
        "data_relation": (
            "reanalysis_of_published_gao2017_vlbi_data_data_selection_unreported"
        ),
        "range_basis": "reported_symmetric_mass_error_span",
        "endpoint_derivation": (
            "reported_mbh_plusminus_error_at_fixed_reported_distance"
        ),
        "sensitivity_direction": "downward_only",
        "accepted_configuration_count": 1,
        "excluded_configuration_ids": None,
        "lower_configuration_id": None,
        "lower_mass_msun": 1_620_000.0,
        "lower_distance_mpc": 68.27,
        "upper_configuration_id": None,
        "upper_mass_msun": 2_080_000.0,
        "upper_distance_mpc": 68.27,
        "lower_primary_ratio_factor": None,
        "upper_primary_ratio_factor": None,
        "assumed_inner_radius_pc": None,
        "assumed_turbulent_velocity_kms": None,
        "assumed_surface_density_power_law_index": -1.0,
        "assumed_rotation_profile": None,
        "mass_bias_eta": None,
        "assumption_profile": "edge_on_coplanar_black_hole_plus_mestel_disk",
        "probability_status": "non_probabilistic",
        "adoption_status": "reported_separately_not_adopted",
        "reference_key": "Deepshikha2024",
        "reference_doi": "10.1017/S1743921323002491",
        "provenance_note": (
            "official_supplement_table_1_abs_alpha_0p303_plusminus0p037_mbh_"
            "1p85_plusminus0p23e6_msun_mdisk_0p82_plusminus0p26e6_msun_"
            "distance_68p27_plusminus0p45_mpc_main_text_edge_on_mestel_"
            "assumption_distance_imported_from_literature_exact_spot_selection_"
            "and_error_confidence_semantics_unreported_ratio_endpoints_hold_"
            "reported_distance_fixed"
        ),
    },
    "Mrk1029": {
        "check_key": "deepshikha_2024_mrk1029_mestel_disk_mass_model",
        "check_type": "disk_mass_model_choice",
        "data_relation": (
            "reanalysis_of_published_gao2017_vlbi_data_data_selection_unreported"
        ),
        "range_basis": "reported_symmetric_mass_error_span",
        "endpoint_derivation": (
            "reported_mbh_plusminus_error_at_fixed_reported_distance"
        ),
        "sensitivity_direction": "two_sided_range",
        "accepted_configuration_count": 1,
        "excluded_configuration_ids": None,
        "lower_configuration_id": None,
        "lower_mass_msun": 410_000.0,
        "lower_distance_mpc": 126.787,
        "upper_configuration_id": None,
        "upper_mass_msun": 730_000.0,
        "upper_distance_mpc": 126.787,
        "lower_primary_ratio_factor": None,
        "upper_primary_ratio_factor": None,
        "assumed_inner_radius_pc": None,
        "assumed_turbulent_velocity_kms": None,
        "assumed_surface_density_power_law_index": -1.0,
        "assumed_rotation_profile": None,
        "mass_bias_eta": None,
        "assumption_profile": "edge_on_coplanar_black_hole_plus_mestel_disk",
        "probability_status": "non_probabilistic",
        "adoption_status": "reported_separately_not_adopted",
        "reference_key": "Deepshikha2024",
        "reference_doi": "10.1017/S1743921323002491",
        "provenance_note": (
            "official_supplement_table_1_abs_alpha_0p310_plusminus0p065_mbh_"
            "0p57_plusminus0p16e6_msun_mdisk_0p93_plusminus0p65e6_msun_"
            "distance_126p787_plusminus0p059_mpc_main_text_edge_on_mestel_"
            "assumption_distance_imported_from_literature_exact_spot_selection_"
            "and_error_confidence_semantics_unreported_ratio_endpoints_hold_"
            "reported_distance_fixed_no_adopted_primary_for_directional_"
            "comparison"
        ),
    },
    "NGC1320": {
        "check_key": "deepshikha_2024_ngc1320_mestel_disk_mass_model",
        "check_type": "disk_mass_model_choice",
        "data_relation": (
            "reanalysis_of_published_gao2017_vlbi_data_data_selection_unreported"
        ),
        "range_basis": "reported_symmetric_mass_error_span",
        "endpoint_derivation": (
            "reported_mbh_plusminus_error_at_fixed_reported_distance"
        ),
        "sensitivity_direction": "two_sided_range",
        "accepted_configuration_count": 1,
        "excluded_configuration_ids": None,
        "lower_configuration_id": None,
        "lower_mass_msun": 4_130_000.0,
        "lower_distance_mpc": 38.4,
        "upper_configuration_id": None,
        "upper_mass_msun": 5_290_000.0,
        "upper_distance_mpc": 38.4,
        "lower_primary_ratio_factor": None,
        "upper_primary_ratio_factor": None,
        "assumed_inner_radius_pc": None,
        "assumed_turbulent_velocity_kms": None,
        "assumed_surface_density_power_law_index": -1.0,
        "assumed_rotation_profile": None,
        "mass_bias_eta": None,
        "assumption_profile": "edge_on_coplanar_black_hole_plus_mestel_disk",
        "probability_status": "non_probabilistic",
        "adoption_status": "reported_separately_not_adopted",
        "reference_key": "Deepshikha2024",
        "reference_doi": "10.1017/S1743921323002491",
        "provenance_note": (
            "official_supplement_table_1_abs_alpha_0p479_plusminus0p078_mbh_"
            "4p71_plusminus0p58e6_msun_mdisk_0p10_plusminus1p14e6_msun_"
            "distance_38p4_plusminus0p44_mpc_main_text_edge_on_mestel_"
            "assumption_distance_imported_from_literature_exact_spot_selection_"
            "and_error_confidence_semantics_unreported_ratio_endpoints_hold_"
            "reported_distance_fixed_no_adopted_primary_for_directional_"
            "comparison"
        ),
    },
}
EXPECTED_MODEL_CHECK_SOURCES = set(REVIEWED_MODEL_CHECKS)
MODEL_CHECK_PRIMARY_POLICIES = {
    "NGC4258": {
        "treatment": "paired_posterior",
        "size_total_status": "quadrature_statistical_plus_systematic",
        "primary_result_status": "retained",
    },
    "NGC3079": {
        "treatment": "shared_distance",
        "size_total_status": "incomplete_unquantified_systematic",
        "primary_result_status": "retained",
    },
    "NGC4388": {
        "treatment": "shared_distance",
        "size_total_status": "incomplete_unquantified_systematic",
        "primary_result_status": "retained",
    },
    "NGC4945": {
        "treatment": "shared_distance",
        "size_total_status": "no_formal_uncertainty",
        "primary_result_status": "retained",
    },
    "Circinus": {
        "treatment": "shared_distance",
        "size_total_status": "incomplete_unquantified_systematic",
        "primary_result_status": "retained",
    },
    "NGC2960": {
        "treatment": "shared_distance",
        "size_total_status": "incomplete_unquantified_systematic",
        "primary_result_status": "retained",
    },
    "NGC6323": {
        "treatment": "paired_posterior",
        "size_total_status": "quadrature_statistical_plus_systematic",
        "primary_result_status": "retained",
    },
    "J0437": {
        "treatment": "shared_distance",
        "size_total_status": "incomplete_unquantified_systematic",
        "primary_result_status": "retained",
    },
    "Mrk1029": {
        "treatment": "no_robust_mass",
        "size_total_status": "not_available",
        "primary_result_status": "absent",
    },
    "NGC1320": {
        "treatment": "no_robust_mass",
        "size_total_status": "not_available",
        "primary_result_status": "absent",
    },
}
if set(MODEL_CHECK_PRIMARY_POLICIES) != EXPECTED_MODEL_CHECK_SOURCES:
    raise AssertionError("model-check policies must match the reviewed source set")

NONADOPTED_FIT_FIELDS = [
    "source_key",
    "fit_key",
    "data_relation",
    "model_basis",
    "ratio_msun_per_mpc",
    "ratio_unc_minus_msun_per_mpc",
    "ratio_unc_plus_msun_per_mpc",
    "interval_status",
    "adoption_status",
    "nonadoption_basis",
    "reference_key",
    "reference_doi",
    "provenance_note",
]
NONADOPTED_FIT_DATA_RELATIONS = {
    "reanalysis_of_published_gao2016_vlbi_spot_data",
    "reanalysis_of_published_gao2017_vlbi_spot_data",
    "reanalysis_of_published_kuo2011_vlbi_spot_data",
    "reanalysis_of_published_pesce2020a_vlbi_spot_data",
    "reanalysis_of_published_reid2013_vlbi_spot_data",
}
NONADOPTED_FIT_MODEL_BASES = {"circular_equatorial_schwarzschild_disk"}
NONADOPTED_FIT_INTERVAL_STATUSES = {"bayesian_formal_conditional_on_model"}
NONADOPTED_FIT_ADOPTION_STATUSES = {"reported_separately_not_adopted"}
NONADOPTED_FIT_BASES = {
    "does_not_resolve_primary_complex_geometry_and_dynamics",
    (
        "restricted_same_data_fit_conditional_on_fixed_gao2017_inclination_"
        "vertical_center_and_selected_scatter_angle_does_not_resolve_primary_"
        "complex_geometry_and_dynamics"
    ),
    (
        "restricted_same_data_fit_conditional_on_fixed_edge_on_inclination_"
        "and_selected_scatter_angle_does_not_resolve_primary_complex_"
        "geometry_and_dynamics"
    ),
    (
        "older_partial_data_reanalysis_conditional_on_fixed_vertical_center_and_"
        "selected_scatter_angle"
    ),
    (
        "same_data_reanalysis_conditional_on_fixed_gao2016_inclination_and_"
        "selected_scatter_angle"
    ),
    (
        "same_data_reanalysis_conditional_on_fixed_gao2017_inclination_and_"
        "selected_scatter_angle"
    ),
    (
        "same_data_reanalysis_conditional_on_fixed_gao2017_inclination_"
        "vertical_center_and_zero_scatter_angle"
    ),
    (
        "older_published_data_reanalysis_conditional_on_fixed_kuo2011_inclination_"
        "vertical_center_and_selected_scatter_angle"
    ),
    (
        "older_partial_data_reanalysis_conditional_on_fixed_kuo2015_inclination_"
        "and_selected_scatter_angle"
    ),
    (
        "same_data_reanalysis_conditional_on_fixed_pesce2020a_inclination_and_"
        "selected_scatter_angle"
    ),
    (
        "same_observation_vlbi_spot_subset_reanalysis_conditional_on_fixed_"
        "reid2013_inclination_and_selected_scatter_angle"
    ),
    "same_data_reanalysis_conditional_on_fixed_vertical_center_and_selected_scatter_angle",
    "same_data_reanalysis_conditional_on_restricted_disk_model_and_discrepant_recession_velocity",
}
NONADOPTED_FIT_REFERENCE_DOIS = {
    "GonzalezJuarez2024": "10.1051/0004-6361/202450098",
    "Villaraos2022": "10.1093/mnras/stac2973",
}
REVIEWED_NONADOPTED_FITS = {
    "NGC1194": {
        "fit_key": "gonzalez_juarez_2024_ngc1194_gr_ratio_fit",
        "data_relation": "reanalysis_of_published_kuo2011_vlbi_spot_data",
        "model_basis": "circular_equatorial_schwarzschild_disk",
        "ratio_msun_per_mpc": 1_310_000.0,
        "ratio_unc_minus_msun_per_mpc": 20_900.0,
        "ratio_unc_plus_msun_per_mpc": 21_100.0,
        "interval_status": "bayesian_formal_conditional_on_model",
        "adoption_status": "reported_separately_not_adopted",
        "nonadoption_basis": (
            "same_data_reanalysis_conditional_on_restricted_disk_model_and_"
            "discrepant_recession_velocity"
        ),
        "reference_key": "GonzalezJuarez2024",
        "reference_doi": "10.1051/0004-6361/202450098",
        "provenance_note": "published_table_2_ratio_and_formal_interval",
    },
    "NGC2273": {
        "fit_key": "villaraos_2022_ngc2273_gr_ratio_fit",
        "data_relation": "reanalysis_of_published_kuo2011_vlbi_spot_data",
        "model_basis": "circular_equatorial_schwarzschild_disk",
        "ratio_msun_per_mpc": 328_300.0,
        "ratio_unc_minus_msun_per_mpc": 5_100.0,
        "ratio_unc_plus_msun_per_mpc": 4_500.0,
        "interval_status": "bayesian_formal_conditional_on_model",
        "adoption_status": "reported_separately_not_adopted",
        "nonadoption_basis": (
            "same_data_reanalysis_conditional_on_fixed_vertical_center_and_"
            "selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
        "reference_doi": "10.1093/mnras/stac2973",
        "provenance_note": "published_table_2_ratio_and_formal_interval",
    },
    "NGC4388": {
        "fit_key": "villaraos_2022_ngc4388_gr_ratio_fit",
        "data_relation": "reanalysis_of_published_kuo2011_vlbi_spot_data",
        "model_basis": "circular_equatorial_schwarzschild_disk",
        "ratio_msun_per_mpc": 423_300.0,
        "ratio_unc_minus_msun_per_mpc": 26_800.0,
        "ratio_unc_plus_msun_per_mpc": 18_500.0,
        "interval_status": "bayesian_formal_conditional_on_model",
        "adoption_status": "reported_separately_not_adopted",
        "nonadoption_basis": (
            "same_data_reanalysis_conditional_on_fixed_vertical_center_and_"
            "selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
        "reference_doi": "10.1093/mnras/stac2973",
        "provenance_note": "published_table_2_ratio_and_formal_interval",
    },
    "CGCG074-064": {
        "fit_key": "villaraos_2022_cgcg074064_gr_ratio_fit",
        "data_relation": "reanalysis_of_published_pesce2020a_vlbi_spot_data",
        "model_basis": "circular_equatorial_schwarzschild_disk",
        "ratio_msun_per_mpc": 274_900.0,
        "ratio_unc_minus_msun_per_mpc": 1_500.0,
        "ratio_unc_plus_msun_per_mpc": 1_700.0,
        "interval_status": "bayesian_formal_conditional_on_model",
        "adoption_status": "reported_separately_not_adopted",
        "nonadoption_basis": (
            "same_data_reanalysis_conditional_on_fixed_pesce2020a_inclination_and_"
            "selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
        "reference_doi": "10.1093/mnras/stac2973",
        "provenance_note": (
            "published_table_2_ratio_and_formal_interval_table_2_scatter_11deg_"
            "chi2red_1p289_conflicts_with_section_4p1p4_scatter_10deg_chi2red_"
            "1p600_no_scatter_or_chi2_value_adopted"
        ),
    },
    "NGC6264": {
        "fit_key": "villaraos_2022_ngc6264_gr_ratio_fit",
        "data_relation": "reanalysis_of_published_kuo2011_vlbi_spot_data",
        "model_basis": "circular_equatorial_schwarzschild_disk",
        "ratio_msun_per_mpc": 212_600.0,
        "ratio_unc_minus_msun_per_mpc": 700.0,
        "ratio_unc_plus_msun_per_mpc": 700.0,
        "interval_status": "bayesian_formal_conditional_on_model",
        "adoption_status": "reported_separately_not_adopted",
        "nonadoption_basis": (
            "older_partial_data_reanalysis_conditional_on_fixed_vertical_center_and_"
            "selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
        "reference_doi": "10.1093/mnras/stac2973",
        "provenance_note": (
            "published_table_2_ratio_and_formal_interval_table_2_scatter_10deg_"
            "chi2red_1p304_conflicts_with_section_4p1p7_chi2red_1p301"
        ),
    },
    "UGC3789": {
        "fit_key": "villaraos_2022_ugc3789_gr_ratio_fit",
        "data_relation": "reanalysis_of_published_reid2013_vlbi_spot_data",
        "model_basis": "circular_equatorial_schwarzschild_disk",
        "ratio_msun_per_mpc": 227_700.0,
        "ratio_unc_minus_msun_per_mpc": 900.0,
        "ratio_unc_plus_msun_per_mpc": 900.0,
        "interval_status": "bayesian_formal_conditional_on_model",
        "adoption_status": "reported_separately_not_adopted",
        "nonadoption_basis": (
            "same_observation_vlbi_spot_subset_reanalysis_conditional_on_fixed_"
            "reid2013_inclination_and_selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
        "reference_doi": "10.1093/mnras/stac2973",
        "provenance_note": (
            "published_table_2_ratio_and_formal_interval_table_2_and_section_"
            "4p1p3_scatter_10deg_chi2red_1p430"
        ),
    },
    "ESO558-G009": {
        "fit_key": "villaraos_2022_eso558g009_gr_ratio_fit",
        "data_relation": "reanalysis_of_published_gao2017_vlbi_spot_data",
        "model_basis": "circular_equatorial_schwarzschild_disk",
        "ratio_msun_per_mpc": 159_400.0,
        "ratio_unc_minus_msun_per_mpc": 2_600.0,
        "ratio_unc_plus_msun_per_mpc": 2_500.0,
        "interval_status": "bayesian_formal_conditional_on_model",
        "adoption_status": "reported_separately_not_adopted",
        "nonadoption_basis": (
            "same_data_reanalysis_conditional_on_fixed_gao2017_inclination_and_"
            "selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
        "reference_doi": "10.1093/mnras/stac2973",
        "provenance_note": (
            "published_table_2_ratio_and_formal_interval_table_2_and_section_"
            "4p1p5_fixed_inclination_98deg_scatter_7deg_chi2red_1p188"
        ),
    },
    "NGC2960": {
        "fit_key": "villaraos_2022_ngc2960_gr_ratio_fit",
        "data_relation": "reanalysis_of_published_kuo2011_vlbi_spot_data",
        "model_basis": "circular_equatorial_schwarzschild_disk",
        "ratio_msun_per_mpc": 173_800.0,
        "ratio_unc_minus_msun_per_mpc": 1_600.0,
        "ratio_unc_plus_msun_per_mpc": 1_500.0,
        "interval_status": "bayesian_formal_conditional_on_model",
        "adoption_status": "reported_separately_not_adopted",
        "nonadoption_basis": (
            "older_published_data_reanalysis_conditional_on_fixed_kuo2011_inclination_"
            "vertical_center_and_selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
        "reference_doi": "10.1093/mnras/stac2973",
        "provenance_note": (
            "published_table_2_ratio_and_formal_interval_fixed_i89_y0_"
            "minus0p020_scatter15deg_table_chi2red1p488_conflicts_with_"
            "section_4p1p6_y0_origin_x0_0p004_chi2red1p511_section_4p1p6_"
            "y0_variation_up_to_4_percent_not_in_formal_interval_exact_spot_"
            "overlap_with_kuo2018_unverifiable"
        ),
    },
    "NGC5765b": {
        "fit_key": "villaraos_2022_ngc5765b_gr_ratio_fit",
        "data_relation": "reanalysis_of_published_gao2016_vlbi_spot_data",
        "model_basis": "circular_equatorial_schwarzschild_disk",
        "ratio_msun_per_mpc": 372_700.0,
        "ratio_unc_minus_msun_per_mpc": 1_300.0,
        "ratio_unc_plus_msun_per_mpc": 1_300.0,
        "interval_status": "bayesian_formal_conditional_on_model",
        "adoption_status": "reported_separately_not_adopted",
        "nonadoption_basis": (
            "same_data_reanalysis_conditional_on_fixed_gao2016_inclination_and_"
            "selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
        "reference_doi": "10.1093/mnras/stac2973",
        "provenance_note": "published_table_2_ratio_and_formal_interval",
    },
    "NGC6323": {
        "fit_key": "villaraos_2022_ngc6323_gr_ratio_fit",
        "data_relation": "reanalysis_of_published_kuo2011_vlbi_spot_data",
        "model_basis": "circular_equatorial_schwarzschild_disk",
        "ratio_msun_per_mpc": 91_600.0,
        "ratio_unc_minus_msun_per_mpc": 500.0,
        "ratio_unc_plus_msun_per_mpc": 500.0,
        "interval_status": "bayesian_formal_conditional_on_model",
        "adoption_status": "reported_separately_not_adopted",
        "nonadoption_basis": (
            "older_partial_data_reanalysis_conditional_on_fixed_kuo2015_"
            "inclination_and_selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
        "reference_doi": "10.1093/mnras/stac2973",
        "provenance_note": (
            "published_table_2_ratio_and_formal_interval_four_parameters_y0_"
            "fitted_table_2_scatter_11deg_chi2red_1p276_conflicts_with_"
            "section_4p1p2_chi2red_1p251"
        ),
    },
    "J0437": {
        "fit_key": "villaraos_2022_j0437_gr_ratio_fit",
        "data_relation": "reanalysis_of_published_gao2017_vlbi_spot_data",
        "model_basis": "circular_equatorial_schwarzschild_disk",
        "ratio_msun_per_mpc": 42_500.0,
        "ratio_unc_minus_msun_per_mpc": 1_200.0,
        "ratio_unc_plus_msun_per_mpc": 1_000.0,
        "interval_status": "bayesian_formal_conditional_on_model",
        "adoption_status": "reported_separately_not_adopted",
        "nonadoption_basis": (
            "same_data_reanalysis_conditional_on_fixed_gao2017_inclination_"
            "vertical_center_and_zero_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
        "reference_doi": "10.1093/mnras/stac2973",
        "provenance_note": (
            "published_table_2_ratio_and_formal_interval_three_parameters_"
            "fixed_i81_y0_minus0p040_zero_scatter_chi2red_1p038_section_"
            "4p1p8_y0_variation_up_to_2_percent_within_formal_interval"
        ),
    },
    "Mrk1029": {
        "fit_key": "gonzalez_juarez_2024_mrk1029_gr_ratio_fit",
        "data_relation": "reanalysis_of_published_gao2017_vlbi_spot_data",
        "model_basis": "circular_equatorial_schwarzschild_disk",
        "ratio_msun_per_mpc": 14_400.0,
        "ratio_unc_minus_msun_per_mpc": 1_100.0,
        "ratio_unc_plus_msun_per_mpc": 1_100.0,
        "interval_status": "bayesian_formal_conditional_on_model",
        "adoption_status": "reported_separately_not_adopted",
        "nonadoption_basis": (
            "restricted_same_data_fit_conditional_on_fixed_gao2017_inclination_"
            "vertical_center_and_selected_scatter_angle_does_not_resolve_primary_"
            "complex_geometry_and_dynamics"
        ),
        "reference_key": "GonzalezJuarez2024",
        "reference_doi": "10.1051/0004-6361/202450098",
        "provenance_note": (
            "published_table_2_ratio_and_formal_interval_three_parameters_m_over_"
            "d_zrec_x0_fitted_y0_fixed_i79_fixed_scatter14deg_chi2red1p385_pa_"
            "plusminus10deg_sensitivity_1p19percent_neglected_within_formal_interval"
        ),
    },
    "NGC1320": {
        "fit_key": "gonzalez_juarez_2024_ngc1320_gr_ratio_fit",
        "data_relation": "reanalysis_of_published_gao2017_vlbi_spot_data",
        "model_basis": "circular_equatorial_schwarzschild_disk",
        "ratio_msun_per_mpc": 149_700.0,
        "ratio_unc_minus_msun_per_mpc": 6_900.0,
        "ratio_unc_plus_msun_per_mpc": 7_100.0,
        "interval_status": "bayesian_formal_conditional_on_model",
        "adoption_status": "reported_separately_not_adopted",
        "nonadoption_basis": (
            "restricted_same_data_fit_conditional_on_fixed_edge_on_inclination_"
            "and_selected_scatter_angle_does_not_resolve_primary_complex_"
            "geometry_and_dynamics"
        ),
        "reference_key": "GonzalezJuarez2024",
        "reference_doi": "10.1051/0004-6361/202450098",
        "provenance_note": (
            "published_table_2_ratio_and_formal_interval_four_parameters_m_over_"
            "d_zrec_x0_y0_fitted_reference_midpoint_between_redshifted_and_"
            "blueshifted_group_geometric_centers_i90_fixed_scatter17deg_"
            "chi2red1p127_pa_plusminus10deg_sensitivity_2p64percent_within_"
            "formal_interval_bh_position_near_blueshifted_group_new_vlbi_"
            "requested"
        ),
    },
}
NONADOPTED_FIT_PRIMARY_POLICIES = {
    "NGC1194": {
        "treatment": "shared_distance",
        "primary_result_status": "retained",
    },
    "NGC2273": {
        "treatment": "shared_distance",
        "primary_result_status": "retained",
    },
    "NGC4388": {
        "treatment": "shared_distance",
        "primary_result_status": "retained",
    },
    "CGCG074-064": {
        "treatment": "paired_posterior",
        "primary_result_status": "retained",
    },
    "NGC6264": {
        "treatment": "paired_posterior",
        "primary_result_status": "retained",
    },
    "UGC3789": {
        "treatment": "paired_posterior",
        "primary_result_status": "retained",
    },
    "ESO558-G009": {
        "treatment": "shared_distance",
        "primary_result_status": "retained",
    },
    "NGC2960": {
        "treatment": "shared_distance",
        "primary_result_status": "retained",
    },
    "NGC5765b": {
        "treatment": "paired_posterior",
        "primary_result_status": "retained",
    },
    "NGC6323": {
        "treatment": "paired_posterior",
        "primary_result_status": "retained",
    },
    "J0437": {
        "treatment": "shared_distance",
        "primary_result_status": "retained",
    },
    "Mrk1029": {
        "treatment": "no_robust_mass",
        "primary_result_status": "absent",
    },
    "NGC1320": {
        "treatment": "no_robust_mass",
        "primary_result_status": "absent",
    },
}
EXPECTED_NONADOPTED_FIT_SOURCES = set(REVIEWED_NONADOPTED_FITS)
if set(NONADOPTED_FIT_PRIMARY_POLICIES) != EXPECTED_NONADOPTED_FIT_SOURCES:
    raise AssertionError("non-adopted fit policies must match the reviewed source set")


def load_paired_quantiles(path: Path = DEFAULT_QUANTILES) -> dict[str, dict[str, object]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        _validate_header(reader.fieldnames, QUANTILE_FIELDS, path)
        raw_rows = list(reader)
    if any(None in row for row in raw_rows):
        raise ValueError(f"one or more rows in {path} contain extra columns")

    rows: dict[str, dict[str, object]] = {}
    for raw in raw_rows:
        key = raw["source_key"]
        if key in rows:
            raise ValueError(f"duplicate paired-quantile source: {key}")
        q16 = float(raw["q16_theta_g_nas"])
        q50 = float(raw["q50_theta_g_nas"])
        q84 = float(raw["q84_theta_g_nas"])
        systematic_fraction = float(raw["systematic_fraction"])
        if not 0.0 < q16 < q50 < q84:
            raise ValueError(f"{key}: paired theta_g quantiles must be positive and ordered")
        if not 0.0 <= systematic_fraction < 1.0:
            raise ValueError(f"{key}: invalid systematic fraction")
        if not raw["provenance"]:
            raise ValueError(f"{key}: provenance is required")
        rows[key] = {
            "q16_theta_g_nas": q16,
            "q50_theta_g_nas": q50,
            "q84_theta_g_nas": q84,
            "systematic_fraction": systematic_fraction,
            "provenance": raw["provenance"],
        }
        reviewed = REVIEWED_PAIRED_QUANTILES.get(key)
        if reviewed is not None:
            for field, expected in reviewed.items():
                if rows[key][field] != expected:
                    raise ValueError(
                        f"{key}: {field} does not match the reviewed paired quantile input"
                    )
    return rows


def _positive_int(value: str, *, field: str, source_key: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"{source_key}: {field} must be an integer") from exc
    if parsed <= 0:
        raise ValueError(f"{source_key}: {field} must be positive")
    return parsed


def _require_literal_na(raw: dict[str, str], fields: tuple[str, ...], source_key: str) -> None:
    for field in fields:
        if raw[field] != "NA":
            raise ValueError(
                f"{source_key}: {field} must be literal NA for this model-check branch"
            )


def load_model_checks(
    path: Path = DEFAULT_MODEL_CHECKS,
    *,
    known_sources: set[str] | None = None,
) -> dict[str, dict[str, object]]:
    """Load non-probabilistic checks that remain outside adopted error budgets."""

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        _validate_header(reader.fieldnames, MODEL_CHECK_FIELDS, path)
        raw_rows = list(reader)
    if any(None in row for row in raw_rows):
        raise ValueError(f"one or more rows in {path} contain extra columns")
    if not raw_rows:
        raise ValueError(f"at least one model-sensitivity check is required in {path}")

    rows: dict[str, dict[str, object]] = {}
    check_keys: set[str] = set()
    for raw in raw_rows:
        source_key = raw["source_key"]
        check_key = raw["check_key"]
        if not source_key or not check_key:
            raise ValueError("model-sensitivity source_key and check_key are required")
        if known_sources is not None and source_key not in known_sources:
            raise ValueError(f"unknown model-sensitivity source: {source_key}")
        if source_key in rows:
            raise ValueError(
                f"multiple model-sensitivity checks for {source_key} require an "
                "explicit aggregation policy"
            )
        if check_key in check_keys:
            raise ValueError(f"duplicate model-sensitivity check_key: {check_key}")

        check_type = raw["check_type"]
        data_relation = raw["data_relation"]
        range_basis = raw["range_basis"]
        endpoint_derivation = raw["endpoint_derivation"]
        sensitivity_direction = raw["sensitivity_direction"]
        probability_status = raw["probability_status"]
        adoption_status = raw["adoption_status"]
        if check_type not in MODEL_CHECK_TYPES:
            raise ValueError(f"{source_key}: unsupported model-check type {check_type!r}")
        if data_relation not in MODEL_CHECK_DATA_RELATIONS:
            raise ValueError(
                f"{source_key}: unsupported model-check data relation {data_relation!r}"
            )
        if range_basis not in MODEL_CHECK_RANGE_BASES:
            raise ValueError(f"{source_key}: unsupported model-check range basis {range_basis!r}")
        if endpoint_derivation not in MODEL_CHECK_ENDPOINT_DERIVATIONS:
            raise ValueError(
                f"{source_key}: unsupported model-check endpoint derivation "
                f"{endpoint_derivation!r}"
            )
        if sensitivity_direction not in MODEL_CHECK_SENSITIVITY_DIRECTIONS:
            raise ValueError(
                f"{source_key}: unsupported model-check sensitivity direction "
                f"{sensitivity_direction!r}"
            )
        if probability_status not in MODEL_CHECK_PROBABILITY_STATUSES:
            raise ValueError(
                f"{source_key}: unsupported model-check probability status "
                f"{probability_status!r}"
            )
        if adoption_status not in MODEL_CHECK_ADOPTION_STATUSES:
            raise ValueError(
                f"{source_key}: unsupported model-check adoption status {adoption_status!r}"
            )

        if endpoint_derivation == "published_mass_distance_endpoints":
            _require_literal_na(
                raw,
                (
                    "lower_primary_ratio_factor",
                    "upper_primary_ratio_factor",
                    "assumed_inner_radius_pc",
                    "assumed_turbulent_velocity_kms",
                    "assumed_surface_density_power_law_index",
                    "assumed_rotation_profile",
                    "mass_bias_eta",
                    "assumption_profile",
                ),
                source_key,
            )
            accepted_count = _positive_int(
                raw["accepted_configuration_count"],
                field="accepted_configuration_count",
                source_key=source_key,
            )
            excluded_tokens = raw["excluded_configuration_ids"].split(";")
            if not excluded_tokens or any(
                not token.isdigit() for token in excluded_tokens
            ):
                raise ValueError(f"{source_key}: malformed excluded_configuration_ids")
            excluded_ids: tuple[int, ...] | None = tuple(
                int(token) for token in excluded_tokens
            )
            if len(set(excluded_ids)) != len(excluded_ids) or any(
                identifier <= 0 for identifier in excluded_ids
            ):
                raise ValueError(f"{source_key}: invalid excluded_configuration_ids")

            lower_configuration = _positive_int(
                raw["lower_configuration_id"],
                field="lower_configuration_id",
                source_key=source_key,
            )
            upper_configuration = _positive_int(
                raw["upper_configuration_id"],
                field="upper_configuration_id",
                source_key=source_key,
            )
            if lower_configuration == upper_configuration:
                raise ValueError(
                    f"{source_key}: model-check endpoints require distinct configurations"
                )
            if (
                lower_configuration in excluded_ids
                or upper_configuration in excluded_ids
            ):
                raise ValueError(
                    f"{source_key}: model-check endpoint cannot be an excluded fit"
                )

            lower_mass = optional_float(raw["lower_mass_msun"])
            lower_distance = optional_float(raw["lower_distance_mpc"])
            upper_mass = optional_float(raw["upper_mass_msun"])
            upper_distance = optional_float(raw["upper_distance_mpc"])
            if any(
                value is None or value <= 0.0
                for value in (lower_mass, lower_distance, upper_mass, upper_distance)
            ):
                raise ValueError(
                    f"{source_key}: model-check endpoint masses/distances must be positive"
                )
            lower_size = shadow_size_uas(float(lower_mass), float(lower_distance))
            upper_size = shadow_size_uas(float(upper_mass), float(upper_distance))
            if lower_size >= upper_size:
                raise ValueError(f"{source_key}: model-check size bounds must be ordered")

            lower_factor = None
            upper_factor = None
            inner_radius = None
            turbulent_velocity = None
            surface_density_index = None
            rotation_profile = None
            mass_bias_eta = None
            assumption_profile = None
        elif endpoint_derivation == (
            "reported_mbh_plusminus_error_at_fixed_reported_distance"
        ):
            _require_literal_na(
                raw,
                (
                    "excluded_configuration_ids",
                    "lower_configuration_id",
                    "upper_configuration_id",
                    "lower_primary_ratio_factor",
                    "upper_primary_ratio_factor",
                    "assumed_inner_radius_pc",
                    "assumed_turbulent_velocity_kms",
                    "assumed_rotation_profile",
                    "mass_bias_eta",
                ),
                source_key,
            )
            accepted_count = _positive_int(
                raw["accepted_configuration_count"],
                field="accepted_configuration_count",
                source_key=source_key,
            )
            if accepted_count != 1:
                raise ValueError(
                    f"{source_key}: a reported mass-error span requires one model"
                )
            excluded_ids = None
            lower_configuration = None
            upper_configuration = None
            lower_mass = optional_float(raw["lower_mass_msun"])
            lower_distance = optional_float(raw["lower_distance_mpc"])
            upper_mass = optional_float(raw["upper_mass_msun"])
            upper_distance = optional_float(raw["upper_distance_mpc"])
            if any(
                value is None or value <= 0.0
                for value in (lower_mass, lower_distance, upper_mass, upper_distance)
            ):
                raise ValueError(
                    f"{source_key}: reported mass-error endpoints must be positive"
                )
            if lower_distance != upper_distance:
                raise ValueError(
                    f"{source_key}: reported mass-error endpoints must hold distance fixed"
                )
            lower_size = shadow_size_uas(float(lower_mass), float(lower_distance))
            upper_size = shadow_size_uas(float(upper_mass), float(upper_distance))
            if lower_size >= upper_size:
                raise ValueError(
                    f"{source_key}: reported mass-error size bounds must be ordered"
                )
            lower_factor = None
            upper_factor = None
            inner_radius = None
            turbulent_velocity = None
            surface_density_index = optional_float(
                raw["assumed_surface_density_power_law_index"]
            )
            if surface_density_index != -1.0:
                raise ValueError(
                    f"{source_key}: the reviewed Mestel-disk index must be -1"
                )
            rotation_profile = None
            mass_bias_eta = None
            assumption_profile = raw["assumption_profile"]
            if assumption_profile != "edge_on_coplanar_black_hole_plus_mestel_disk":
                raise ValueError(
                    f"{source_key}: the reviewed Mestel-disk assumptions do not match"
                )
        else:
            _require_literal_na(
                raw,
                (
                    "accepted_configuration_count",
                    "excluded_configuration_ids",
                    "lower_configuration_id",
                    "lower_mass_msun",
                    "lower_distance_mpc",
                    "upper_configuration_id",
                    "upper_mass_msun",
                    "upper_distance_mpc",
                ),
                source_key,
            )
            accepted_count = None
            excluded_ids = None
            lower_configuration = None
            lower_mass = None
            lower_distance = None
            upper_configuration = None
            upper_mass = None
            upper_distance = None
            lower_size = None
            upper_size = None

            lower_factor = optional_float(raw["lower_primary_ratio_factor"])
            upper_factor = optional_float(raw["upper_primary_ratio_factor"])
            if (
                lower_factor is None
                or upper_factor is None
                or lower_factor <= 0.0
                or lower_factor >= upper_factor
            ):
                raise ValueError(
                    f"{source_key}: primary-ratio factors must be positive and ordered"
                )
            if check_type == "turbulent_support_mass_bias":
                inner_radius = optional_float(raw["assumed_inner_radius_pc"])
                turbulent_velocity = optional_float(
                    raw["assumed_turbulent_velocity_kms"]
                )
                surface_density_index = optional_float(
                    raw["assumed_surface_density_power_law_index"]
                )
                if inner_radius is None or inner_radius <= 0.0:
                    raise ValueError(
                        f"{source_key}: assumed_inner_radius_pc must be positive"
                    )
                if turbulent_velocity is None or turbulent_velocity <= 0.0:
                    raise ValueError(
                        f"{source_key}: assumed_turbulent_velocity_kms must be positive"
                    )
                if surface_density_index is None:
                    raise ValueError(
                        f"{source_key}: surface-density power-law index is required"
                    )
                rotation_profile = raw["assumed_rotation_profile"]
                if rotation_profile in {"", "NA"}:
                    raise ValueError(
                        f"{source_key}: assumed_rotation_profile is required"
                    )
                if rotation_profile not in MODEL_CHECK_ROTATION_PROFILES:
                    raise ValueError(
                        f"{source_key}: unsupported model-check rotation profile"
                    )
                mass_bias_eta = optional_float(raw["mass_bias_eta"])
                if mass_bias_eta is None:
                    raise ValueError(f"{source_key}: mass_bias_eta is required")
                if mass_bias_eta <= 0.0:
                    raise ValueError(f"{source_key}: mass_bias_eta must be positive")
                assumption_profile = raw["assumption_profile"]
                if assumption_profile not in MODEL_CHECK_ASSUMPTION_PROFILES:
                    raise ValueError(
                        f"{source_key}: unsupported model-check assumption profile"
                    )
            elif check_type == "disk_mass_model_choice":
                _require_literal_na(
                    raw,
                    (
                        "assumed_inner_radius_pc",
                        "assumed_turbulent_velocity_kms",
                        "assumed_surface_density_power_law_index",
                        "assumed_rotation_profile",
                        "mass_bias_eta",
                        "assumption_profile",
                    ),
                    source_key,
                )
                inner_radius = None
                turbulent_velocity = None
                surface_density_index = None
                rotation_profile = None
                mass_bias_eta = None
                assumption_profile = None
            else:
                raise ValueError(
                    f"{source_key}: ratio-factor endpoint is unsupported for "
                    f"model-check type {check_type!r}"
                )

        reference_key = raw["reference_key"]
        reference_doi = raw["reference_doi"]
        if reference_key not in MODEL_CHECK_REFERENCE_DOIS:
            raise ValueError(f"{source_key}: unsupported model-check reference_key")
        if reference_doi != MODEL_CHECK_REFERENCE_DOIS[reference_key]:
            raise ValueError(
                f"{source_key}: model-check reference key does not match its DOI"
            )
        if not raw["provenance_note"]:
            raise ValueError(f"{source_key}: model-check provenance_note is required")

        reviewed = REVIEWED_MODEL_CHECKS.get(source_key)
        if reviewed is None:
            raise ValueError(f"unsupported reviewed model-sensitivity source: {source_key}")
        reviewed_values: dict[str, object] = {
            "check_key": check_key,
            "check_type": check_type,
            "data_relation": data_relation,
            "range_basis": range_basis,
            "endpoint_derivation": endpoint_derivation,
            "sensitivity_direction": sensitivity_direction,
            "accepted_configuration_count": accepted_count,
            "excluded_configuration_ids": excluded_ids,
            "lower_configuration_id": lower_configuration,
            "lower_mass_msun": lower_mass,
            "lower_distance_mpc": lower_distance,
            "upper_configuration_id": upper_configuration,
            "upper_mass_msun": upper_mass,
            "upper_distance_mpc": upper_distance,
            "lower_primary_ratio_factor": lower_factor,
            "upper_primary_ratio_factor": upper_factor,
            "assumed_inner_radius_pc": inner_radius,
            "assumed_turbulent_velocity_kms": turbulent_velocity,
            "assumed_surface_density_power_law_index": surface_density_index,
            "assumed_rotation_profile": rotation_profile,
            "mass_bias_eta": mass_bias_eta,
            "assumption_profile": assumption_profile,
            "probability_status": probability_status,
            "adoption_status": adoption_status,
            "reference_key": reference_key,
            "reference_doi": reference_doi,
            "provenance_note": raw["provenance_note"],
        }
        for field, expected in reviewed.items():
            if reviewed_values[field] != expected:
                raise ValueError(
                    f"{source_key}: {field} does not match the reviewed model check"
                )

        rows[source_key] = {
            "source_key": source_key,
            "check_key": check_key,
            "check_type": check_type,
            "data_relation": data_relation,
            "range_basis": range_basis,
            "endpoint_derivation": endpoint_derivation,
            "sensitivity_direction": sensitivity_direction,
            "accepted_configuration_count": accepted_count,
            "excluded_configuration_ids": (
                ";".join(str(value) for value in excluded_ids)
                if excluded_ids is not None
                else None
            ),
            "lower_configuration_id": lower_configuration,
            "lower_mass_msun": lower_mass,
            "lower_distance_mpc": lower_distance,
            "upper_configuration_id": upper_configuration,
            "upper_mass_msun": upper_mass,
            "upper_distance_mpc": upper_distance,
            "lower_primary_ratio_factor": lower_factor,
            "upper_primary_ratio_factor": upper_factor,
            "assumed_inner_radius_pc": inner_radius,
            "assumed_turbulent_velocity_kms": turbulent_velocity,
            "assumed_surface_density_power_law_index": surface_density_index,
            "assumed_rotation_profile": rotation_profile,
            "mass_bias_eta": mass_bias_eta,
            "assumption_profile": assumption_profile,
            "size_lower_uas": lower_size,
            "size_upper_uas": upper_size,
            "probability_status": probability_status,
            "adoption_status": adoption_status,
            "reference_key": reference_key,
            "reference_doi": reference_doi,
            "provenance_note": raw["provenance_note"],
        }
        check_keys.add(check_key)
    return rows


def load_audit_only_model_checks(
    path: Path = DEFAULT_AUDIT_ONLY_MODEL_CHECKS,
    *,
    known_sources: set[str] | None = None,
) -> dict[str, dict[str, object]]:
    """Validate checks retained as evidence but forbidden from output attachment."""

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        _validate_header(reader.fieldnames, AUDIT_ONLY_MODEL_CHECK_FIELDS, path)
        raw_rows = list(reader)
    if any(None in row for row in raw_rows):
        raise ValueError(f"one or more rows in {path} contain extra columns")

    numeric_fields = (
        "lower_reference_ratio_factor",
        "upper_reference_ratio_factor",
        "assumed_inner_radius_pc",
        "assumed_turbulent_velocity_kms",
        "assumed_surface_density_power_law_index",
        "mass_bias_eta",
    )
    rows: dict[str, dict[str, object]] = {}
    sources: set[str] = set()
    for raw in raw_rows:
        source_key = raw["source_key"]
        check_key = raw["check_key"]
        if not source_key or not check_key:
            raise ValueError(
                "audit-only model-sensitivity source_key and check_key are required"
            )
        if known_sources is not None and source_key not in known_sources:
            raise ValueError(f"unknown audit-only model-sensitivity source: {source_key}")
        if source_key in sources:
            raise ValueError(
                f"multiple audit-only model-sensitivity checks for {source_key} "
                "require an explicit aggregation policy"
            )
        if check_key in rows:
            raise ValueError(f"duplicate audit-only model-sensitivity check_key: {check_key}")

        reviewed = REVIEWED_AUDIT_ONLY_MODEL_CHECKS.get(check_key)
        if reviewed is None:
            raise ValueError(f"unsupported audit-only model-sensitivity check: {check_key}")
        candidate: dict[str, object] = {
            field: raw[field]
            for field in AUDIT_ONLY_MODEL_CHECK_FIELDS
            if field not in numeric_fields and field != "check_key"
        }
        for field in numeric_fields:
            try:
                value = float(raw[field])
            except ValueError as error:
                raise ValueError(
                    f"{check_key}: {field} must be a finite number"
                ) from error
            if not math.isfinite(value):
                raise ValueError(f"{check_key}: {field} must be a finite number")
            candidate[field] = value

        for field, expected in reviewed.items():
            if candidate[field] != expected:
                raise ValueError(
                    f"{check_key}: {field} does not match the reviewed audit-only check"
                )
        if not (
            0.0
            < float(candidate["lower_reference_ratio_factor"])
            < float(candidate["upper_reference_ratio_factor"])
        ):
            raise ValueError(f"{check_key}: conditional ratio factors must be ordered")
        if candidate["attachment_status"] != (
            "audit_only_unattached_no_robust_primary"
        ):
            raise ValueError(f"{check_key}: audit-only attachment status is required")

        rows[check_key] = {"check_key": check_key, **candidate}
        sources.add(source_key)

    if set(rows) != set(REVIEWED_AUDIT_ONLY_MODEL_CHECKS):
        raise ValueError(
            "audit-only model-sensitivity checks must exactly match the reviewed "
            f"set: expected {sorted(REVIEWED_AUDIT_ONLY_MODEL_CHECKS)!r}, "
            f"found {sorted(rows)!r}"
        )
    return rows


def load_nonadopted_ratio_fits(
    path: Path = DEFAULT_NONADOPTED_FITS,
    *,
    source_treatments: dict[str, str] | None = None,
) -> dict[str, dict[str, object]]:
    """Load published conditional fits that are explicitly not adopted."""

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        _validate_header(reader.fieldnames, NONADOPTED_FIT_FIELDS, path)
        raw_rows = list(reader)
    if any(None in row for row in raw_rows):
        raise ValueError(f"one or more rows in {path} contain extra columns")
    if not raw_rows:
        raise ValueError(f"at least one non-adopted ratio fit is required in {path}")

    rows: dict[str, dict[str, object]] = {}
    fit_keys: set[str] = set()
    for raw in raw_rows:
        source_key = raw["source_key"]
        fit_key = raw["fit_key"]
        if not source_key or not fit_key:
            raise ValueError("non-adopted source_key and fit_key are required")
        if source_key in rows:
            raise ValueError(f"duplicate non-adopted ratio-fit source: {source_key}")
        if fit_key in fit_keys:
            raise ValueError(f"duplicate non-adopted ratio fit_key: {fit_key}")

        reviewed = REVIEWED_NONADOPTED_FITS.get(source_key)
        if reviewed is None:
            if source_treatments is not None and source_key not in source_treatments:
                raise ValueError(f"unknown non-adopted ratio-fit source: {source_key}")
            raise ValueError(f"unsupported reviewed non-adopted source: {source_key}")
        primary_policy = NONADOPTED_FIT_PRIMARY_POLICIES[source_key]
        if source_treatments is not None:
            if source_key not in source_treatments:
                raise ValueError(f"unknown non-adopted ratio-fit source: {source_key}")
            expected_treatment = primary_policy["treatment"]
            if source_treatments[source_key] != expected_treatment:
                raise ValueError(
                    f"{source_key}: non-adopted ratio fit requires "
                    f"{expected_treatment} primary treatment"
                )

        data_relation = raw["data_relation"]
        model_basis = raw["model_basis"]
        interval_status = raw["interval_status"]
        adoption_status = raw["adoption_status"]
        nonadoption_basis = raw["nonadoption_basis"]
        if data_relation not in NONADOPTED_FIT_DATA_RELATIONS:
            raise ValueError(
                f"{source_key}: unsupported non-adopted data relation {data_relation!r}"
            )
        if model_basis not in NONADOPTED_FIT_MODEL_BASES:
            raise ValueError(
                f"{source_key}: unsupported non-adopted model basis {model_basis!r}"
            )
        if interval_status not in NONADOPTED_FIT_INTERVAL_STATUSES:
            raise ValueError(
                f"{source_key}: unsupported non-adopted interval status "
                f"{interval_status!r}"
            )
        if adoption_status not in NONADOPTED_FIT_ADOPTION_STATUSES:
            raise ValueError(
                f"{source_key}: unsupported non-adopted adoption status "
                f"{adoption_status!r}"
            )
        if nonadoption_basis not in NONADOPTED_FIT_BASES:
            raise ValueError(
                f"{source_key}: unsupported non-adoption basis {nonadoption_basis!r}"
            )

        ratio = optional_float(raw["ratio_msun_per_mpc"])
        ratio_minus = optional_float(raw["ratio_unc_minus_msun_per_mpc"])
        ratio_plus = optional_float(raw["ratio_unc_plus_msun_per_mpc"])
        if any(
            value is None or value <= 0.0
            for value in (ratio, ratio_minus, ratio_plus)
        ):
            raise ValueError(
                f"{source_key}: non-adopted ratio and uncertainties must be positive"
            )
        if ratio_minus >= ratio or ratio_plus >= ratio:
            raise ValueError(
                f"{source_key}: non-adopted ratio uncertainties must be smaller "
                "than the ratio"
            )

        reference_key = raw["reference_key"]
        reference_doi = raw["reference_doi"]
        if reference_key not in NONADOPTED_FIT_REFERENCE_DOIS:
            raise ValueError(f"{source_key}: unsupported non-adopted reference_key")
        if reference_doi != NONADOPTED_FIT_REFERENCE_DOIS[reference_key]:
            raise ValueError(
                f"{source_key}: non-adopted reference key does not match its DOI"
            )
        if not raw["provenance_note"]:
            raise ValueError(f"{source_key}: non-adopted provenance_note is required")

        reviewed_values = {
            "fit_key": fit_key,
            "data_relation": data_relation,
            "model_basis": model_basis,
            "ratio_msun_per_mpc": ratio,
            "ratio_unc_minus_msun_per_mpc": ratio_minus,
            "ratio_unc_plus_msun_per_mpc": ratio_plus,
            "interval_status": interval_status,
            "adoption_status": adoption_status,
            "nonadoption_basis": nonadoption_basis,
            "reference_key": reference_key,
            "reference_doi": reference_doi,
            "provenance_note": raw["provenance_note"],
        }
        for field, expected in reviewed.items():
            if reviewed_values[field] != expected:
                raise ValueError(
                    f"{source_key}: {field} does not match the reviewed "
                    "non-adopted fit"
                )

        rows[source_key] = {
            "source_key": source_key,
            "fit_key": fit_key,
            "data_relation": data_relation,
            "model_basis": model_basis,
            "ratio_msun_per_mpc": ratio,
            "ratio_unc_minus_msun_per_mpc": ratio_minus,
            "ratio_unc_plus_msun_per_mpc": ratio_plus,
            "size_uas": shadow_size_uas(float(ratio), 1.0),
            "size_unc_minus_uas": shadow_size_uas(float(ratio_minus), 1.0),
            "size_unc_plus_uas": shadow_size_uas(float(ratio_plus), 1.0),
            "interval_status": interval_status,
            "adoption_status": adoption_status,
            "nonadoption_basis": nonadoption_basis,
            "reference_key": reference_key,
            "reference_doi": reference_doi,
            "provenance_note": raw["provenance_note"],
            "primary_treatment": primary_policy["treatment"],
            "primary_result_status": primary_policy["primary_result_status"],
        }
        fit_keys.add(fit_key)

    if set(rows) != EXPECTED_NONADOPTED_FIT_SOURCES:
        raise ValueError(
            "non-adopted ratio fits must exactly match the reviewed source set: "
            f"expected {sorted(EXPECTED_NONADOPTED_FIT_SOURCES)!r}, "
            f"found {sorted(rows)!r}"
        )
    return rows


def log_uncertainties(
    mass: float | None,
    minus: float | None,
    plus: float | None,
) -> tuple[float | None, float | None, float | None, float | None]:
    if mass is None:
        return None, None, None, None
    log_mass = math.log10(mass)
    if minus is None and plus is None:
        return log_mass, None, None, None
    if minus is None or plus is None or minus >= mass:
        raise ValueError("invalid mass uncertainty")
    log_minus = math.log10(mass / (mass - minus))
    log_plus = math.log10((mass + plus) / mass)
    return log_mass, (log_minus + log_plus) / 2.0, log_minus, log_plus


def distance_components(row: dict[str, str]) -> dict[str, float | None]:
    measurement_minus = optional_float(row["distance_measurement_unc_minus_mpc"])
    measurement_plus = optional_float(row["distance_measurement_unc_plus_mpc"])
    flow = optional_float(row["distance_flow_unc_mpc"])
    systematic = optional_float(row["distance_shared_systematic_unc_mpc"])
    stat_minus = positive_quadrature(measurement_minus, flow)
    stat_plus = positive_quadrature(measurement_plus, flow)
    systematic = systematic if systematic is not None and systematic > 0.0 else None
    return {
        "stat_minus": stat_minus,
        "stat_plus": stat_plus,
        "systematic_minus": systematic,
        "systematic_plus": systematic,
        "total_minus": positive_quadrature(stat_minus, systematic),
        "total_plus": positive_quadrature(stat_plus, systematic),
    }


def mass_nondistance_fractions(
    row: dict[str, str],
) -> tuple[float | None, float | None]:
    """Return the retained non-distance interval, whatever its evidence type.

    Most rows supply a formal fit interval.  NGC 1194 and NGC 2273 instead
    carry the undecomposed mass uncertainties reported by Gao et al. (2017).
    Their compilation does not identify how its values differ from the smaller
    fixed-distance fit-plus-position errors in Kuo et al. (2011), so the
    neutral function/input names avoid assigning an unsupported error type.
    """

    override_minus = optional_float(row["size_stat_fraction_minus_override"])
    override_plus = optional_float(row["size_stat_fraction_plus_override"])
    if override_minus is not None or override_plus is not None:
        if override_minus is None or override_plus is None:
            raise ValueError(f"{row['source_key']}: both statistical overrides are required")
        return override_minus, override_plus

    mass = optional_float(row["reference_mass_msun"])
    error_minus = optional_float(row["reference_mass_unc_minus_msun"])
    error_plus = optional_float(row["reference_mass_unc_plus_msun"])
    if error_minus is None and error_plus is None:
        return None, None
    if mass is None or error_minus is None or error_plus is None:
        raise ValueError(f"{row['source_key']}: both mass-error bounds are required")
    return error_minus / mass, error_plus / mass


def _base_output(row: dict[str, str], distance: dict[str, float | None]) -> dict[str, object]:
    adopted_distance = float(row["adopted_distance_mpc"])
    distance_total_minus = distance["total_minus"]
    distance_total_plus = distance["total_plus"]
    return {
        "source_key": row["source_key"],
        "source_label": row["source_label"],
        "size_unc_independent_uas": None,
        "size_unc_correlated_uas": None,
        "distance_mpc": adopted_distance,
        "distance_unc_mpc": mean_optional(distance_total_minus, distance_total_plus),
        "distance_unc_minus_mpc": distance_total_minus,
        "distance_unc_plus_mpc": distance_total_plus,
        "distance_stat_unc_minus_mpc": distance["stat_minus"],
        "distance_stat_unc_plus_mpc": distance["stat_plus"],
        "distance_stat_component_status": row["distance_stat_component_status"],
        "distance_systematic_unc_minus_mpc": distance["systematic_minus"],
        "distance_systematic_unc_plus_mpc": distance["systematic_plus"],
        "parameter_reference_doi": row["mass_reference_doi"],
        "mass_uncertainty_reference_doi": (
            None
            if row["mass_uncertainty_reference_doi"] == "NA"
            else row["mass_uncertainty_reference_doi"]
        ),
        "distance_reference_doi": row["distance_reference_doi"],
        "treatment": row["treatment"],
        "source_order": int(row["source_order"]),
        "size_model_check_key": None,
        "size_model_check_type": None,
        "size_model_check_data_relation": None,
        "size_model_check_range_basis": None,
        "size_model_check_accepted_configuration_count": None,
        "size_model_check_excluded_configuration_ids": None,
        "size_model_check_lower_uas": None,
        "size_model_check_upper_uas": None,
        "size_model_check_probability_status": None,
        "size_model_check_adoption_status": None,
        "size_model_check_reference_key": None,
        "size_model_check_reference_doi": None,
        "nonadopted_fit_key": None,
        "nonadopted_fit_ratio_msun_per_mpc": None,
        "nonadopted_fit_ratio_unc_minus_msun_per_mpc": None,
        "nonadopted_fit_ratio_unc_plus_msun_per_mpc": None,
        "nonadopted_fit_size_uas": None,
        "nonadopted_fit_size_unc_minus_uas": None,
        "nonadopted_fit_size_unc_plus_uas": None,
        "nonadopted_fit_interval_status": None,
        "nonadopted_fit_adoption_status": None,
        "nonadopted_fit_basis": None,
        "nonadopted_fit_reference_key": None,
        "nonadopted_fit_reference_doi": None,
    }


def build_no_mass_row(row: dict[str, str]) -> dict[str, object]:
    distance = distance_components(row)
    output = _base_output(row, distance)
    output.update(
        {
            "size_uas": None,
            "size_unc_adopted_uas": None,
            "adopted_uncertainty_model": "not_available",
            "log_mass": None,
            "log_mass_unc": None,
            "mass_frac_percent": None,
            "mass_msun": None,
            "mass_unc_minus_msun": None,
            "mass_unc_plus_msun": None,
            "log_mass_unc_minus": None,
            "log_mass_unc_plus": None,
            "mass_systematic_unc_minus_msun": None,
            "mass_systematic_unc_plus_msun": None,
            "log_mass_systematic_unc_minus": None,
            "log_mass_systematic_unc_plus": None,
            "mass_systematic_status": "not_available",
            "uncertainty_status": "no robust primary black-hole mass",
            "size_posterior_q16_uas": None,
            "size_posterior_q50_uas": None,
            "size_posterior_q84_uas": None,
            "size_unc_minus_uas": None,
            "size_unc_plus_uas": None,
            "mass_reference_distance_mpc": None,
            "mass_fit_unc_minus_msun": None,
            "mass_fit_unc_plus_msun": None,
            "mass_shared_distance_unc_minus_msun": None,
            "mass_shared_distance_unc_plus_msun": None,
            "mass_total_unc_minus_msun": None,
            "mass_total_unc_plus_msun": None,
            "mass_display_unc_minus_msun": None,
            "mass_display_unc_plus_msun": None,
            "mass_display_interval_status": "not_available",
            "size_stat_unc_minus_uas": None,
            "size_stat_unc_plus_uas": None,
            "size_stat_component_status": "not_available",
            "size_systematic_unc_minus_uas": None,
            "size_systematic_unc_plus_uas": None,
            "size_total_unc_minus_uas": None,
            "size_total_unc_plus_uas": None,
            "size_systematic_status": "not_applicable",
            "size_total_status": "not_available",
            "legacy_uncertainty_basis": "not_available",
            "ratio_provenance": "no_robust_primary_mass",
        }
    )
    return output


def build_mass_record(
    row: dict[str, str], distance: dict[str, float | None]
) -> dict[str, object]:
    reference_mass = float(row["reference_mass_msun"])
    reference_distance = float(row["reference_distance_mpc"])
    adopted_distance = float(row["adopted_distance_mpc"])
    scale = adopted_distance / reference_distance if row["treatment"] == "shared_distance" else 1.0
    mass = reference_mass * scale
    fit_fraction_minus, fit_fraction_plus = mass_nondistance_fractions(row)
    fit_minus = mass * fit_fraction_minus if fit_fraction_minus is not None else None
    fit_plus = mass * fit_fraction_plus if fit_fraction_plus is not None else None
    mass_systematic_minus_input = optional_float(
        row["reference_mass_systematic_unc_minus_msun"]
    )
    mass_systematic_plus_input = optional_float(
        row["reference_mass_systematic_unc_plus_msun"]
    )
    mass_systematic_minus = (
        mass_systematic_minus_input * scale
        if mass_systematic_minus_input is not None
        else None
    )
    mass_systematic_plus = (
        mass_systematic_plus_input * scale
        if mass_systematic_plus_input is not None
        else None
    )

    shared_minus = None
    shared_plus = None
    if row["treatment"] == "shared_distance":
        if distance["total_minus"] is not None:
            shared_minus = mass * distance["total_minus"] / adopted_distance
        if distance["total_plus"] is not None:
            shared_plus = mass * distance["total_plus"] / adopted_distance

    total_minus = positive_quadrature(
        fit_minus, shared_minus, mass_systematic_minus
    )
    total_plus = positive_quadrature(fit_plus, shared_plus, mass_systematic_plus)
    # ``log_unc`` is a backward-compatible scalar derived only from the
    # fit/statistical mass interval.  A source-reported model-selection
    # systematic is retained in its own exact machine-readable log bounds and
    # never folded into the paired M/D shadow-size posterior.
    log_mass, log_unc, log_minus, log_plus = log_uncertainties(mass, fit_minus, fit_plus)
    (
        _systematic_log_mass,
        _systematic_log_unc,
        log_systematic_minus,
        log_systematic_plus,
    ) = log_uncertainties(mass, mass_systematic_minus, mass_systematic_plus)
    return {
        "mass": mass,
        "fit_fraction_minus": fit_fraction_minus,
        "fit_fraction_plus": fit_fraction_plus,
        "fit_minus": fit_minus,
        "fit_plus": fit_plus,
        "systematic_minus": mass_systematic_minus,
        "systematic_plus": mass_systematic_plus,
        "log_systematic_minus": log_systematic_minus,
        "log_systematic_plus": log_systematic_plus,
        "systematic_status": row["mass_systematic_status"],
        "shared_minus": shared_minus,
        "shared_plus": shared_plus,
        "total_minus": total_minus,
        "total_plus": total_plus,
        "log_mass": log_mass,
        "log_unc": log_unc,
        "log_minus": log_minus,
        "log_plus": log_plus,
        "reference_distance": reference_distance,
    }


def mass_display_interval(
    row: dict[str, str], mass: dict[str, float | None]
) -> tuple[float | None, float | None, str]:
    """Select an explicit publication mass interval without changing diagnostics.

    Most rows display the available decomposed interval.  A
    ``separate_decomposed_and_position_systematic_components`` row displays
    the quadrature of its formal fit and adopted-distance terms while
    retaining a source-reported absolute-position systematic separately.  A
    ``separate_decomposed_and_geometry_systematic_components`` row uses the
    same base interval while retaining a source-reported disk-geometry
    systematic separately.  A
    ``separate_statistical_and_systematic_components`` row displays the
    fit/statistical interval while retaining its source-reported mass
    systematic as a distinct component.  A
    ``source_reported_total`` row instead retains the uncertainty printed by the
    primary source, rescaled with the mass normalization when necessary.  The
    decomposed fit, shared-distance, and quadrature-total fields remain intact
    for auditing and are never silently substituted into this presentation
    interval.
    """

    status = row["mass_display_interval_status"]
    if status == "decomposed_components":
        lower = mass["total_minus"]
        upper = mass["total_plus"]
        if lower is None or upper is None:
            raise ValueError(
                f"{row['source_key']}: decomposed mass display interval is incomplete"
            )
        return lower, upper, status
    if status == "separate_decomposed_and_position_systematic_components":
        lower = positive_quadrature(mass["fit_minus"], mass["shared_minus"])
        upper = positive_quadrature(mass["fit_plus"], mass["shared_plus"])
        if lower is None or upper is None:
            raise ValueError(
                f"{row['source_key']}: decomposed mass base interval is incomplete"
            )
        if mass["systematic_minus"] is None or mass["systematic_plus"] is None:
            raise ValueError(
                f"{row['source_key']}: separate position systematic is incomplete"
            )
        return lower, upper, status
    if status == "separate_decomposed_and_geometry_systematic_components":
        lower = positive_quadrature(mass["fit_minus"], mass["shared_minus"])
        upper = positive_quadrature(mass["fit_plus"], mass["shared_plus"])
        if lower is None or upper is None:
            raise ValueError(
                f"{row['source_key']}: decomposed mass base interval is incomplete"
            )
        if mass["systematic_minus"] is None or mass["systematic_plus"] is None:
            raise ValueError(
                f"{row['source_key']}: separate geometry systematic is incomplete"
            )
        return lower, upper, status
    if status == "separate_statistical_and_systematic_components":
        lower = mass["fit_minus"]
        upper = mass["fit_plus"]
        if lower is None or upper is None:
            raise ValueError(
                f"{row['source_key']}: separate statistical mass interval is incomplete"
            )
        if mass["systematic_minus"] is None or mass["systematic_plus"] is None:
            raise ValueError(
                f"{row['source_key']}: separate mass systematic is incomplete"
            )
        return lower, upper, status
    if status == "source_reported_total":
        reference_mass = float(row["reference_mass_msun"])
        reported_minus = optional_float(row["reference_mass_unc_minus_msun"])
        reported_plus = optional_float(row["reference_mass_unc_plus_msun"])
        if reported_minus is None or reported_plus is None:
            raise ValueError(
                f"{row['source_key']}: source-reported mass interval is incomplete"
            )
        scale = float(mass["mass"]) / reference_mass
        return reported_minus * scale, reported_plus * scale, status
    if status == "not_available":
        return None, None, status
    raise ValueError(
        f"{row['source_key']}: unsupported mass display interval status {status!r}"
    )


def build_paired_row(
    row: dict[str, str], paired: dict[str, object]
) -> dict[str, object]:
    distance = distance_components(row)
    mass = build_mass_record(row, distance)
    mass_display_minus, mass_display_plus, mass_display_status = (
        mass_display_interval(row, mass)
    )
    q16 = theta_g_nas_to_shadow_uas(float(paired["q16_theta_g_nas"]))
    q50 = theta_g_nas_to_shadow_uas(float(paired["q50_theta_g_nas"]))
    q84 = theta_g_nas_to_shadow_uas(float(paired["q84_theta_g_nas"]))
    stat_minus = q50 - q16
    stat_plus = q84 - q50
    systematic_minus = q50 * float(paired["systematic_fraction"])
    systematic_plus = systematic_minus
    total_minus = math.hypot(stat_minus, systematic_minus)
    total_plus = math.hypot(stat_plus, systematic_plus)
    adopted_uncertainty = (total_minus + total_plus) / 2.0

    output = _base_output(row, distance)
    output.update(
        {
            "size_uas": q50,
            "size_unc_adopted_uas": adopted_uncertainty,
            "adopted_uncertainty_model": "paired_m_over_d_posterior",
            "log_mass": mass["log_mass"],
            "log_mass_unc": mass["log_unc"],
            "mass_frac_percent": mean_optional(
                mass["fit_fraction_minus"], mass["fit_fraction_plus"]
            )
            * 100.0,
            "mass_msun": mass["mass"],
            "mass_unc_minus_msun": mass["fit_minus"],
            "mass_unc_plus_msun": mass["fit_plus"],
            "log_mass_unc_minus": mass["log_minus"],
            "log_mass_unc_plus": mass["log_plus"],
            "mass_systematic_unc_minus_msun": mass["systematic_minus"],
            "mass_systematic_unc_plus_msun": mass["systematic_plus"],
            "log_mass_systematic_unc_minus": mass["log_systematic_minus"],
            "log_mass_systematic_unc_plus": mass["log_systematic_plus"],
            "mass_systematic_status": mass["systematic_status"],
            "uncertainty_status": (
                "paired M/D posterior; asymmetric statistical interval and "
                "separate 1% non-Keplerian allowance"
            ),
            "size_posterior_q16_uas": q16,
            "size_posterior_q50_uas": q50,
            "size_posterior_q84_uas": q84,
            "size_unc_minus_uas": stat_minus,
            "size_unc_plus_uas": stat_plus,
            "mass_reference_distance_mpc": mass["reference_distance"],
            "mass_fit_unc_minus_msun": mass["fit_minus"],
            "mass_fit_unc_plus_msun": mass["fit_plus"],
            "mass_shared_distance_unc_minus_msun": None,
            "mass_shared_distance_unc_plus_msun": None,
            "mass_total_unc_minus_msun": mass["total_minus"],
            "mass_total_unc_plus_msun": mass["total_plus"],
            "mass_display_unc_minus_msun": mass_display_minus,
            "mass_display_unc_plus_msun": mass_display_plus,
            "mass_display_interval_status": mass_display_status,
            "size_stat_unc_minus_uas": stat_minus,
            "size_stat_unc_plus_uas": stat_plus,
            "size_stat_component_status": "paired_chain_quantile_interval",
            "size_systematic_unc_minus_uas": systematic_minus,
            "size_systematic_unc_plus_uas": systematic_plus,
            "size_total_unc_minus_uas": total_minus,
            "size_total_unc_plus_uas": total_plus,
            "size_systematic_status": "numeric_one_percent_non_keplerian_allowance",
            "size_total_status": "quadrature_statistical_plus_systematic",
            "legacy_uncertainty_basis": "symmetrized_quadrature_total",
            "ratio_provenance": paired["provenance"],
        }
    )
    return output


def build_shared_row(row: dict[str, str]) -> dict[str, object]:
    distance = distance_components(row)
    mass = build_mass_record(row, distance)
    mass_display_minus, mass_display_plus, mass_display_status = (
        mass_display_interval(row, mass)
    )
    size = shadow_size_uas(float(mass["mass"]), float(row["adopted_distance_mpc"]))
    reference_size = shadow_size_uas(
        float(row["reference_mass_msun"]), float(row["reference_distance_mpc"])
    )
    if not math.isclose(size, reference_size, rel_tol=2e-15, abs_tol=0.0):
        raise AssertionError(f"{row['source_key']}: shared distance failed to cancel")

    stat_minus = (
        size * mass["fit_fraction_minus"]
        if mass["fit_fraction_minus"] is not None
        else None
    )
    stat_plus = (
        size * mass["fit_fraction_plus"]
        if mass["fit_fraction_plus"] is not None
        else None
    )
    systematic_fraction_minus = optional_float(row["size_systematic_fraction_minus"])
    systematic_fraction_plus = optional_float(row["size_systematic_fraction_plus"])
    systematic_minus = (
        size * systematic_fraction_minus if systematic_fraction_minus is not None else None
    )
    systematic_plus = (
        size * systematic_fraction_plus if systematic_fraction_plus is not None else None
    )

    status = row["systematic_status"]
    if status in {"numeric_systematic", "model_sensitivity_envelope"}:
        if systematic_minus is None or systematic_plus is None:
            raise ValueError(f"{row['source_key']}: numerical systematic is missing")
        total_minus = positive_quadrature(stat_minus, systematic_minus)
        total_plus = positive_quadrature(stat_plus, systematic_plus)
        total_status = (
            "quadrature_for_display_model_envelope_not_probability_interval"
            if status == "model_sensitivity_envelope"
            else (
                "quadrature_mixed_fit_position_plus_systematic"
                if row["source_key"] == "ESO558-G009"
                else "quadrature_statistical_plus_systematic"
            )
        )
    elif status in {"none_adopted", "reported_undecomposed_uncertainty"}:
        total_minus, total_plus = stat_minus, stat_plus
        total_status = (
            "reported_undecomposed_uncertainty_component_only"
            if status == "reported_undecomposed_uncertainty"
            else "formal_component_only"
        )
    elif status in {
        "unquantified_astrometric",
        "unquantified_model",
        "pending_nakai",
    }:
        total_minus, total_plus = None, None
        total_status = "incomplete_unquantified_systematic"
    elif status == "no_formal_uncertainty":
        total_minus, total_plus = None, None
        total_status = "no_formal_uncertainty"
    else:
        raise ValueError(f"{row['source_key']}: invalid shared-distance status {status!r}")

    if (
        status == "reported_undecomposed_uncertainty"
        and total_minus is not None
        and total_plus is not None
    ):
        adopted_uncertainty = (total_minus + total_plus) / 2.0
        legacy_basis = "reported_undecomposed_uncertainty_scalar"
    elif total_minus is not None and total_plus is not None:
        adopted_uncertainty = (total_minus + total_plus) / 2.0
        legacy_basis = "symmetrized_available_total"
    elif stat_minus is not None and stat_plus is not None:
        adopted_uncertainty = (stat_minus + stat_plus) / 2.0
        if row["source_key"] == "Circinus":
            legacy_basis = (
                "symmetrized_angular_radius_only_component_total_incomplete"
            )
        elif row["source_key"] in FIT_PLUS_ABSOLUTE_POSITION_SOURCES:
            legacy_basis = (
                "symmetrized_fit_plus_absolute_position_component_total_incomplete"
            )
        else:
            legacy_basis = "symmetrized_formal_statistical_term_total_incomplete"
    else:
        adopted_uncertainty = None
        legacy_basis = "not_available"

    if status == "pending_nakai":
        uncertainty_status = (
            "shared distance canceled; 0.17/4.14 formal fit term; additional "
            "systematic pending Nakai"
        )
    elif status == "unquantified_astrometric":
        uncertainty_status = (
            "shared distance canceled; numerical formal term only; residual "
            "absolute-coordinate and error-floor completeness unreported"
        )
    elif status == "unquantified_model":
        if row["source_key"] == "Circinus":
            uncertainty_status = (
                "shared distance canceled; angular-radius-only component retained; "
                "model systematic unquantified"
            )
        elif stat_minus is None and stat_plus is None:
            uncertainty_status = (
                "shared distance canceled; approximate enclosed-mass black-hole "
                "proxy has no formal fit interval; model systematic unquantified"
            )
        elif row["source_key"] in FIT_PLUS_ABSOLUTE_POSITION_SOURCES:
            uncertainty_status = (
                "shared distance canceled; fixed-distance fitting plus "
                "absolute-position component retained; model systematic unquantified"
            )
        else:
            uncertainty_status = (
                "shared distance canceled; numerical formal term only; "
                "model systematic unquantified"
            )
    elif status == "model_sensitivity_envelope":
        uncertainty_status = "shared distance canceled; formal interval and separate model-sensitivity envelope"
    elif status == "numeric_systematic":
        uncertainty_status = (
            "shared distance canceled; mixed fit-plus-absolute-position component "
            "and systematic component retained separately"
            if row["source_key"] == "ESO558-G009"
            else "shared distance canceled; statistical and systematic components retained separately"
        )
    elif status == "no_formal_uncertainty":
        uncertainty_status = "shared distance canceled; no formal mass-fit uncertainty published"
    elif status == "reported_undecomposed_uncertainty":
        uncertainty_status = (
            "shared distance canceled; reported undecomposed mass uncertainty propagated"
        )
    else:
        uncertainty_status = "shared distance canceled; fixed-distance fit uncertainty only"

    output = _base_output(row, distance)
    output.update(
        {
            "size_uas": size,
            "size_unc_adopted_uas": adopted_uncertainty,
            "adopted_uncertainty_model": (
                "shared_distance_cancellation_no_formal_uncertainty"
                if status == "no_formal_uncertainty"
                else "shared_distance_cancellation"
            ),
            "log_mass": mass["log_mass"],
            "log_mass_unc": mass["log_unc"],
            "mass_frac_percent": (
                mean_optional(mass["fit_fraction_minus"], mass["fit_fraction_plus"])
                * 100.0
                if mass["fit_fraction_minus"] is not None
                else None
            ),
            "mass_msun": mass["mass"],
            "mass_unc_minus_msun": mass["fit_minus"],
            "mass_unc_plus_msun": mass["fit_plus"],
            "log_mass_unc_minus": mass["log_minus"],
            "log_mass_unc_plus": mass["log_plus"],
            "mass_systematic_unc_minus_msun": mass["systematic_minus"],
            "mass_systematic_unc_plus_msun": mass["systematic_plus"],
            "log_mass_systematic_unc_minus": mass["log_systematic_minus"],
            "log_mass_systematic_unc_plus": mass["log_systematic_plus"],
            "mass_systematic_status": mass["systematic_status"],
            "uncertainty_status": uncertainty_status,
            "size_posterior_q16_uas": None,
            "size_posterior_q50_uas": None,
            "size_posterior_q84_uas": None,
            "size_unc_minus_uas": stat_minus,
            "size_unc_plus_uas": stat_plus,
            "mass_reference_distance_mpc": mass["reference_distance"],
            "mass_fit_unc_minus_msun": mass["fit_minus"],
            "mass_fit_unc_plus_msun": mass["fit_plus"],
            "mass_shared_distance_unc_minus_msun": mass["shared_minus"],
            "mass_shared_distance_unc_plus_msun": mass["shared_plus"],
            "mass_total_unc_minus_msun": (
                None
                if mass_display_status == "not_available"
                else mass["total_minus"]
            ),
            "mass_total_unc_plus_msun": (
                None
                if mass_display_status == "not_available"
                else mass["total_plus"]
            ),
            "mass_display_unc_minus_msun": mass_display_minus,
            "mass_display_unc_plus_msun": mass_display_plus,
            "mass_display_interval_status": mass_display_status,
            "size_stat_unc_minus_uas": stat_minus,
            "size_stat_unc_plus_uas": stat_plus,
            "size_stat_component_status": (
                "no_formal_uncertainty"
                if stat_minus is None and stat_plus is None
                else (
                    "fit_plus_absolute_position_mixed_component"
                    if row["source_key"] in FIT_PLUS_ABSOLUTE_POSITION_SOURCES
                    else (
                        "reported_undecomposed_uncertainty"
                        if status == "reported_undecomposed_uncertainty"
                        else (
                            "angular_radius_only_component"
                            if row["source_key"] == "Circinus"
                            else "fixed_distance_formal_fit_interval"
                        )
                    )
                )
            ),
            "size_systematic_unc_minus_uas": systematic_minus,
            "size_systematic_unc_plus_uas": systematic_plus,
            "size_total_unc_minus_uas": total_minus,
            "size_total_unc_plus_uas": total_plus,
            "size_systematic_status": status,
            "size_total_status": total_status,
            "legacy_uncertainty_basis": legacy_basis,
            "ratio_provenance": "published_mass_normalization_over_its_reference_distance",
        }
    )
    return output


def attach_model_check(
    output: dict[str, object], check: dict[str, object]
) -> dict[str, object]:
    """Attach a separately reported check without changing adopted uncertainties."""

    source_key = str(output["source_key"])
    if check["source_key"] != source_key:
        raise ValueError(
            f"{source_key}: model-check source does not match output source"
        )
    primary_policy = MODEL_CHECK_PRIMARY_POLICIES.get(source_key)
    if primary_policy is None:
        raise ValueError(f"{source_key}: unsupported model-check source")
    if output["treatment"] != primary_policy["treatment"]:
        raise ValueError(f"{source_key}: model-check primary treatment mismatch")
    if output["size_total_status"] != primary_policy["size_total_status"]:
        raise ValueError(f"{source_key}: model-check primary uncertainty policy mismatch")
    primary_size = output["size_uas"]
    primary_result_status = primary_policy["primary_result_status"]
    if primary_result_status == "retained":
        if (
            output["mass_msun"] is None
            or primary_size is None
            or float(output["mass_msun"]) <= 0.0
            or float(primary_size) <= 0.0
        ):
            raise ValueError(
                f"{source_key}: model check requires the retained primary mass and size"
            )
    elif primary_result_status == "absent":
        if output["mass_msun"] is not None or primary_size is not None:
            raise ValueError(
                f"{source_key}: standalone diagnostic cannot coexist with an adopted primary"
            )
        if check["sensitivity_direction"] != "two_sided_range":
            raise ValueError(
                f"{source_key}: standalone diagnostic must remain a two-sided reported span"
            )
        if check["endpoint_derivation"] != (
            "reported_mbh_plusminus_error_at_fixed_reported_distance"
        ):
            raise ValueError(
                f"{source_key}: standalone diagnostic requires absolute reported endpoints"
            )
    else:
        raise ValueError(f"{source_key}: unsupported primary-result policy")
    if output["size_model_check_key"] is not None:
        raise ValueError(f"{source_key}: model-sensitivity check is already attached")

    if check["endpoint_derivation"] in {
        "published_mass_distance_endpoints",
        "reported_mbh_plusminus_error_at_fixed_reported_distance",
    }:
        lower_size = check["size_lower_uas"]
        upper_size = check["size_upper_uas"]
    elif check["endpoint_derivation"] == "adopted_primary_ratio_factors":
        if primary_result_status != "retained":
            raise ValueError(
                f"{source_key}: a no-primary diagnostic cannot use primary-ratio factors"
            )
        lower_size = float(primary_size) * float(check["lower_primary_ratio_factor"])
        upper_size = float(primary_size) * float(check["upper_primary_ratio_factor"])
    else:
        raise ValueError(f"{source_key}: unsupported model-check endpoint derivation")
    if (
        lower_size is None
        or upper_size is None
        or float(lower_size) <= 0.0
        or float(lower_size) >= float(upper_size)
    ):
        raise ValueError(f"{source_key}: model-check size bounds must be ordered")

    output.update(
        {
            "size_model_check_key": check["check_key"],
            "size_model_check_type": check["check_type"],
            "size_model_check_data_relation": check["data_relation"],
            "size_model_check_range_basis": check["range_basis"],
            "size_model_check_accepted_configuration_count": check[
                "accepted_configuration_count"
            ],
            "size_model_check_excluded_configuration_ids": check[
                "excluded_configuration_ids"
            ],
            "size_model_check_lower_uas": lower_size,
            "size_model_check_upper_uas": upper_size,
            "size_model_check_probability_status": check["probability_status"],
            "size_model_check_adoption_status": check["adoption_status"],
            "size_model_check_reference_key": check["reference_key"],
            "size_model_check_reference_doi": check["reference_doi"],
        }
    )
    if primary_result_status == "absent":
        sensitivity_label = "standalone disk-model diagnostic"
    else:
        sensitivity_labels = {
            "data_selection_component_tracking": "data-selection sensitivity",
            "disk_mass_model_choice": "disk-model sensitivity",
            "turbulent_support_mass_bias": "one-sided model sensitivity",
        }
        sensitivity_label = sensitivity_labels.get(str(check["check_type"]))
        if sensitivity_label is None:
            raise ValueError(f"{source_key}: unsupported model-check type")
    output["uncertainty_status"] = (
        f"{output['uncertainty_status']}; separate non-probabilistic "
        f"{sensitivity_label} not combined"
    )
    return output


def attach_nonadopted_ratio_fit(
    output: dict[str, object], fit: dict[str, object]
) -> dict[str, object]:
    """Attach a conditional literature fit without replacing the primary result."""

    source_key = str(output["source_key"])
    if fit["source_key"] != source_key:
        raise ValueError(
            f"{source_key}: non-adopted ratio-fit source does not match output source"
        )
    primary_policy = NONADOPTED_FIT_PRIMARY_POLICIES.get(source_key)
    if primary_policy is None:
        raise ValueError(f"{source_key}: unsupported non-adopted ratio-fit source")
    if (
        fit["primary_treatment"] != primary_policy["treatment"]
        or fit["primary_result_status"] != primary_policy["primary_result_status"]
    ):
        raise ValueError(f"{source_key}: non-adopted ratio-fit policy does not match")

    expected_treatment = primary_policy["treatment"]
    if output["treatment"] != expected_treatment:
        raise ValueError(
            f"{source_key}: non-adopted ratio fit requires "
            f"{expected_treatment} primary treatment"
        )
    primary_result_status = primary_policy["primary_result_status"]
    if primary_result_status == "absent" and (
        output["mass_msun"] is not None or output["size_uas"] is not None
    ):
        raise ValueError(
            f"{source_key}: non-adopted ratio fit cannot coexist with "
            "an adopted mass or shadow size"
        )
    if primary_result_status == "retained":
        mass = output["mass_msun"]
        size = output["size_uas"]
        if mass is None or size is None or float(mass) <= 0.0 or float(size) <= 0.0:
            raise ValueError(
                f"{source_key}: non-adopted ratio fit requires the retained "
                "primary mass and shadow size"
            )
    elif primary_result_status != "absent":
        raise ValueError(f"{source_key}: unsupported primary-result policy")
    if output["nonadopted_fit_key"] is not None:
        raise ValueError(f"{source_key}: non-adopted ratio fit is already attached")

    output.update(
        {
            "nonadopted_fit_key": fit["fit_key"],
            "nonadopted_fit_ratio_msun_per_mpc": fit["ratio_msun_per_mpc"],
            "nonadopted_fit_ratio_unc_minus_msun_per_mpc": fit[
                "ratio_unc_minus_msun_per_mpc"
            ],
            "nonadopted_fit_ratio_unc_plus_msun_per_mpc": fit[
                "ratio_unc_plus_msun_per_mpc"
            ],
            "nonadopted_fit_size_uas": fit["size_uas"],
            "nonadopted_fit_size_unc_minus_uas": fit["size_unc_minus_uas"],
            "nonadopted_fit_size_unc_plus_uas": fit["size_unc_plus_uas"],
            "nonadopted_fit_interval_status": fit["interval_status"],
            "nonadopted_fit_adoption_status": fit["adoption_status"],
            "nonadopted_fit_basis": fit["nonadoption_basis"],
            "nonadopted_fit_reference_key": fit["reference_key"],
            "nonadopted_fit_reference_doi": fit["reference_doi"],
        }
    )
    if primary_result_status == "absent":
        output["uncertainty_status"] = (
            f"{output['uncertainty_status']}; conditional model-dependent literature "
            "ratio reported separately and not adopted"
        )
    return output


def build_rows(
    input_path: Path = DEFAULT_INPUT,
    quantile_path: Path = DEFAULT_QUANTILES,
    model_check_path: Path = DEFAULT_MODEL_CHECKS,
    nonadopted_fit_path: Path = DEFAULT_NONADOPTED_FITS,
    audit_only_model_check_path: Path = DEFAULT_AUDIT_ONLY_MODEL_CHECKS,
) -> list[dict[str, object]]:
    source_rows = load_source_inputs(input_path)
    paired_quantiles = load_paired_quantiles(quantile_path)
    source_treatments = {
        row["source_key"]: row["treatment"] for row in source_rows
    }
    model_checks = load_model_checks(
        model_check_path,
        known_sources=set(source_treatments),
    )
    load_audit_only_model_checks(
        audit_only_model_check_path,
        known_sources=set(source_treatments),
    )
    nonadopted_fits = load_nonadopted_ratio_fits(
        nonadopted_fit_path,
        source_treatments=source_treatments,
    )
    if set(model_checks) != EXPECTED_MODEL_CHECK_SOURCES:
        raise ValueError(
            "model-sensitivity checks must exactly match the reviewed source set: "
            f"expected {sorted(EXPECTED_MODEL_CHECK_SOURCES)!r}, "
            f"found {sorted(model_checks)!r}"
        )
    expected_paired = {
        row["source_key"] for row in source_rows if row["treatment"] == "paired_posterior"
    }
    if set(paired_quantiles) != expected_paired:
        raise ValueError(
            "paired quantile keys must exactly match paired source inputs: "
            f"expected {sorted(expected_paired)!r}, found {sorted(paired_quantiles)!r}"
        )

    outputs: list[dict[str, object]] = []
    for row in source_rows:
        treatment = row["treatment"]
        if treatment == "paired_posterior":
            output = build_paired_row(row, paired_quantiles[row["source_key"]])
        elif treatment == "shared_distance":
            output = build_shared_row(row)
        elif treatment == "no_robust_mass":
            output = build_no_mass_row(row)
        else:  # Defensive: load_source_inputs validates the treatment counts.
            raise AssertionError(f"unhandled treatment: {treatment}")
        check = model_checks.get(row["source_key"])
        if check is not None:
            output = attach_model_check(output, check)
        nonadopted_fit = nonadopted_fits.get(row["source_key"])
        if nonadopted_fit is not None:
            output = attach_nonadopted_ratio_fit(output, nonadopted_fit)
        outputs.append(output)
    return outputs


def format_value(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("cannot serialize non-finite output")
        return format(value, ".17g")
    return str(value)


def render_csv(rows: list[dict[str, object]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=OUTPUT_FIELDS, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        missing = set(OUTPUT_FIELDS) - set(row)
        extra = set(row) - set(OUTPUT_FIELDS)
        if missing or extra:
            raise ValueError(f"output schema mismatch: missing={missing}, extra={extra}")
        if row["size_stat_component_status"] not in SIZE_STAT_COMPONENT_STATUSES:
            raise ValueError(
                f"{row['source_key']}: unsupported size statistical component status"
            )
        if row["mass_display_interval_status"] not in MASS_DISPLAY_INTERVAL_STATUSES:
            raise ValueError(
                f"{row['source_key']}: unsupported mass display interval status"
            )
        if row["mass_systematic_status"] not in MASS_SYSTEMATIC_STATUSES:
            raise ValueError(
                f"{row['source_key']}: unsupported mass systematic status"
            )
        if row["distance_stat_component_status"] not in DISTANCE_STAT_COMPONENT_STATUSES:
            raise ValueError(
                f"{row['source_key']}: unsupported distance statistical component status"
            )
        if row["source_key"] == "NGC4388":
            if (
                row["size_stat_component_status"]
                != "fit_plus_absolute_position_mixed_component"
            ):
                raise ValueError(
                    "NGC4388: reviewed size component must remain fit plus position"
                )
            if (
                row["legacy_uncertainty_basis"]
                != "symmetrized_fit_plus_absolute_position_component_total_incomplete"
            ):
                raise ValueError(
                    "NGC4388: reviewed incomplete-component basis does not match"
                )
            if row["distance_stat_component_status"] != "not_available":
                raise ValueError(
                    "NGC4388: Tully--Fisher distance term must remain systematic"
                )
            if (
                row["distance_stat_unc_minus_mpc"] is not None
                or row["distance_stat_unc_plus_mpc"] is not None
                or row["distance_systematic_unc_minus_mpc"] != 2.09
                or row["distance_systematic_unc_plus_mpc"] != 2.09
            ):
                raise ValueError(
                    "NGC4388: reviewed distance-component decomposition does not match"
                )
        if row["source_key"] == "IC2560":
            reviewed_ic2560 = {
                "size_uas": 0.009543265887703845,
                "size_stat_unc_minus_uas": 0.00039187323693469896,
                "size_stat_unc_plus_uas": 0.00039187323693469896,
                "mass_msun": 4_140_000.0,
                "mass_fit_unc_minus_msun": 170_000.0,
                "mass_fit_unc_plus_msun": 170_000.0,
                "mass_shared_distance_unc_minus_msun": 567_505.6179775281,
                "mass_shared_distance_unc_plus_msun": 567_505.6179775281,
                "mass_total_unc_minus_msun": 592_420.9875047103,
                "mass_total_unc_plus_msun": 592_420.9875047103,
                "mass_display_unc_minus_msun": 592_420.9875047103,
                "mass_display_unc_plus_msun": 592_420.9875047103,
                "distance_stat_unc_minus_mpc": 6.1,
                "distance_stat_unc_plus_mpc": 6.1,
                "distance_unc_minus_mpc": 6.1,
                "distance_unc_plus_mpc": 6.1,
            }
            for field, expected in reviewed_ic2560.items():
                actual = row[field]
                if actual is None or not math.isclose(
                    float(actual), expected, rel_tol=1e-15, abs_tol=1e-15
                ):
                    raise ValueError(
                        f"IC2560: {field} does not match the reviewed value"
                    )
            reviewed_ic2560_text = {
                "treatment": "shared_distance",
                "mass_systematic_status": "not_available",
                "mass_display_interval_status": "decomposed_components",
                "distance_stat_component_status": "statistical",
                "size_stat_component_status": "fixed_distance_formal_fit_interval",
                "size_systematic_status": "unquantified_model",
                "size_total_status": "incomplete_unquantified_systematic",
                "legacy_uncertainty_basis": (
                    "symmetrized_formal_statistical_term_total_incomplete"
                ),
            }
            for field, expected in reviewed_ic2560_text.items():
                if row[field] != expected:
                    raise ValueError(
                        f"IC2560: {field} does not match the reviewed value"
                    )
            for field in (
                "size_systematic_unc_minus_uas",
                "size_systematic_unc_plus_uas",
                "size_total_unc_minus_uas",
                "size_total_unc_plus_uas",
                "mass_systematic_unc_minus_msun",
                "mass_systematic_unc_plus_msun",
                "size_model_check_key",
                "nonadopted_fit_key",
            ):
                if row[field] is not None:
                    raise ValueError(f"IC2560: {field} must remain blank")
        if row["source_key"] == "J0437":
            reviewed_j0437 = {
                "size_uas": 0.004191067803754549,
                "size_stat_unc_minus_uas": 0.0002930816645982202,
                "size_stat_unc_plus_uas": 0.0002930816645982202,
                "mass_msun": 2_635_479.2190218437,
                "mass_fit_unc_minus_msun": 184_299.24608544362,
                "mass_fit_unc_plus_msun": 184_299.24608544362,
                "mass_shared_distance_unc_minus_msun": 176_180.85157776746,
                "mass_shared_distance_unc_plus_msun": 176_180.85157776746,
                "mass_total_unc_minus_msun": 254_962.55523180307,
                "mass_total_unc_plus_msun": 254_962.55523180307,
                "mass_display_unc_minus_msun": 254_962.55523180307,
                "mass_display_unc_plus_msun": 254_962.55523180307,
                "distance_stat_unc_minus_mpc": 3.4259780517380682,
                "distance_stat_unc_plus_mpc": 3.4259780517380682,
                "distance_systematic_unc_minus_mpc": 2.618595512716046,
                "distance_systematic_unc_plus_mpc": 2.618595512716046,
                "distance_unc_minus_mpc": 4.312118744910393,
                "distance_unc_plus_mpc": 4.312118744910393,
                "size_model_check_accepted_configuration_count": 1.0,
                "size_model_check_lower_uas": 0.002434118995564536,
                "size_model_check_upper_uas": 0.0031252885868976762,
                "nonadopted_fit_ratio_msun_per_mpc": 42_500.0,
                "nonadopted_fit_ratio_unc_minus_msun_per_mpc": 1_200.0,
                "nonadopted_fit_ratio_unc_plus_msun_per_mpc": 1_000.0,
                "nonadopted_fit_size_uas": 0.004359589760898525,
                "nonadopted_fit_size_unc_minus_uas": 0.00012309429913125248,
                "nonadopted_fit_size_unc_plus_uas": 0.00010257858260937706,
            }
            for field, expected in reviewed_j0437.items():
                actual = row[field]
                if actual is None or not math.isclose(
                    float(actual), expected, rel_tol=0.0, abs_tol=1e-15
                ):
                    raise ValueError(
                        f"J0437: {field} does not match the reviewed value"
                    )
            reviewed_j0437_text = {
                "treatment": "shared_distance",
                "mass_systematic_status": "not_available",
                "mass_display_interval_status": "decomposed_components",
                "distance_stat_component_status": "statistical",
                "size_stat_component_status": "fixed_distance_formal_fit_interval",
                "size_systematic_status": "unquantified_model",
                "size_total_status": "incomplete_unquantified_systematic",
                "legacy_uncertainty_basis": (
                    "symmetrized_formal_statistical_term_total_incomplete"
                ),
                "size_model_check_key": (
                    "deepshikha_2024_j0437_mestel_disk_mass_model"
                ),
                "size_model_check_type": "disk_mass_model_choice",
                "size_model_check_data_relation": (
                    "reanalysis_of_published_gao2017_vlbi_data_data_selection_"
                    "unreported"
                ),
                "size_model_check_range_basis": (
                    "reported_symmetric_mass_error_span"
                ),
                "size_model_check_probability_status": "non_probabilistic",
                "size_model_check_adoption_status": (
                    "reported_separately_not_adopted"
                ),
                "size_model_check_reference_key": "Deepshikha2024",
                "size_model_check_reference_doi": (
                    "10.1017/S1743921323002491"
                ),
                "nonadopted_fit_key": "villaraos_2022_j0437_gr_ratio_fit",
                "nonadopted_fit_interval_status": (
                    "bayesian_formal_conditional_on_model"
                ),
                "nonadopted_fit_adoption_status": (
                    "reported_separately_not_adopted"
                ),
                "nonadopted_fit_basis": (
                    "same_data_reanalysis_conditional_on_fixed_gao2017_"
                    "inclination_vertical_center_and_zero_scatter_angle"
                ),
                "nonadopted_fit_reference_key": "Villaraos2022",
                "nonadopted_fit_reference_doi": "10.1093/mnras/stac2973",
            }
            for field, expected in reviewed_j0437_text.items():
                if row[field] != expected:
                    raise ValueError(
                        f"J0437: {field} does not match the reviewed value"
                    )
            for field in (
                "size_systematic_unc_minus_uas",
                "size_systematic_unc_plus_uas",
                "size_total_unc_minus_uas",
                "size_total_unc_plus_uas",
                "mass_systematic_unc_minus_msun",
                "mass_systematic_unc_plus_msun",
                "size_model_check_excluded_configuration_ids",
            ):
                if row[field] is not None:
                    raise ValueError(f"J0437: {field} must remain blank")
        if row["source_key"] == "Mrk1029":
            reviewed_mrk1029 = {
                "distance_mpc": 119.32882273342354,
                "distance_stat_unc_minus_mpc": 3.4105504657071966,
                "distance_stat_unc_plus_mpc": 3.4105504657071966,
                "distance_systematic_unc_minus_mpc": 4.844201193508398,
                "distance_systematic_unc_plus_mpc": 4.844201193508398,
                "distance_unc_minus_mpc": 5.924368293946939,
                "distance_unc_plus_mpc": 5.924368293946939,
                "size_model_check_accepted_configuration_count": 1.0,
                "size_model_check_lower_uas": 0.0003317155455200028,
                "size_model_check_upper_uas": 0.0005906154834868343,
                "nonadopted_fit_ratio_msun_per_mpc": 14_400.0,
                "nonadopted_fit_ratio_unc_minus_msun_per_mpc": 1_100.0,
                "nonadopted_fit_ratio_unc_plus_msun_per_mpc": 1_100.0,
                "nonadopted_fit_size_uas": 0.0014771315895750299,
                "nonadopted_fit_size_unc_minus_uas": 0.00011283644087031476,
                "nonadopted_fit_size_unc_plus_uas": 0.00011283644087031476,
            }
            for field, expected in reviewed_mrk1029.items():
                actual = row[field]
                if actual is None or not math.isclose(
                    float(actual), expected, rel_tol=0.0, abs_tol=1e-15
                ):
                    raise ValueError(
                        f"Mrk1029: {field} does not match the reviewed value"
                    )
            reviewed_mrk1029_text = {
                "treatment": "no_robust_mass",
                "adopted_uncertainty_model": "not_available",
                "mass_systematic_status": "not_available",
                "mass_display_interval_status": "not_available",
                "distance_stat_component_status": "statistical",
                "size_stat_component_status": "not_available",
                "size_systematic_status": "not_applicable",
                "size_total_status": "not_available",
                "legacy_uncertainty_basis": "not_available",
                "ratio_provenance": "no_robust_primary_mass",
                "uncertainty_status": (
                    "no robust primary black-hole mass; separate non-probabilistic "
                    "standalone disk-model diagnostic not combined; conditional "
                    "model-dependent literature ratio reported separately and not "
                    "adopted"
                ),
                "size_model_check_key": (
                    "deepshikha_2024_mrk1029_mestel_disk_mass_model"
                ),
                "size_model_check_type": "disk_mass_model_choice",
                "size_model_check_data_relation": (
                    "reanalysis_of_published_gao2017_vlbi_data_data_selection_"
                    "unreported"
                ),
                "size_model_check_range_basis": (
                    "reported_symmetric_mass_error_span"
                ),
                "size_model_check_probability_status": "non_probabilistic",
                "size_model_check_adoption_status": (
                    "reported_separately_not_adopted"
                ),
                "size_model_check_reference_key": "Deepshikha2024",
                "size_model_check_reference_doi": (
                    "10.1017/S1743921323002491"
                ),
                "nonadopted_fit_key": (
                    "gonzalez_juarez_2024_mrk1029_gr_ratio_fit"
                ),
                "nonadopted_fit_interval_status": (
                    "bayesian_formal_conditional_on_model"
                ),
                "nonadopted_fit_adoption_status": (
                    "reported_separately_not_adopted"
                ),
                "nonadopted_fit_basis": (
                    "restricted_same_data_fit_conditional_on_fixed_gao2017_"
                    "inclination_vertical_center_and_selected_scatter_angle_"
                    "does_not_resolve_primary_complex_geometry_and_dynamics"
                ),
                "nonadopted_fit_reference_key": "GonzalezJuarez2024",
                "nonadopted_fit_reference_doi": (
                    "10.1051/0004-6361/202450098"
                ),
            }
            for field, expected in reviewed_mrk1029_text.items():
                if row[field] != expected:
                    raise ValueError(
                        f"Mrk1029: {field} does not match the reviewed value"
                    )
            for field in (
                "size_uas",
                "size_unc_adopted_uas",
                "log_mass",
                "mass_msun",
                "mass_fit_unc_minus_msun",
                "mass_fit_unc_plus_msun",
                "mass_shared_distance_unc_minus_msun",
                "mass_shared_distance_unc_plus_msun",
                "mass_total_unc_minus_msun",
                "mass_total_unc_plus_msun",
                "mass_display_unc_minus_msun",
                "mass_display_unc_plus_msun",
                "size_stat_unc_minus_uas",
                "size_stat_unc_plus_uas",
                "size_systematic_unc_minus_uas",
                "size_systematic_unc_plus_uas",
                "size_total_unc_minus_uas",
                "size_total_unc_plus_uas",
                "size_model_check_excluded_configuration_ids",
            ):
                if row[field] is not None:
                    raise ValueError(f"Mrk1029: {field} must remain blank")
        if row["source_key"] == "NGC1320":
            reviewed_ngc1320 = {
                "distance_mpc": 33.783491204330176,
                "distance_stat_unc_minus_mpc": 3.389871133734789,
                "distance_stat_unc_plus_mpc": 3.389871133734789,
                "distance_systematic_unc_minus_mpc": 1.371454311407179,
                "distance_systematic_unc_plus_mpc": 1.371454311407179,
                "distance_unc_minus_mpc": 3.656790017434105,
                "distance_unc_plus_mpc": 3.656790017434105,
                "size_model_check_accepted_configuration_count": 1.0,
                "size_model_check_lower_uas": 0.01103254026501894,
                "size_model_check_upper_uas": 0.014131268281343874,
                "nonadopted_fit_ratio_msun_per_mpc": 149_700.0,
                "nonadopted_fit_ratio_unc_minus_msun_per_mpc": 6_900.0,
                "nonadopted_fit_ratio_unc_plus_msun_per_mpc": 7_100.0,
                "nonadopted_fit_size_uas": 0.015356013816623746,
                "nonadopted_fit_size_unc_minus_uas": 0.0007077922200047018,
                "nonadopted_fit_size_unc_plus_uas": 0.0007283079365265772,
            }
            for field, expected in reviewed_ngc1320.items():
                actual = row[field]
                if actual is None or not math.isclose(
                    float(actual), expected, rel_tol=0.0, abs_tol=1e-15
                ):
                    raise ValueError(
                        f"NGC1320: {field} does not match the reviewed value"
                    )
            reviewed_ngc1320_text = {
                "treatment": "no_robust_mass",
                "adopted_uncertainty_model": "not_available",
                "mass_systematic_status": "not_available",
                "mass_display_interval_status": "not_available",
                "distance_stat_component_status": "statistical",
                "size_stat_component_status": "not_available",
                "size_systematic_status": "not_applicable",
                "size_total_status": "not_available",
                "legacy_uncertainty_basis": "not_available",
                "ratio_provenance": "no_robust_primary_mass",
                "uncertainty_status": (
                    "no robust primary black-hole mass; separate non-probabilistic "
                    "standalone disk-model diagnostic not combined; conditional "
                    "model-dependent literature ratio reported separately and not "
                    "adopted"
                ),
                "size_model_check_key": (
                    "deepshikha_2024_ngc1320_mestel_disk_mass_model"
                ),
                "size_model_check_type": "disk_mass_model_choice",
                "size_model_check_data_relation": (
                    "reanalysis_of_published_gao2017_vlbi_data_data_selection_"
                    "unreported"
                ),
                "size_model_check_range_basis": (
                    "reported_symmetric_mass_error_span"
                ),
                "size_model_check_probability_status": "non_probabilistic",
                "size_model_check_adoption_status": (
                    "reported_separately_not_adopted"
                ),
                "size_model_check_reference_key": "Deepshikha2024",
                "size_model_check_reference_doi": (
                    "10.1017/S1743921323002491"
                ),
                "nonadopted_fit_key": (
                    "gonzalez_juarez_2024_ngc1320_gr_ratio_fit"
                ),
                "nonadopted_fit_interval_status": (
                    "bayesian_formal_conditional_on_model"
                ),
                "nonadopted_fit_adoption_status": (
                    "reported_separately_not_adopted"
                ),
                "nonadopted_fit_basis": (
                    "restricted_same_data_fit_conditional_on_fixed_edge_on_"
                    "inclination_and_selected_scatter_angle_does_not_resolve_"
                    "primary_complex_geometry_and_dynamics"
                ),
                "nonadopted_fit_reference_key": "GonzalezJuarez2024",
                "nonadopted_fit_reference_doi": (
                    "10.1051/0004-6361/202450098"
                ),
            }
            for field, expected in reviewed_ngc1320_text.items():
                if row[field] != expected:
                    raise ValueError(
                        f"NGC1320: {field} does not match the reviewed value"
                    )
            for field in (
                "size_uas",
                "size_unc_adopted_uas",
                "log_mass",
                "mass_msun",
                "mass_fit_unc_minus_msun",
                "mass_fit_unc_plus_msun",
                "mass_shared_distance_unc_minus_msun",
                "mass_shared_distance_unc_plus_msun",
                "mass_total_unc_minus_msun",
                "mass_total_unc_plus_msun",
                "mass_display_unc_minus_msun",
                "mass_display_unc_plus_msun",
                "size_stat_unc_minus_uas",
                "size_stat_unc_plus_uas",
                "size_systematic_unc_minus_uas",
                "size_systematic_unc_plus_uas",
                "size_total_unc_minus_uas",
                "size_total_unc_plus_uas",
                "size_model_check_excluded_configuration_ids",
            ):
                if row[field] is not None:
                    raise ValueError(f"NGC1320: {field} must remain blank")
        if row["source_key"] == "NGC3079":
            reviewed_ngc3079 = {
                "size_uas": 0.011858795677384631,
                "mass_msun": 2_000_000.0,
                "distance_mpc": 17.3,
                "distance_stat_unc_minus_mpc": 1.5,
                "distance_stat_unc_plus_mpc": 1.5,
                "distance_unc_minus_mpc": 1.5,
                "distance_unc_plus_mpc": 1.5,
                "mass_shared_distance_unc_minus_msun": 173_410.40462427746,
                "mass_shared_distance_unc_plus_msun": 173_410.40462427746,
                "size_model_check_lower_uas": 0.011858795677384631,
                "size_model_check_upper_uas": 0.015890786207695407,
            }
            for field, expected in reviewed_ngc3079.items():
                actual = row[field]
                if actual is None or not math.isclose(
                    float(actual), expected, rel_tol=1e-15, abs_tol=1e-15
                ):
                    raise ValueError(
                        f"NGC3079: {field} does not match the reviewed value"
                    )
            reviewed_ngc3079_text = {
                "treatment": "shared_distance",
                "mass_display_interval_status": "not_available",
                "mass_systematic_status": "not_available",
                "distance_stat_component_status": "reported_legacy_tully_fisher",
                "size_stat_component_status": "no_formal_uncertainty",
                "size_systematic_status": "unquantified_model",
                "size_total_status": "incomplete_unquantified_systematic",
                "legacy_uncertainty_basis": "not_available",
                "size_model_check_probability_status": "non_probabilistic",
                "size_model_check_adoption_status": (
                    "reported_separately_not_adopted"
                ),
            }
            for field, expected in reviewed_ngc3079_text.items():
                if row[field] != expected:
                    raise ValueError(
                        f"NGC3079: {field} does not match the reviewed value"
                    )
            for field in (
                "size_unc_adopted_uas",
                "size_stat_unc_minus_uas",
                "size_stat_unc_plus_uas",
                "size_systematic_unc_minus_uas",
                "size_systematic_unc_plus_uas",
                "size_total_unc_minus_uas",
                "size_total_unc_plus_uas",
                "mass_fit_unc_minus_msun",
                "mass_fit_unc_plus_msun",
                "mass_systematic_unc_minus_msun",
                "mass_systematic_unc_plus_msun",
                "mass_total_unc_minus_msun",
                "mass_total_unc_plus_msun",
                "mass_display_unc_minus_msun",
                "mass_display_unc_plus_msun",
                "mass_uncertainty_reference_doi",
                "nonadopted_fit_key",
            ):
                if row[field] is not None:
                    raise ValueError(f"NGC3079: {field} must remain blank")
        if row["source_key"] == "NGC6323":
            reviewed_ngc6323 = {
                "size_uas": 0.00954948892245025,
                "size_posterior_q16_uas": 0.009493370476285018,
                "size_posterior_q50_uas": 0.00954948892245025,
                "size_posterior_q84_uas": 0.009610803521038186,
                "size_stat_unc_minus_uas": 5.611844616523233e-05,
                "size_stat_unc_plus_uas": 6.131459858793598e-05,
                "size_systematic_unc_minus_uas": 9.54948892245025e-05,
                "size_systematic_unc_plus_uas": 9.54948892245025e-05,
                "size_total_unc_minus_uas": 0.00011076350422409036,
                "size_total_unc_plus_uas": 0.00011348459749234575,
                "mass_msun": 10_200_000.0,
                "mass_unc_minus_msun": 2_200_000.0,
                "mass_unc_plus_msun": 3_200_000.0,
                "distance_mpc": 109.4,
                "distance_unc_minus_mpc": 23.4,
                "distance_unc_plus_mpc": 34.2,
                "size_model_check_lower_uas": 0.009459918021333369,
                "size_model_check_upper_uas": 0.009536827923945834,
                "nonadopted_fit_ratio_msun_per_mpc": 91_600.0,
                "nonadopted_fit_ratio_unc_minus_msun_per_mpc": 500.0,
                "nonadopted_fit_ratio_unc_plus_msun_per_mpc": 500.0,
                "nonadopted_fit_size_uas": 0.00939619816701894,
                "nonadopted_fit_size_unc_minus_uas": 5.128929130468853e-05,
                "nonadopted_fit_size_unc_plus_uas": 5.128929130468853e-05,
            }
            for field, expected in reviewed_ngc6323.items():
                actual = row[field]
                if actual is None or not math.isclose(
                    float(actual), expected, rel_tol=0.0, abs_tol=1e-15
                ):
                    raise ValueError(
                        f"NGC6323: {field} does not match the reviewed value"
                    )
            reviewed_ngc6323_text = {
                "treatment": "paired_posterior",
                "mass_display_interval_status": "decomposed_components",
                "distance_stat_component_status": "statistical",
                "size_stat_component_status": "paired_chain_quantile_interval",
                "size_systematic_status": (
                    "numeric_one_percent_non_keplerian_allowance"
                ),
                "size_total_status": "quadrature_statistical_plus_systematic",
                "size_model_check_key": (
                    "kuo_2018_ngc6323_physically_allowed_disk_mass_models"
                ),
                "size_model_check_probability_status": "non_probabilistic",
                "size_model_check_adoption_status": (
                    "reported_separately_not_adopted"
                ),
                "nonadopted_fit_key": "villaraos_2022_ngc6323_gr_ratio_fit",
                "nonadopted_fit_interval_status": (
                    "bayesian_formal_conditional_on_model"
                ),
                "nonadopted_fit_adoption_status": (
                    "reported_separately_not_adopted"
                ),
                "nonadopted_fit_reference_key": "Villaraos2022",
                "nonadopted_fit_reference_doi": "10.1093/mnras/stac2973",
            }
            for field, expected in reviewed_ngc6323_text.items():
                if row[field] != expected:
                    raise ValueError(
                        f"NGC6323: {field} does not match the reviewed value"
                    )
        if row["source_key"] == "CGCG074-064":
            reviewed_cgcg = {
                "size_uas": 0.02834708992683376,
                "size_posterior_q16_uas": 0.02816418536155449,
                "size_posterior_q50_uas": 0.02834708992683376,
                "size_posterior_q84_uas": 0.0285518183322884,
                "mass_msun": 24_200_000.0,
                "mass_unc_minus_msun": 2_000_000.0,
                "mass_unc_plus_msun": 2_200_000.0,
                "log_mass_unc": 0.03762547620959619,
                "mass_systematic_unc_minus_msun": 500_000.0,
                "mass_systematic_unc_plus_msun": 500_000.0,
                "log_mass_systematic_unc_minus": 0.009067019970327397,
                "log_mass_systematic_unc_plus": 0.008881587279234472,
                "mass_total_unc_minus_msun": 2_061_552.8128088303,
                "mass_total_unc_plus_msun": 2_256_102.8345356956,
                "mass_display_unc_minus_msun": 2_000_000.0,
                "mass_display_unc_plus_msun": 2_200_000.0,
                "nonadopted_fit_ratio_msun_per_mpc": 274_900.0,
                "nonadopted_fit_ratio_unc_minus_msun_per_mpc": 1_500.0,
                "nonadopted_fit_ratio_unc_plus_msun_per_mpc": 1_700.0,
                "nonadopted_fit_size_uas": 0.02819885235931776,
                "nonadopted_fit_size_unc_minus_uas": 0.00015386787391406562,
                "nonadopted_fit_size_unc_plus_uas": 0.00017438359043594104,
            }
            for field, expected in reviewed_cgcg.items():
                actual = row[field]
                if actual is None or not math.isclose(
                    float(actual), expected, rel_tol=0.0, abs_tol=1e-15
                ):
                    raise ValueError(
                        f"CGCG074-064: {field} does not match the reviewed value"
                    )
            reviewed_cgcg_text = {
                "mass_systematic_status": (
                    "source_reported_model_selection_systematic"
                ),
                "mass_display_interval_status": (
                    "separate_statistical_and_systematic_components"
                ),
                "nonadopted_fit_key": "villaraos_2022_cgcg074064_gr_ratio_fit",
                "nonadopted_fit_interval_status": (
                    "bayesian_formal_conditional_on_model"
                ),
                "nonadopted_fit_adoption_status": (
                    "reported_separately_not_adopted"
                ),
                "nonadopted_fit_reference_key": "Villaraos2022",
                "nonadopted_fit_reference_doi": "10.1093/mnras/stac2973",
            }
            for field, expected in reviewed_cgcg_text.items():
                if row[field] != expected:
                    raise ValueError(
                        f"CGCG074-064: {field} does not match the reviewed value"
                    )
        elif row["source_key"] == "UGC6093":
            reviewed_ugc6093 = {
                "size_uas": 0.017445797174172238,
                "size_stat_unc_minus_uas": 0.0006085743200292641,
                "size_stat_unc_plus_uas": 0.000743813057813545,
                "size_systematic_unc_minus_uas": 0.000872289858708612,
                "size_systematic_unc_plus_uas": 0.000872289858708612,
                "size_total_unc_minus_uas": 0.0010636034508241178,
                "size_total_unc_plus_uas": 0.0011463627098697107,
                "mass_msun": 24_443_424.32138069,
                "mass_fit_unc_minus_msun": 852_677.5926063032,
                "mass_fit_unc_plus_msun": 1_042_161.5020743705,
                "distance_stat_unc_minus_mpc": 3.386942772210658,
                "distance_stat_unc_plus_mpc": 3.386105159165424,
                "distance_unc_minus_mpc": 6.746328117856665,
                "distance_unc_plus_mpc": 6.7459076394863375,
                "mass_shared_distance_unc_minus_msun": 1_147_364.9666493207,
                "mass_shared_distance_unc_plus_msun": 1_147_293.4548368326,
                "mass_systematic_unc_minus_msun": 1_222_171.2160690345,
                "mass_systematic_unc_plus_msun": 1_222_171.2160690345,
                "log_mass_systematic_unc_minus": 0.02227639471115221,
                "log_mass_systematic_unc_plus": 0.021189299069938092,
                "mass_total_unc_minus_msun": 1_880_746.6403039885,
                "mass_total_unc_plus_msun": 1_973_850.386758076,
                "mass_display_unc_minus_msun": 1_429_512.309715127,
                "mass_display_unc_plus_msun": 1_549_962.214996657,
            }
            for field, expected in reviewed_ugc6093.items():
                actual = row[field]
                if actual is None or not math.isclose(
                    float(actual), expected, rel_tol=1e-15, abs_tol=1e-15
                ):
                    raise ValueError(
                        f"UGC6093: {field} does not match the reviewed value"
                    )
            reviewed_ugc6093_text = {
                "mass_systematic_status": (
                    "source_reported_absolute_position_systematic"
                ),
                "mass_display_interval_status": (
                    "separate_decomposed_and_position_systematic_components"
                ),
                "size_stat_component_status": "fixed_distance_formal_fit_interval",
                "size_systematic_status": "numeric_systematic",
                "size_total_status": "quadrature_statistical_plus_systematic",
            }
            for field, expected in reviewed_ugc6093_text.items():
                if row[field] != expected:
                    raise ValueError(
                        f"UGC6093: {field} does not match the reviewed value"
                    )
        elif row["source_key"] == "ESO558-G009":
            reviewed_eso558 = {
                "size_uas": 0.01620665338623987,
                "size_stat_unc_minus_uas": 0.00037275302788351703,
                "size_stat_unc_plus_uas": 0.00037275302788351703,
                "size_systematic_unc_minus_uas": 0.0009729392749756379,
                "size_systematic_unc_plus_uas": 0.0009729392749756379,
                "size_total_unc_minus_uas": 0.0010419000204369179,
                "size_total_unc_plus_uas": 0.0010419000204369179,
                "distance_stat_unc_minus_mpc": 3.4026220986146067,
                "distance_stat_unc_plus_mpc": 3.4026220986146067,
                "distance_unc_minus_mpc": 5.49508167148691,
                "distance_unc_plus_mpc": 5.49508167148691,
                "mass_msun": 16_792_963.464140728,
                "mass_fit_unc_minus_msun": 386_238.15967523673,
                "mass_fit_unc_plus_msun": 386_238.15967523673,
                "mass_shared_distance_unc_minus_msun": 868_182.0484691216,
                "mass_shared_distance_unc_plus_msun": 868_182.0484691216,
                "mass_systematic_unc_minus_msun": 1_008_137.4178931694,
                "mass_systematic_unc_plus_msun": 1_008_137.4178931694,
                "log_mass_systematic_unc_minus": 0.02688754291121219,
                "log_mass_systematic_unc_plus": 0.02531951831781824,
                "mass_total_unc_minus_msun": 1_385_373.9706771097,
                "mass_total_unc_plus_msun": 1_385_373.9706771097,
                "mass_display_unc_minus_msun": 950_221.0191704632,
                "mass_display_unc_plus_msun": 950_221.0191704632,
                "nonadopted_fit_ratio_msun_per_mpc": 159_400.0,
                "nonadopted_fit_ratio_unc_minus_msun_per_mpc": 2_600.0,
                "nonadopted_fit_ratio_unc_plus_msun_per_mpc": 2_500.0,
                "nonadopted_fit_size_uas": 0.016351026067934704,
                "nonadopted_fit_size_unc_minus_uas": 0.00026670431478438035,
                "nonadopted_fit_size_unc_plus_uas": 0.00025644645652344266,
            }
            for field, expected in reviewed_eso558.items():
                actual = row[field]
                if actual is None or not math.isclose(
                    float(actual), expected, rel_tol=1e-15, abs_tol=1e-15
                ):
                    raise ValueError(
                        f"ESO558-G009: {field} does not match the reviewed value"
                    )
            reviewed_eso558_text = {
                "mass_systematic_status": (
                    "source_reported_midline_inclination_systematic"
                ),
                "mass_display_interval_status": (
                    "separate_decomposed_and_geometry_systematic_components"
                ),
                "size_stat_component_status": (
                    "fit_plus_absolute_position_mixed_component"
                ),
                "size_systematic_status": "numeric_systematic",
                "size_total_status": (
                    "quadrature_mixed_fit_position_plus_systematic"
                ),
                "nonadopted_fit_key": (
                    "villaraos_2022_eso558g009_gr_ratio_fit"
                ),
                "nonadopted_fit_interval_status": (
                    "bayesian_formal_conditional_on_model"
                ),
                "nonadopted_fit_adoption_status": (
                    "reported_separately_not_adopted"
                ),
                "nonadopted_fit_reference_key": "Villaraos2022",
                "nonadopted_fit_reference_doi": "10.1093/mnras/stac2973",
            }
            for field, expected in reviewed_eso558_text.items():
                if row[field] != expected:
                    raise ValueError(
                        f"ESO558-G009: {field} does not match the reviewed value"
                    )
        elif (
            row["mass_systematic_status"] != "not_available"
            or row["mass_systematic_unc_minus_msun"] is not None
            or row["mass_systematic_unc_plus_msun"] is not None
            or row["log_mass_systematic_unc_minus"] is not None
            or row["log_mass_systematic_unc_plus"] is not None
        ):
            raise ValueError(
                f"{row['source_key']}: unexpected source-reported mass systematic"
            )
        writer.writerow({field: format_value(row[field]) for field in OUTPUT_FIELDS})
    return buffer.getvalue()


def require_final_release(rows: list[dict[str, object]]) -> None:
    """Reject a final release while an explicitly requested term is pending."""

    pending = sorted(
        str(row["source_label"])
        for row in rows
        if row["size_systematic_status"] == "pending_nakai"
    )
    if pending:
        raise RuntimeError(
            "final black-hole parameter release awaits an explicit formal-only "
            f"or systematic decision for: {', '.join(pending)}"
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--quantiles", type=Path, default=DEFAULT_QUANTILES)
    parser.add_argument("--model-checks", type=Path, default=DEFAULT_MODEL_CHECKS)
    parser.add_argument(
        "--audit-only-model-checks",
        type=Path,
        default=DEFAULT_AUDIT_ONLY_MODEL_CHECKS,
    )
    parser.add_argument(
        "--nonadopted-fits", type=Path, default=DEFAULT_NONADOPTED_FITS
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if the existing output is not the exact deterministic build",
    )
    parser.add_argument(
        "--stdout",
        action="store_true",
        help="write the generated CSV to stdout instead of a file",
    )
    parser.add_argument(
        "--require-final",
        action="store_true",
        help="fail before output/check while an explicitly requested term is pending",
    )
    args = parser.parse_args(argv)

    rows = build_rows(
        args.input,
        args.quantiles,
        args.model_checks,
        args.nonadopted_fits,
        args.audit_only_model_checks,
    )
    if args.require_final:
        try:
            require_final_release(rows)
        except RuntimeError as error:
            print(f"final-release gate failed: {error}", file=sys.stderr)
            return 2

    rendered = render_csv(rows)
    if args.stdout:
        sys.stdout.write(rendered)
        return 0
    if args.check:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != rendered:
            print(f"stale generated table: {args.output}", file=sys.stderr)
            return 1
        print(f"verified deterministic table: {args.output}")
        return 0
    args.output.write_text(rendered, encoding="utf-8", newline="")
    print(f"wrote {len(rendered.splitlines()) - 1} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
