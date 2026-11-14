#!/usr/bin/env python3
"""Render manuscript Table 1 deterministically from ``BHS_params.csv``.

The calculation pipeline owns every mass, distance, and shadow-size number.
This module owns only presentation metadata (source typography, inclination
cells, and BibTeX keys).  The comparison rows are fixed literature values and
are deliberately kept separate from the 21-source machine table.
"""

from __future__ import annotations

import argparse
import csv
import math
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "analysis" / "CountingBHS_and_Spectra"
DEFAULT_PARAMS = DATA_DIR / "BHS_params.csv"
DEFAULT_METADATA = DATA_DIR / "bhs_latex_metadata.csv"
DEFAULT_MANUSCRIPT_OUTPUT = (
    REPO_ROOT.parent / "revision_package" / "working_clean" / "BHS_Params.tex"
)

DISPLAY_SOURCE_KEYS = (
    "NGC4258",
    "NGC1194",
    "NGC1068",
    "NGC3393",
    "NGC4388",
    "Circinus",
    "NGC4945",
    "NGC5765b",
    "NGC2273",
    "CGCG074-064",
    "UGC3789",
    "NGC6264",
    "NGC2960",
    "UGC6093",
    "ESO558-G009",
    "NGC3079",
    "NGC6323",
    "IC2560",
    "J0437",
    "Mrk1029",
    "NGC1320",
)

METADATA_FIELDS = [
    "source_order",
    "source_key",
    "source_tex",
    "inclination_tex",
    "mass_reference_key",
    "distance_reference_key",
    "inclination_reference_key",
    "size_uncertainty_marker",
]

EXPECTED_DAGGERED_SAMPLE_KEYS = frozenset(
    {"NGC1194", "NGC2273", "ESO558-G009"}
)

REQUIRED_PARAM_FIELDS = {
    "source_order",
    "source_key",
    "size_uas",
    "mass_msun",
    "mass_total_unc_minus_msun",
    "mass_total_unc_plus_msun",
    "mass_display_unc_minus_msun",
    "mass_display_unc_plus_msun",
    "mass_display_interval_status",
    "mass_systematic_unc_minus_msun",
    "mass_systematic_unc_plus_msun",
    "log_mass_systematic_unc_minus",
    "log_mass_systematic_unc_plus",
    "mass_systematic_status",
    "distance_mpc",
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
    "size_systematic_status",
    "treatment",
    "parameter_reference_doi",
    "mass_uncertainty_reference_doi",
    "distance_reference_doi",
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
}

# This mapping makes the reference-cell metadata auditable against the DOI
# fields produced by build_bhs_params.py.
REFERENCE_DOIS = {
    "Reid2019": "10.3847/2041-8213/ab552d",
    "Kuo2010": "10.1088/0004-637X/727/1/20",
    "Gallimore2024": "10.3847/2041-8213/ad864f",
    "Markham2026": "10.3847/1538-4357/ae47d8",
    "Kondratko2008": "10.1086/586879",
    "Greenhill2003": "10.1086/374862",
    "Pesce2020b": "10.3847/2041-8213/ab75f0",
    "Greenhill1997": "10.1086/310643",
    "Pesce2020a": "10.3847/1538-4357/ab6bcd",
    "Zhao2018": "10.3847/1538-4357/aaa95c",
    "Gao2016b": "10.3847/1538-4357/834/1/52",
    "Kuo2018Disk": "10.3847/1538-4357/aabff1",
    "Nakai2024": "10.1017/S1743921323002971",
    "Kondratko2005": "10.1086/426101",
    "Pesce2021J0437": "10.3847/1538-4357/abde3d",
    "vanDerBoom2026": "10.1093/mnras/stag1327",
    "GonzalezJuarez2024": "10.1051/0004-6361/202450098",
    "Villaraos2022": "10.1093/mnras/stac2973",
    "Kuo2020": "10.1093/mnras/staa2260",
    "Deepshikha2024": "10.1017/S1743921323002491",
}

NGC4258_MODEL_CHECK_LOWER_UAS = 0.53679114913867576
NGC4258_MODEL_CHECK_UPPER_UAS = 0.54888539817298243
NGC3079_MODEL_CHECK_LOWER_UAS = 0.011858795677384631
NGC3079_MODEL_CHECK_UPPER_UAS = 0.015890786207695407
NGC4945_MODEL_CHECK_LOWER_UAS = 0.03881351774408862
NGC4945_MODEL_CHECK_UPPER_UAS = 0.04812876200266989
NGC4388_MODEL_CHECK_LOWER_UAS = 0.04535053125888249
NGC4388_MODEL_CHECK_UPPER_UAS = 0.04539588179014137
CIRCINUS_MODEL_CHECK_LOWER_UAS = 0.04151990248474786
CIRCINUS_MODEL_CHECK_UPPER_UAS = 0.0427654995592903
NGC2960_MODEL_CHECK_LOWER_UAS = 0.017503713819217805
NGC2960_MODEL_CHECK_UPPER_UAS = 0.017574751618808785
NGC6323_MODEL_CHECK_LOWER_UAS = 0.009459918021333369
NGC6323_MODEL_CHECK_UPPER_UAS = 0.009536827923945834
J0437_MODEL_CHECK_LOWER_UAS = 0.002434118995564536
J0437_MODEL_CHECK_UPPER_UAS = 0.0031252885868976762
MRK1029_MODEL_CHECK_LOWER_UAS = 0.0003317155455200028
MRK1029_MODEL_CHECK_UPPER_UAS = 0.0005906154834868343
NGC1320_MODEL_CHECK_LOWER_UAS = 0.01103254026501894
NGC1320_MODEL_CHECK_UPPER_UAS = 0.014131268281343874

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
DISTANCE_STAT_LABELS = {
    "statistical": "stat.",
    "reported_legacy_multimethod": "reported",
    "reported_legacy_tully_fisher": "reported",
    "not_available": None,
}

REVIEWED_NONADOPTED_FITS = {
    "NGC1194": {
        "fit_key": "gonzalez_juarez_2024_ngc1194_gr_ratio_fit",
        "ratio": 1_310_000.0,
        "ratio_minus": 20_900.0,
        "ratio_plus": 21_100.0,
        "size": 0.13437794321828395,
        "size_minus": 0.002143892376535981,
        "size_plus": 0.002164408093057856,
        "primary_treatment": "shared_distance",
        "primary_result_status": "retained",
        "basis": (
            "same_data_reanalysis_conditional_on_restricted_disk_model_and_"
            "discrepant_recession_velocity"
        ),
        "reference_key": "GonzalezJuarez2024",
    },
    "NGC2273": {
        "fit_key": "villaraos_2022_ngc2273_gr_ratio_fit",
        "ratio": 328_300.0,
        "ratio_minus": 5_100.0,
        "ratio_plus": 4_500.0,
        "size": 0.033676548670658495,
        "size_minus": 0.0005231507713078231,
        "size_plus": 0.0004616036217421969,
        "primary_treatment": "shared_distance",
        "primary_result_status": "retained",
        "basis": (
            "same_data_reanalysis_conditional_on_fixed_vertical_center_and_"
            "selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
    },
    "NGC4388": {
        "fit_key": "villaraos_2022_ngc4388_gr_ratio_fit",
        "ratio": 423_300.0,
        "ratio_minus": 26_800.0,
        "ratio_plus": 18_500.0,
        "size": 0.043421514018549316,
        "size_minus": 0.0027491060139313056,
        "size_plus": 0.0018977037782734758,
        "primary_treatment": "shared_distance",
        "primary_result_status": "retained",
        "basis": (
            "same_data_reanalysis_conditional_on_fixed_vertical_center_and_"
            "selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
    },
    "CGCG074-064": {
        "fit_key": "villaraos_2022_cgcg074064_gr_ratio_fit",
        "ratio": 274_900.0,
        "ratio_minus": 1_500.0,
        "ratio_plus": 1_700.0,
        "size": 0.02819885235931776,
        "size_minus": 0.00015386787391406562,
        "size_plus": 0.00017438359043594104,
        "primary_treatment": "paired_posterior",
        "primary_result_status": "retained",
        "basis": (
            "same_data_reanalysis_conditional_on_fixed_pesce2020a_inclination_and_"
            "selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
    },
    "NGC6264": {
        "fit_key": "villaraos_2022_ngc6264_gr_ratio_fit",
        "ratio": 212_600.0,
        "ratio_minus": 700.0,
        "ratio_plus": 700.0,
        "size": 0.021808206662753567,
        "size_minus": 0.00007180500782656395,
        "size_plus": 0.00007180500782656395,
        "primary_treatment": "paired_posterior",
        "primary_result_status": "retained",
        "basis": (
            "older_partial_data_reanalysis_conditional_on_fixed_vertical_center_and_"
            "selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
    },
    "UGC3789": {
        "fit_key": "villaraos_2022_ugc3789_gr_ratio_fit",
        "ratio": 227_700.0,
        "ratio_minus": 900.0,
        "ratio_plus": 900.0,
        "size": 0.02335714326015516,
        "size_minus": 0.00009232072434843937,
        "size_plus": 0.00009232072434843937,
        "primary_treatment": "paired_posterior",
        "primary_result_status": "retained",
        "basis": (
            "same_observation_vlbi_spot_subset_reanalysis_conditional_on_fixed_"
            "reid2013_inclination_and_selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
    },
    "ESO558-G009": {
        "fit_key": "villaraos_2022_eso558g009_gr_ratio_fit",
        "ratio": 159_400.0,
        "ratio_minus": 2_600.0,
        "ratio_plus": 2_500.0,
        "size": 0.016351026067934704,
        "size_minus": 0.00026670431478438035,
        "size_plus": 0.00025644645652344266,
        "primary_treatment": "shared_distance",
        "primary_result_status": "retained",
        "basis": (
            "same_data_reanalysis_conditional_on_fixed_gao2017_inclination_and_"
            "selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
    },
    "NGC2960": {
        "fit_key": "villaraos_2022_ngc2960_gr_ratio_fit",
        "ratio": 173_800.0,
        "ratio_minus": 1_600.0,
        "ratio_plus": 1_500.0,
        "size": 0.017828157657509735,
        "size_minus": 0.0001641257321750033,
        "size_plus": 0.0001538678739140656,
        "primary_treatment": "shared_distance",
        "primary_result_status": "retained",
        "basis": (
            "older_published_data_reanalysis_conditional_on_fixed_kuo2011_"
            "inclination_vertical_center_and_selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
    },
    "NGC5765b": {
        "fit_key": "villaraos_2022_ngc5765b_gr_ratio_fit",
        "ratio": 372_700.0,
        "ratio_minus": 1_300.0,
        "ratio_plus": 1_300.0,
        "size": 0.038231037738514835,
        "size_minus": 0.00013335215739219018,
        "size_plus": 0.00013335215739219018,
        "primary_treatment": "paired_posterior",
        "primary_result_status": "retained",
        "basis": (
            "same_data_reanalysis_conditional_on_fixed_gao2016_inclination_and_"
            "selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
    },
    "NGC6323": {
        "fit_key": "villaraos_2022_ngc6323_gr_ratio_fit",
        "ratio": 91_600.0,
        "ratio_minus": 500.0,
        "ratio_plus": 500.0,
        "size": 0.00939619816701894,
        "size_minus": 5.128929130468853e-05,
        "size_plus": 5.128929130468853e-05,
        "primary_treatment": "paired_posterior",
        "primary_result_status": "retained",
        "basis": (
            "older_partial_data_reanalysis_conditional_on_fixed_kuo2015_"
            "inclination_and_selected_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
    },
    "J0437": {
        "fit_key": "villaraos_2022_j0437_gr_ratio_fit",
        "ratio": 42_500.0,
        "ratio_minus": 1_200.0,
        "ratio_plus": 1_000.0,
        "size": 0.004359589760898526,
        "size_minus": 0.0001230942991312525,
        "size_plus": 0.0001025785826093771,
        "primary_treatment": "shared_distance",
        "primary_result_status": "retained",
        "basis": (
            "same_data_reanalysis_conditional_on_fixed_gao2017_inclination_"
            "vertical_center_and_zero_scatter_angle"
        ),
        "reference_key": "Villaraos2022",
    },
    "Mrk1029": {
        "fit_key": "gonzalez_juarez_2024_mrk1029_gr_ratio_fit",
        "ratio": 14_400.0,
        "ratio_minus": 1_100.0,
        "ratio_plus": 1_100.0,
        "size": 0.0014771315895750299,
        "size_minus": 0.00011283644087031476,
        "size_plus": 0.00011283644087031476,
        "primary_treatment": "no_robust_mass",
        "primary_result_status": "absent",
        "basis": (
            "restricted_same_data_fit_conditional_on_fixed_gao2017_inclination_"
            "vertical_center_and_selected_scatter_angle_does_not_resolve_primary_"
            "complex_geometry_and_dynamics"
        ),
        "reference_key": "GonzalezJuarez2024",
    },
    "NGC1320": {
        "fit_key": "gonzalez_juarez_2024_ngc1320_gr_ratio_fit",
        "ratio": 149_700.0,
        "ratio_minus": 6_900.0,
        "ratio_plus": 7_100.0,
        "size": 0.015356013816623746,
        "size_minus": 0.0007077922200047018,
        "size_plus": 0.0007283079365265772,
        "primary_treatment": "no_robust_mass",
        "primary_result_status": "absent",
        "basis": (
            "restricted_same_data_fit_conditional_on_fixed_edge_on_inclination_"
            "and_selected_scatter_angle_does_not_resolve_primary_complex_"
            "geometry_and_dynamics"
        ),
        "reference_key": "GonzalezJuarez2024",
    },
}

# Reference order is generated from first appearance in the table.  These two
# rows are literature comparisons rather than members of BHS_params.csv.
COMPARISON_ROWS = [
    {
        "source_tex": "Sgr~A*",
        "size_tex": r"51.8 \pm 2.3",
        "size_uncertainty_marker": "",
        "size_reference_key": "EventHorizonTelescopeCollaboration2022d",
        "mass_tex": r"6.6332 \pm 0.0012",
        "mass_reference_key": "GRAVITYCollaboration2022",
        "distance_tex": r"0.00815 \pm 0.00015",
        "distance_reference_key": "EventHorizonTelescopeCollaboration2022d",
        "inclination_tex": r"\leq 30",
        "inclination_reference_key": "EventHorizonTelescopeCollaboration2022f",
    },
    {
        "source_tex": "M87*",
        "size_tex": r"42 \pm 3",
        "size_uncertainty_marker": "dagger",
        "size_reference_key": "EventHorizonTelescopeCollaboration2019f",
        "mass_tex": r"9.74 \pm 0.03",
        "mass_reference_key": "Simon2024",
        "distance_tex": r"16.8 \pm 0.8",
        "distance_reference_key": "EventHorizonTelescopeCollaboration2019f",
        "inclination_tex": r"\sim 17",
        "inclination_reference_key": "EventHorizonTelescopeCollaboration2021b",
    },
]

INDIVIDUALLY_DISCUSSED_SOURCE_KEYS = frozenset(
    {
        "NGC4258",
        "NGC1194",
        "NGC1068",
        "NGC3079",
        "NGC4945",
        "Circinus",
    }
)


def optional_float(value: str | None) -> float | None:
    if value is None or value.strip() in {"", "NA"}:
        return None
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"non-finite value: {value!r}")
    return parsed


def _load_dict_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        rows = list(reader)
    if any(None in row for row in rows):
        raise ValueError(f"one or more rows in {path} contain extra columns")
    return fields, rows


def load_params(path: Path = DEFAULT_PARAMS) -> list[dict[str, str]]:
    fields, rows = _load_dict_rows(path)
    missing = REQUIRED_PARAM_FIELDS.difference(fields)
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)!r}")
    if len(rows) != 21:
        raise ValueError(f"expected 21 parameter rows, found {len(rows)}")
    rows.sort(key=lambda row: int(row["source_order"]))
    if [int(row["source_order"]) for row in rows] != list(range(1, 22)):
        raise ValueError("parameter source_order must be the sequence 1..21")
    return rows


def load_metadata(path: Path = DEFAULT_METADATA) -> list[dict[str, str]]:
    fields, rows = _load_dict_rows(path)
    if fields != METADATA_FIELDS:
        raise ValueError(
            f"unexpected metadata columns: expected {METADATA_FIELDS!r}, found {fields!r}"
        )
    if len(rows) != 21:
        raise ValueError(f"expected 21 metadata rows, found {len(rows)}")
    rows.sort(key=lambda row: int(row["source_order"]))
    if [int(row["source_order"]) for row in rows] != list(range(1, 22)):
        raise ValueError("metadata source_order must be the sequence 1..21")
    return rows


def validate_inputs(
    params: list[dict[str, str]], metadata: list[dict[str, str]]
) -> None:
    if [row["source_key"] for row in params] != [
        row["source_key"] for row in metadata
    ]:
        raise ValueError("parameter and metadata source inventories do not match")

    daggered = {
        row["source_key"]
        for row in metadata
        if row["size_uncertainty_marker"] == "dagger"
    }
    invalid_markers = {
        row["size_uncertainty_marker"]
        for row in metadata
        if row["size_uncertainty_marker"] not in {"", "dagger"}
    }
    if invalid_markers:
        raise ValueError(
            f"unsupported size uncertainty markers: {sorted(invalid_markers)!r}"
        )
    if daggered != EXPECTED_DAGGERED_SAMPLE_KEYS:
        raise ValueError(
            "size uncertainty dagger set must match the accepted presentation: "
            f"expected {sorted(EXPECTED_DAGGERED_SAMPLE_KEYS)!r}, "
            f"found {sorted(daggered)!r}"
        )

    for param, meta in zip(params, metadata, strict=True):
        key = param["source_key"]
        mass_display_status = param["mass_display_interval_status"]
        if mass_display_status not in MASS_DISPLAY_INTERVAL_STATUSES:
            raise ValueError(f"{key}: unsupported mass display interval status")
        display_mass_minus = optional_float(param["mass_display_unc_minus_msun"])
        display_mass_plus = optional_float(param["mass_display_unc_plus_msun"])
        if (display_mass_minus is None) != (display_mass_plus is None):
            raise ValueError(f"{key}: incomplete mass display interval")
        mass_systematic_status = param["mass_systematic_status"]
        if mass_systematic_status not in MASS_SYSTEMATIC_STATUSES:
            raise ValueError(f"{key}: unsupported mass systematic status")
        mass_systematic_minus = optional_float(
            param["mass_systematic_unc_minus_msun"]
        )
        mass_systematic_plus = optional_float(
            param["mass_systematic_unc_plus_msun"]
        )
        log_mass_systematic_minus = optional_float(
            param["log_mass_systematic_unc_minus"]
        )
        log_mass_systematic_plus = optional_float(
            param["log_mass_systematic_unc_plus"]
        )
        if (mass_systematic_minus is None) != (mass_systematic_plus is None) or (
            (log_mass_systematic_minus is None)
            != (log_mass_systematic_plus is None)
        ):
            raise ValueError(f"{key}: incomplete mass systematic interval")
        if mass_systematic_status == "not_available":
            if any(
                value is not None
                for value in (
                    mass_systematic_minus,
                    mass_systematic_plus,
                    log_mass_systematic_minus,
                    log_mass_systematic_plus,
                )
            ):
                raise ValueError(
                    f"{key}: unavailable mass systematic must remain blank"
                )
        elif mass_systematic_status == "source_reported_model_selection_systematic":
            mass = optional_float(param["mass_msun"])
            total_minus = optional_float(param["mass_total_unc_minus_msun"])
            total_plus = optional_float(param["mass_total_unc_plus_msun"])
            if (
                key != "CGCG074-064"
                or mass_display_status
                != "separate_statistical_and_systematic_components"
                or mass != 24_200_000.0
                or mass_systematic_minus != 500_000.0
                or mass_systematic_plus != 500_000.0
                or log_mass_systematic_minus is None
                or log_mass_systematic_plus is None
                or not math.isclose(
                    log_mass_systematic_minus,
                    math.log10(24_200_000.0 / 23_700_000.0),
                    rel_tol=0.0,
                    abs_tol=1e-15,
                )
                or not math.isclose(
                    log_mass_systematic_plus,
                    math.log10(24_700_000.0 / 24_200_000.0),
                    rel_tol=0.0,
                    abs_tol=1e-15,
                )
                or total_minus is None
                or total_plus is None
                or not math.isclose(
                    total_minus,
                    math.hypot(2_000_000.0, 500_000.0),
                    rel_tol=0.0,
                    abs_tol=1e-9,
                )
                or not math.isclose(
                    total_plus,
                    math.hypot(2_200_000.0, 500_000.0),
                    rel_tol=0.0,
                    abs_tol=1e-9,
                )
            ):
                raise ValueError(
                    f"{key}: mass model-selection systematic does not match "
                    "the reviewed source values"
                )
        elif mass_systematic_status == "source_reported_absolute_position_systematic":
            mass = optional_float(param["mass_msun"])
            total_minus = optional_float(param["mass_total_unc_minus_msun"])
            total_plus = optional_float(param["mass_total_unc_plus_msun"])
            if (
                key != "UGC6093"
                or mass_display_status
                != "separate_decomposed_and_position_systematic_components"
                or mass is None
                or not math.isclose(
                    mass, 24_443_424.32138069, rel_tol=1e-15, abs_tol=1e-9
                )
                or mass_systematic_minus is None
                or mass_systematic_plus is None
                or not math.isclose(
                    mass_systematic_minus,
                    1_222_171.2160690345,
                    rel_tol=1e-15,
                    abs_tol=1e-9,
                )
                or not math.isclose(
                    mass_systematic_plus,
                    1_222_171.2160690345,
                    rel_tol=1e-15,
                    abs_tol=1e-9,
                )
                or log_mass_systematic_minus is None
                or log_mass_systematic_plus is None
                or not math.isclose(
                    log_mass_systematic_minus,
                    0.02227639471115221,
                    rel_tol=0.0,
                    abs_tol=1e-15,
                )
                or not math.isclose(
                    log_mass_systematic_plus,
                    0.021189299069938092,
                    rel_tol=0.0,
                    abs_tol=1e-15,
                )
                or total_minus is None
                or total_plus is None
                or not math.isclose(
                    total_minus,
                    1_880_746.6403039885,
                    rel_tol=1e-15,
                    abs_tol=1e-9,
                )
                or not math.isclose(
                    total_plus,
                    1_973_850.386758076,
                    rel_tol=1e-15,
                    abs_tol=1e-9,
                )
            ):
                raise ValueError(
                    f"{key}: absolute-position mass systematic does not match "
                    "the reviewed source values"
                )
        elif (
            mass_systematic_status
            == "source_reported_midline_inclination_systematic"
        ):
            mass = optional_float(param["mass_msun"])
            total_minus = optional_float(param["mass_total_unc_minus_msun"])
            total_plus = optional_float(param["mass_total_unc_plus_msun"])
            if (
                key != "ESO558-G009"
                or mass_display_status
                != "separate_decomposed_and_geometry_systematic_components"
                or mass is None
                or not math.isclose(
                    mass, 16_792_963.464140728, rel_tol=1e-15, abs_tol=1e-9
                )
                or mass_systematic_minus is None
                or mass_systematic_plus is None
                or not math.isclose(
                    mass_systematic_minus,
                    1_008_137.4178931694,
                    rel_tol=1e-15,
                    abs_tol=1e-9,
                )
                or not math.isclose(
                    mass_systematic_plus,
                    1_008_137.4178931694,
                    rel_tol=1e-15,
                    abs_tol=1e-9,
                )
                or log_mass_systematic_minus is None
                or log_mass_systematic_plus is None
                or not math.isclose(
                    log_mass_systematic_minus,
                    0.02688754291121219,
                    rel_tol=0.0,
                    abs_tol=1e-15,
                )
                or not math.isclose(
                    log_mass_systematic_plus,
                    0.02531951831781824,
                    rel_tol=0.0,
                    abs_tol=1e-15,
                )
                or total_minus is None
                or total_plus is None
                or not math.isclose(
                    total_minus,
                    1_385_373.9706771097,
                    rel_tol=1e-15,
                    abs_tol=1e-9,
                )
                or not math.isclose(
                    total_plus,
                    1_385_373.9706771097,
                    rel_tol=1e-15,
                    abs_tol=1e-9,
                )
            ):
                raise ValueError(
                    f"{key}: disk-geometry mass systematic does not match "
                    "the reviewed source values"
                )
        else:  # Defensive: the status inventory above is exhaustive.
            raise AssertionError(
                f"unhandled mass systematic status: {mass_systematic_status}"
            )
        if mass_display_status == "decomposed_components":
            total_minus = optional_float(param["mass_total_unc_minus_msun"])
            total_plus = optional_float(param["mass_total_unc_plus_msun"])
            if (
                display_mass_minus is None
                or display_mass_plus is None
                or total_minus is None
                or total_plus is None
                or not math.isclose(
                    display_mass_minus, total_minus, rel_tol=0.0, abs_tol=1e-9
                )
                or not math.isclose(
                    display_mass_plus, total_plus, rel_tol=0.0, abs_tol=1e-9
                )
            ):
                raise ValueError(
                    f"{key}: decomposed mass display must match diagnostic total"
                )
        elif (
            mass_display_status
            == "separate_decomposed_and_position_systematic_components"
        ):
            if (
                key != "UGC6093"
                or mass_systematic_status
                != "source_reported_absolute_position_systematic"
                or display_mass_minus is None
                or display_mass_plus is None
                or not math.isclose(
                    display_mass_minus,
                    1_429_512.309715127,
                    rel_tol=1e-15,
                    abs_tol=1e-9,
                )
                or not math.isclose(
                    display_mass_plus,
                    1_549_962.214996657,
                    rel_tol=1e-15,
                    abs_tol=1e-9,
                )
            ):
                raise ValueError(
                    f"{key}: separate position-systematic mass components do "
                    "not match reviewed bounds"
                )
        elif (
            mass_display_status
            == "separate_decomposed_and_geometry_systematic_components"
        ):
            if (
                key != "ESO558-G009"
                or mass_systematic_status
                != "source_reported_midline_inclination_systematic"
                or display_mass_minus is None
                or display_mass_plus is None
                or not math.isclose(
                    display_mass_minus,
                    950_221.0191704632,
                    rel_tol=1e-15,
                    abs_tol=1e-9,
                )
                or not math.isclose(
                    display_mass_plus,
                    950_221.0191704632,
                    rel_tol=1e-15,
                    abs_tol=1e-9,
                )
            ):
                raise ValueError(
                    f"{key}: separate geometry-systematic mass components do "
                    "not match reviewed bounds"
                )
        elif mass_display_status == "separate_statistical_and_systematic_components":
            if (
                key != "CGCG074-064"
                or mass_systematic_status
                != "source_reported_model_selection_systematic"
                or display_mass_minus != 2_000_000.0
                or display_mass_plus != 2_200_000.0
            ):
                raise ValueError(
                    f"{key}: separate mass components do not match reviewed bounds"
                )
        elif mass_display_status == "source_reported_total":
            if (
                key != "Circinus"
                or display_mass_minus != 300_000.0
                or display_mass_plus != 300_000.0
            ):
                raise ValueError(
                    f"{key}: source-reported mass display does not match reviewed bounds"
                )
        elif display_mass_minus is not None or display_mass_plus is not None:
            raise ValueError(f"{key}: unavailable mass display must be blank")

        distance_component_status = param["distance_stat_component_status"]
        if distance_component_status not in DISTANCE_STAT_LABELS:
            raise ValueError(f"{key}: unsupported distance component status")
        distance_stat_minus = optional_float(param["distance_stat_unc_minus_mpc"])
        distance_stat_plus = optional_float(param["distance_stat_unc_plus_mpc"])
        if (distance_stat_minus is None) != (distance_stat_plus is None):
            raise ValueError(f"{key}: incomplete distance interval")
        if distance_component_status == "not_available":
            if distance_stat_minus is not None or distance_stat_plus is not None:
                raise ValueError(f"{key}: unavailable distance component must be blank")
        elif distance_stat_minus is None or distance_stat_plus is None:
            raise ValueError(f"{key}: distance component requires numerical bounds")
        if key == "Circinus":
            if mass_display_status != "source_reported_total":
                raise ValueError("Circinus: published mass total must control display")
            if distance_component_status != "reported_legacy_multimethod":
                raise ValueError("Circinus: distance component must be labeled reported")
        elif key == "NGC3079":
            if distance_component_status != "reported_legacy_tully_fisher":
                raise ValueError(
                    "NGC3079: Tully--Fisher distance component must be labeled reported"
                )
        elif mass_display_status == "source_reported_total":
            raise ValueError(f"{key}: unexpected source-reported mass display")
        elif distance_component_status in {
            "reported_legacy_multimethod",
            "reported_legacy_tully_fisher",
        }:
            raise ValueError(f"{key}: unexpected reported legacy distance status")

        if key == "NGC4388":
            distance_sys_minus = optional_float(
                param["distance_systematic_unc_minus_mpc"]
            )
            distance_sys_plus = optional_float(
                param["distance_systematic_unc_plus_mpc"]
            )
            if distance_component_status != "not_available":
                raise ValueError(
                    "NGC4388: Tully--Fisher distance term must not be labeled statistical"
                )
            if distance_sys_minus != 2.09 or distance_sys_plus != 2.09:
                raise ValueError(
                    "NGC4388: Tully--Fisher Cepheid-calibration systematic must be 2.09 Mpc"
                )
            if (
                param["size_stat_component_status"]
                != "fit_plus_absolute_position_mixed_component"
            ):
                raise ValueError(
                    "NGC4388: size interval must be labeled fit plus position"
                )
            if meta["inclination_tex"] != r"\mathrm{unconstr.}":
                raise ValueError(
                    "NGC4388: inclination must remain unconstrained, not measured edge-on"
                )

        if key == "IC2560":
            reviewed_metadata = {
                "mass_reference_key": "Nakai2024",
                "distance_reference_key": "Nakai2024",
                "inclination_reference_key": "Yamauchi2012",
                "inclination_tex": r"\mathrm{e.o.}",
            }
            for field, expected in reviewed_metadata.items():
                if meta[field] != expected:
                    raise ValueError(
                        f"IC2560: {field} does not match the reviewed evidence"
                    )

        if key == "NGC3079":
            reviewed_metadata = {
                "source_tex": "NGC~3079",
                "mass_reference_key": "Kondratko2005",
                "distance_reference_key": "Kondratko2005",
                "inclination_reference_key": "Kondratko2005",
                "inclination_tex": r"\mathrm{t/s}",
            }
            for field, expected in reviewed_metadata.items():
                if meta[field] != expected:
                    raise ValueError(
                        f"NGC3079: {field} does not match the reviewed evidence"
                    )

        if key == "NGC5765b":
            reviewed_metadata = {
                "mass_reference_key": "Pesce2020b",
                "distance_reference_key": "Pesce2020b",
                "inclination_reference_key": "Pesce2020b",
                "inclination_tex": r"(0):72.4\pm0.5",
            }
            for field, expected in reviewed_metadata.items():
                if meta[field] != expected:
                    raise ValueError(
                        f"NGC5765b: {field} does not match the reviewed updated fit"
                    )

        if key == "NGC6264":
            reviewed_metadata = {
                "mass_reference_key": "Pesce2020b",
                "distance_reference_key": "Pesce2020b",
                "inclination_reference_key": "Pesce2020b",
                "inclination_tex": r"(0):91.3\pm2.3",
            }
            for field, expected in reviewed_metadata.items():
                if meta[field] != expected:
                    raise ValueError(
                        f"NGC6264: {field} does not match the reviewed updated fit"
                    )

        if key == "UGC3789":
            reviewed_metadata = {
                "mass_reference_key": "Pesce2020b",
                "distance_reference_key": "Pesce2020b",
                "inclination_reference_key": "Pesce2020b",
                "inclination_tex": r"(0):84.9\pm0.6",
            }
            for field, expected in reviewed_metadata.items():
                if meta[field] != expected:
                    raise ValueError(
                        f"UGC3789: {field} does not match the reviewed updated fit"
                    )

        if key == "NGC6323":
            reviewed_metadata = {
                "source_tex": "NGC~6323",
                "mass_reference_key": "Pesce2020b",
                "distance_reference_key": "Pesce2020b",
                "inclination_reference_key": "Pesce2020b",
                "inclination_tex": r"(0):91.5\pm0.3",
            }
            for field, expected in reviewed_metadata.items():
                if meta[field] != expected:
                    raise ValueError(
                        f"NGC6323: {field} does not match the reviewed updated fit"
                    )

        if key == "J0437":
            reviewed_values = {
                "size_uas": 0.004191067803754549,
                "size_stat_unc_minus_uas": 0.0002930816645982202,
                "size_stat_unc_plus_uas": 0.0002930816645982202,
                "mass_msun": 2_635_479.2190218437,
                "mass_fit_unc_minus_msun": 184_299.24608544362,
                "mass_fit_unc_plus_msun": 184_299.24608544362,
                "mass_shared_distance_unc_minus_msun": 176_180.85157776746,
                "mass_shared_distance_unc_plus_msun": 176_180.85157776746,
                "mass_display_unc_minus_msun": 254_962.55523180307,
                "mass_display_unc_plus_msun": 254_962.55523180307,
                "distance_stat_unc_minus_mpc": 3.4259780517380682,
                "distance_stat_unc_plus_mpc": 3.4259780517380682,
                "distance_unc_minus_mpc": 4.312118744910393,
                "distance_unc_plus_mpc": 4.312118744910393,
                "size_model_check_accepted_configuration_count": 1.0,
                "size_model_check_lower_uas": J0437_MODEL_CHECK_LOWER_UAS,
                "size_model_check_upper_uas": J0437_MODEL_CHECK_UPPER_UAS,
                "nonadopted_fit_ratio_msun_per_mpc": 42_500.0,
                "nonadopted_fit_ratio_unc_minus_msun_per_mpc": 1_200.0,
                "nonadopted_fit_ratio_unc_plus_msun_per_mpc": 1_000.0,
                "nonadopted_fit_size_uas": 0.004359589760898525,
                "nonadopted_fit_size_unc_minus_uas": 0.00012309429913125248,
                "nonadopted_fit_size_unc_plus_uas": 0.00010257858260937706,
            }
            for field, expected in reviewed_values.items():
                actual = optional_float(param[field])
                if actual is None or not math.isclose(
                    actual, expected, rel_tol=0.0, abs_tol=1e-15
                ):
                    raise ValueError(
                        f"J0437: {field} does not match the reviewed value"
                    )
            reviewed_statuses = {
                "treatment": "shared_distance",
                "mass_display_interval_status": "decomposed_components",
                "distance_stat_component_status": "statistical",
                "size_stat_component_status": "fixed_distance_formal_fit_interval",
                "size_systematic_status": "unquantified_model",
                "size_total_status": "incomplete_unquantified_systematic",
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
                "nonadopted_fit_reference_key": "Villaraos2022",
                "nonadopted_fit_reference_doi": "10.1093/mnras/stac2973",
            }
            for field, expected in reviewed_statuses.items():
                if param[field] != expected:
                    raise ValueError(
                        f"J0437: {field} does not match the reviewed value"
                    )
            for field in (
                "size_systematic_unc_minus_uas",
                "size_systematic_unc_plus_uas",
                "size_total_unc_minus_uas",
                "size_total_unc_plus_uas",
            ):
                if optional_float(param[field]) is not None:
                    raise ValueError(f"J0437: {field} must remain NA")
            if param["size_model_check_excluded_configuration_ids"] != "NA":
                raise ValueError(
                    "J0437: size_model_check_excluded_configuration_ids must remain NA"
                )
            reviewed_metadata = {
                "source_tex": "J0437+2456",
                "mass_reference_key": "Pesce2021J0437",
                "distance_reference_key": "Gao2016b",
                "inclination_reference_key": "Pesce2021J0437",
                "inclination_tex": r"\mathrm{unconstr.}",
            }
            for field, expected in reviewed_metadata.items():
                if meta[field] != expected:
                    raise ValueError(
                        f"J0437: {field} does not match the reviewed updated fit"
                    )

        if key == "Mrk1029":
            reviewed_values = {
                "distance_mpc": 119.32882273342354,
                "distance_stat_unc_minus_mpc": 3.4105504657071966,
                "distance_stat_unc_plus_mpc": 3.4105504657071966,
                "distance_systematic_unc_minus_mpc": 4.844201193508398,
                "distance_systematic_unc_plus_mpc": 4.844201193508398,
                "distance_unc_minus_mpc": 5.924368293946939,
                "distance_unc_plus_mpc": 5.924368293946939,
                "size_model_check_accepted_configuration_count": 1.0,
                "size_model_check_lower_uas": MRK1029_MODEL_CHECK_LOWER_UAS,
                "size_model_check_upper_uas": MRK1029_MODEL_CHECK_UPPER_UAS,
                "nonadopted_fit_ratio_msun_per_mpc": 14_400.0,
                "nonadopted_fit_ratio_unc_minus_msun_per_mpc": 1_100.0,
                "nonadopted_fit_ratio_unc_plus_msun_per_mpc": 1_100.0,
                "nonadopted_fit_size_uas": 0.0014771315895750299,
                "nonadopted_fit_size_unc_minus_uas": 0.00011283644087031476,
                "nonadopted_fit_size_unc_plus_uas": 0.00011283644087031476,
            }
            for field, expected in reviewed_values.items():
                actual = optional_float(param[field])
                if actual is None or not math.isclose(
                    actual, expected, rel_tol=0.0, abs_tol=1e-15
                ):
                    raise ValueError(
                        f"Mrk1029: {field} does not match the reviewed value"
                    )
            reviewed_statuses = {
                "treatment": "no_robust_mass",
                "mass_display_interval_status": "not_available",
                "distance_stat_component_status": "statistical",
                "size_stat_component_status": "not_available",
                "size_systematic_status": "not_applicable",
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
            for field, expected in reviewed_statuses.items():
                if param[field] != expected:
                    raise ValueError(
                        f"Mrk1029: {field} does not match the reviewed value"
                    )
            for field in (
                "size_uas",
                "mass_msun",
                "mass_total_unc_minus_msun",
                "mass_total_unc_plus_msun",
                "mass_display_unc_minus_msun",
                "mass_display_unc_plus_msun",
                "size_stat_unc_minus_uas",
                "size_stat_unc_plus_uas",
                "size_systematic_unc_minus_uas",
                "size_systematic_unc_plus_uas",
                "size_model_check_excluded_configuration_ids",
            ):
                if optional_float(param[field]) is not None:
                    raise ValueError(f"Mrk1029: {field} must remain NA")
            reviewed_metadata = {
                "source_tex": "Mrk~1029",
                "mass_reference_key": "Gao2016b",
                "distance_reference_key": "Gao2016b",
                "inclination_reference_key": "Gao2016b",
                "inclination_tex": r"79\pm2\;\mathrm{(est.)}",
            }
            for field, expected in reviewed_metadata.items():
                if meta[field] != expected:
                    raise ValueError(
                        f"Mrk1029: {field} does not match the reviewed evidence"
                    )

        if key == "NGC1320":
            reviewed_values = {
                "distance_mpc": 33.783491204330176,
                "distance_stat_unc_minus_mpc": 3.389871133734789,
                "distance_stat_unc_plus_mpc": 3.389871133734789,
                "distance_systematic_unc_minus_mpc": 1.371454311407179,
                "distance_systematic_unc_plus_mpc": 1.371454311407179,
                "distance_unc_minus_mpc": 3.656790017434105,
                "distance_unc_plus_mpc": 3.656790017434105,
                "size_model_check_accepted_configuration_count": 1.0,
                "size_model_check_lower_uas": NGC1320_MODEL_CHECK_LOWER_UAS,
                "size_model_check_upper_uas": NGC1320_MODEL_CHECK_UPPER_UAS,
                "nonadopted_fit_ratio_msun_per_mpc": 149_700.0,
                "nonadopted_fit_ratio_unc_minus_msun_per_mpc": 6_900.0,
                "nonadopted_fit_ratio_unc_plus_msun_per_mpc": 7_100.0,
                "nonadopted_fit_size_uas": 0.015356013816623746,
                "nonadopted_fit_size_unc_minus_uas": 0.0007077922200047018,
                "nonadopted_fit_size_unc_plus_uas": 0.0007283079365265772,
            }
            for field, expected in reviewed_values.items():
                actual = optional_float(param[field])
                if actual is None or not math.isclose(
                    actual, expected, rel_tol=0.0, abs_tol=1e-15
                ):
                    raise ValueError(
                        f"NGC1320: {field} does not match the reviewed value"
                    )
            reviewed_statuses = {
                "treatment": "no_robust_mass",
                "mass_display_interval_status": "not_available",
                "distance_stat_component_status": "statistical",
                "size_stat_component_status": "not_available",
                "size_systematic_status": "not_applicable",
                "size_total_status": "not_available",
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
            for field, expected in reviewed_statuses.items():
                if param[field] != expected:
                    raise ValueError(
                        f"NGC1320: {field} does not match the reviewed value"
                    )
            for field in (
                "size_uas",
                "mass_msun",
                "mass_total_unc_minus_msun",
                "mass_total_unc_plus_msun",
                "mass_display_unc_minus_msun",
                "mass_display_unc_plus_msun",
                "size_stat_unc_minus_uas",
                "size_stat_unc_plus_uas",
                "size_systematic_unc_minus_uas",
                "size_systematic_unc_plus_uas",
                "size_model_check_excluded_configuration_ids",
            ):
                if optional_float(param[field]) is not None:
                    raise ValueError(f"NGC1320: {field} must remain NA")
            reviewed_metadata = {
                "source_tex": "NGC~1320",
                "mass_reference_key": "Gao2016b",
                "distance_reference_key": "Gao2016b",
                "inclination_reference_key": "Gao2016b",
                "inclination_tex": r"\mathrm{unconstr.}",
            }
            for field, expected in reviewed_metadata.items():
                if meta[field] != expected:
                    raise ValueError(
                        f"NGC1320: {field} does not match the reviewed evidence"
                    )

        if key == "ESO558-G009":
            reviewed_metadata = {
                "mass_reference_key": "Gao2016b",
                "distance_reference_key": "Gao2016b",
                "inclination_reference_key": "Gao2016b",
                "inclination_tex": r"98\pm1",
            }
            for field, expected in reviewed_metadata.items():
                if meta[field] != expected:
                    raise ValueError(
                        f"ESO558-G009: {field} does not match the reviewed fit"
                    )

        if key == "CGCG074-064":
            reviewed_metadata = {
                "mass_reference_key": "Pesce2020a",
                "distance_reference_key": "Pesce2020a",
                "inclination_reference_key": "Pesce2020a",
                "inclination_tex": (
                    r"(0):90.8\pm0.6\,\mathrm{(stat.)}\;\pm1.1\,\mathrm{(sys.)}"
                ),
            }
            for field, expected in reviewed_metadata.items():
                if meta[field] != expected:
                    raise ValueError(
                        f"CGCG074-064: {field} does not match the reviewed fit"
                    )

        mass_key = meta["mass_reference_key"]
        distance_key = meta["distance_reference_key"]
        if mass_key not in REFERENCE_DOIS or distance_key not in REFERENCE_DOIS:
            raise ValueError(f"{key}: mass/distance BibTeX key lacks a DOI mapping")
        if param["parameter_reference_doi"] != REFERENCE_DOIS[mass_key]:
            raise ValueError(f"{key}: mass reference key does not match parameter DOI")
        if param["distance_reference_doi"] != REFERENCE_DOIS[distance_key]:
            raise ValueError(f"{key}: distance reference key does not match distance DOI")
        if not meta["inclination_reference_key"]:
            raise ValueError(f"{key}: inclination reference key is required")
        if key in {"NGC1194", "NGC2273"}:
            if param["parameter_reference_doi"] != REFERENCE_DOIS["Kuo2010"]:
                raise ValueError(f"{key}: central mass normalization must cite Kuo 2011")
            if param["mass_uncertainty_reference_doi"] != REFERENCE_DOIS["Gao2016b"]:
                raise ValueError(f"{key}: reported uncertainty must cite Gao 2017")

        model_key = param["size_model_check_key"]
        model_fields = (
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
        )
        if key == "NGC4258":
            if model_key != "van_der_boom_kushnir_2026_component_tracking":
                raise ValueError("NGC4258: required data-selection check is missing")
            if param["size_model_check_type"] != "data_selection_component_tracking":
                raise ValueError("NGC4258: unsupported data-selection check type")
            if (
                param["size_model_check_data_relation"]
                != "subset_reanalysis_of_argon2007_vlbi_epochs_used_by_primary"
            ):
                raise ValueError("NGC4258: unsupported model-check data relation")
            if param["size_model_check_range_basis"] != "accepted_configuration_medians":
                raise ValueError("NGC4258: unsupported model-check range basis")
            if param["size_model_check_accepted_configuration_count"] != "15":
                raise ValueError("NGC4258: accepted model-check count must be 15")
            if param["size_model_check_excluded_configuration_ids"] != "2;17":
                raise ValueError("NGC4258: excluded model-check configurations must be 2;17")
            if param["size_model_check_probability_status"] != "non_probabilistic":
                raise ValueError("NGC4258: data-selection range must be non-probabilistic")
            if (
                param["size_model_check_adoption_status"]
                != "reported_separately_not_adopted"
            ):
                raise ValueError("NGC4258: data-selection range must remain separate")
            reference_key = param["size_model_check_reference_key"]
            if (
                reference_key != "vanDerBoom2026"
                or param["size_model_check_reference_doi"]
                != REFERENCE_DOIS[reference_key]
            ):
                raise ValueError("NGC4258: model-check citation does not match its DOI")
            lower = optional_float(param["size_model_check_lower_uas"])
            upper = optional_float(param["size_model_check_upper_uas"])
            if lower is None or upper is None or not 0.0 < lower < upper:
                raise ValueError("NGC4258: invalid data-selection range")
            if not math.isclose(
                lower,
                NGC4258_MODEL_CHECK_LOWER_UAS,
                rel_tol=0.0,
                abs_tol=1e-15,
            ) or not math.isclose(
                upper,
                NGC4258_MODEL_CHECK_UPPER_UAS,
                rel_tol=0.0,
                abs_tol=1e-15,
            ):
                raise ValueError("NGC4258: model-check bounds do not match reviewed endpoints")
        elif key == "NGC3079":
            reviewed_numeric = {
                "size_uas": 0.011858795677384631,
                "mass_msun": 2_000_000.0,
                "distance_mpc": 17.3,
                "distance_stat_unc_minus_mpc": 1.5,
                "distance_stat_unc_plus_mpc": 1.5,
                "distance_unc_minus_mpc": 1.5,
                "distance_unc_plus_mpc": 1.5,
                "mass_shared_distance_unc_minus_msun": 173_410.40462427746,
                "mass_shared_distance_unc_plus_msun": 173_410.40462427746,
            }
            for field, expected in reviewed_numeric.items():
                actual = optional_float(param[field])
                if actual is None or not math.isclose(
                    actual, expected, rel_tol=1e-15, abs_tol=1e-15
                ):
                    raise ValueError(
                        f"NGC3079: {field} does not match the reviewed evidence"
                    )
            reviewed_statuses = {
                "size_stat_component_status": "no_formal_uncertainty",
                "size_systematic_status": "unquantified_model",
                "size_total_status": "incomplete_unquantified_systematic",
                "distance_stat_component_status": "reported_legacy_tully_fisher",
                "mass_display_interval_status": "not_available",
            }
            for field, expected in reviewed_statuses.items():
                if param[field] != expected:
                    raise ValueError(
                        f"NGC3079: {field} does not match the reviewed evidence"
                    )
            for field in (
                "size_unc_adopted_uas",
                "size_stat_unc_minus_uas",
                "size_stat_unc_plus_uas",
                "size_systematic_unc_minus_uas",
                "size_systematic_unc_plus_uas",
                "size_total_unc_minus_uas",
                "size_total_unc_plus_uas",
                "mass_total_unc_minus_msun",
                "mass_total_unc_plus_msun",
                "mass_fit_unc_minus_msun",
                "mass_fit_unc_plus_msun",
                "mass_systematic_unc_minus_msun",
                "mass_systematic_unc_plus_msun",
                "mass_display_unc_minus_msun",
                "mass_display_unc_plus_msun",
            ):
                if optional_float(param[field]) is not None:
                    raise ValueError(f"NGC3079: {field} must remain unavailable")
            if model_key != "kuo_2020_ngc3079_turbulent_support_bias":
                raise ValueError("NGC3079: required turbulent-support check is missing")
            if param["size_model_check_type"] != "turbulent_support_mass_bias":
                raise ValueError("NGC3079: unsupported turbulent-support check type")
            if (
                param["size_model_check_data_relation"]
                != "published_analytic_model_applied_to_ngc3079_disk_parameters"
            ):
                raise ValueError("NGC3079: unsupported model-check data relation")
            if (
                param["size_model_check_range_basis"]
                != "one_sided_fractional_mass_underestimate"
            ):
                raise ValueError("NGC3079: unsupported model-check range basis")
            if param["size_model_check_accepted_configuration_count"] != "NA":
                raise ValueError("NGC3079: configuration count must remain NA")
            if param["size_model_check_excluded_configuration_ids"] != "NA":
                raise ValueError("NGC3079: excluded configurations must remain NA")
            if param["size_model_check_probability_status"] != "non_probabilistic":
                raise ValueError("NGC3079: model-sensitivity range must be non-probabilistic")
            if (
                param["size_model_check_adoption_status"]
                != "reported_separately_not_adopted"
            ):
                raise ValueError("NGC3079: model-sensitivity range must remain separate")
            reference_key = param["size_model_check_reference_key"]
            if (
                reference_key != "Kuo2020"
                or param["size_model_check_reference_doi"]
                != REFERENCE_DOIS[reference_key]
            ):
                raise ValueError("NGC3079: model-check citation does not match its DOI")
            lower = optional_float(param["size_model_check_lower_uas"])
            upper = optional_float(param["size_model_check_upper_uas"])
            if lower is None or upper is None or not 0.0 < lower < upper:
                raise ValueError("NGC3079: invalid one-sided model-sensitivity range")
            if not math.isclose(
                lower,
                NGC3079_MODEL_CHECK_LOWER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ) or not math.isclose(
                upper,
                NGC3079_MODEL_CHECK_UPPER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ):
                raise ValueError("NGC3079: model-check bounds do not match reviewed endpoints")
        elif key == "NGC4945":
            if model_key != "kuo_2020_ngc4945_turbulent_support_bias":
                raise ValueError("NGC4945: required turbulent-support check is missing")
            if param["size_model_check_type"] != "turbulent_support_mass_bias":
                raise ValueError("NGC4945: unsupported turbulent-support check type")
            if (
                param["size_model_check_data_relation"]
                != "published_analytic_model_applied_to_ngc4945_disk_parameters"
            ):
                raise ValueError("NGC4945: unsupported model-check data relation")
            if (
                param["size_model_check_range_basis"]
                != "one_sided_fractional_mass_underestimate"
            ):
                raise ValueError("NGC4945: unsupported model-check range basis")
            if param["size_model_check_accepted_configuration_count"] != "NA":
                raise ValueError("NGC4945: configuration count must remain NA")
            if param["size_model_check_excluded_configuration_ids"] != "NA":
                raise ValueError("NGC4945: excluded configurations must remain NA")
            if param["size_model_check_probability_status"] != "non_probabilistic":
                raise ValueError("NGC4945: model-sensitivity range must be non-probabilistic")
            if (
                param["size_model_check_adoption_status"]
                != "reported_separately_not_adopted"
            ):
                raise ValueError("NGC4945: model-sensitivity range must remain separate")
            reference_key = param["size_model_check_reference_key"]
            if (
                reference_key != "Kuo2020"
                or param["size_model_check_reference_doi"]
                != REFERENCE_DOIS[reference_key]
            ):
                raise ValueError("NGC4945: model-check citation does not match its DOI")
            lower = optional_float(param["size_model_check_lower_uas"])
            upper = optional_float(param["size_model_check_upper_uas"])
            if lower is None or upper is None or not 0.0 < lower < upper:
                raise ValueError("NGC4945: invalid one-sided model-sensitivity range")
            if not math.isclose(
                lower,
                NGC4945_MODEL_CHECK_LOWER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ) or not math.isclose(
                upper,
                NGC4945_MODEL_CHECK_UPPER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ):
                raise ValueError("NGC4945: model-check bounds do not match reviewed endpoints")
        elif key == "NGC4388":
            if model_key != "kuo_2020_ngc4388_turbulent_support_bias":
                raise ValueError("NGC4388: required turbulent-support check is missing")
            if param["size_model_check_type"] != "turbulent_support_mass_bias":
                raise ValueError("NGC4388: unsupported turbulent-support check type")
            if (
                param["size_model_check_data_relation"]
                != "published_analytic_model_applied_to_ngc4388_disk_parameters"
            ):
                raise ValueError("NGC4388: unsupported model-check data relation")
            if (
                param["size_model_check_range_basis"]
                != "one_sided_fractional_mass_underestimate"
            ):
                raise ValueError("NGC4388: unsupported model-check range basis")
            if param["size_model_check_accepted_configuration_count"] != "NA":
                raise ValueError("NGC4388: configuration count must remain NA")
            if param["size_model_check_excluded_configuration_ids"] != "NA":
                raise ValueError("NGC4388: excluded configurations must remain NA")
            if param["size_model_check_probability_status"] != "non_probabilistic":
                raise ValueError(
                    "NGC4388: model-sensitivity range must be non-probabilistic"
                )
            if (
                param["size_model_check_adoption_status"]
                != "reported_separately_not_adopted"
            ):
                raise ValueError("NGC4388: model-sensitivity range must remain separate")
            reference_key = param["size_model_check_reference_key"]
            if (
                reference_key != "Kuo2020"
                or param["size_model_check_reference_doi"]
                != REFERENCE_DOIS[reference_key]
            ):
                raise ValueError("NGC4388: model-check citation does not match its DOI")
            lower = optional_float(param["size_model_check_lower_uas"])
            upper = optional_float(param["size_model_check_upper_uas"])
            if lower is None or upper is None or not 0.0 < lower < upper:
                raise ValueError("NGC4388: invalid one-sided model-sensitivity range")
            if not math.isclose(
                lower,
                NGC4388_MODEL_CHECK_LOWER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ) or not math.isclose(
                upper,
                NGC4388_MODEL_CHECK_UPPER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ):
                raise ValueError(
                    "NGC4388: model-check bounds do not match reviewed endpoints"
                )
        elif key == "NGC6323":
            if model_key != "kuo_2018_ngc6323_physically_allowed_disk_mass_models":
                raise ValueError("NGC6323: required disk-model check is missing")
            if param["size_model_check_type"] != "disk_mass_model_choice":
                raise ValueError("NGC6323: unsupported disk-model check type")
            if (
                param["size_model_check_data_relation"]
                != "predecessor_fit_data_alternative_accretion_disk_models"
            ):
                raise ValueError("NGC6323: unsupported model-check data relation")
            if (
                param["size_model_check_range_basis"]
                != "published_physically_allowed_model_envelope"
            ):
                raise ValueError("NGC6323: unsupported model-check range basis")
            if param["size_model_check_accepted_configuration_count"] != "3":
                raise ValueError("NGC6323: accepted model count must remain three")
            if param["size_model_check_excluded_configuration_ids"] != "2;4":
                raise ValueError("NGC6323: excluded unphysical models must remain 2;4")
            if param["size_model_check_probability_status"] != "non_probabilistic":
                raise ValueError("NGC6323: model envelope must be non-probabilistic")
            if (
                param["size_model_check_adoption_status"]
                != "reported_separately_not_adopted"
            ):
                raise ValueError("NGC6323: model envelope must remain separate")
            reference_key = param["size_model_check_reference_key"]
            if (
                reference_key != "Kuo2018Disk"
                or param["size_model_check_reference_doi"]
                != REFERENCE_DOIS[reference_key]
            ):
                raise ValueError("NGC6323: model-check citation does not match its DOI")
            lower = optional_float(param["size_model_check_lower_uas"])
            upper = optional_float(param["size_model_check_upper_uas"])
            if lower is None or upper is None or not 0.0 < lower < upper:
                raise ValueError("NGC6323: invalid one-sided model-sensitivity range")
            if not math.isclose(
                lower,
                NGC6323_MODEL_CHECK_LOWER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ) or not math.isclose(
                upper,
                NGC6323_MODEL_CHECK_UPPER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ):
                raise ValueError(
                    "NGC6323: model-check bounds do not match reviewed endpoints"
                )
            reviewed_statuses = {
                "treatment": "paired_posterior",
                "size_stat_component_status": "paired_chain_quantile_interval",
                "size_systematic_status": (
                    "numeric_one_percent_non_keplerian_allowance"
                ),
                "size_total_status": "quadrature_statistical_plus_systematic",
            }
            for field, expected in reviewed_statuses.items():
                if param[field] != expected:
                    raise ValueError(
                        f"NGC6323: {field} does not match the reviewed primary"
                    )
            reviewed_numeric = {
                "size_uas": 0.00954948892245025,
                "size_stat_unc_minus_uas": 5.611844616523233e-05,
                "size_stat_unc_plus_uas": 6.131459858793598e-05,
                "size_systematic_unc_minus_uas": 9.54948892245025e-05,
                "size_systematic_unc_plus_uas": 9.54948892245025e-05,
                "mass_msun": 10_200_000.0,
                "distance_mpc": 109.4,
            }
            for field, expected in reviewed_numeric.items():
                actual = optional_float(param[field])
                if actual is None or not math.isclose(
                    actual, expected, rel_tol=0.0, abs_tol=1e-15
                ):
                    raise ValueError(
                        f"NGC6323: {field} does not match the reviewed primary"
                    )
        elif key == "J0437":
            if model_key != "deepshikha_2024_j0437_mestel_disk_mass_model":
                raise ValueError("J0437: required Mestel-disk check is missing")
            if param["size_model_check_type"] != "disk_mass_model_choice":
                raise ValueError("J0437: unsupported disk-model check type")
            if (
                param["size_model_check_data_relation"]
                != "reanalysis_of_published_gao2017_vlbi_data_data_selection_unreported"
            ):
                raise ValueError("J0437: unsupported model-check data relation")
            if (
                param["size_model_check_range_basis"]
                != "reported_symmetric_mass_error_span"
            ):
                raise ValueError("J0437: unsupported model-check range basis")
            if param["size_model_check_accepted_configuration_count"] != "1":
                raise ValueError("J0437: model count must remain one")
            if param["size_model_check_excluded_configuration_ids"] != "NA":
                raise ValueError("J0437: excluded configurations must remain NA")
            if param["size_model_check_probability_status"] != "non_probabilistic":
                raise ValueError("J0437: model-sensitivity range must be non-probabilistic")
            if (
                param["size_model_check_adoption_status"]
                != "reported_separately_not_adopted"
            ):
                raise ValueError("J0437: model-sensitivity range must remain separate")
            reference_key = param["size_model_check_reference_key"]
            if (
                reference_key != "Deepshikha2024"
                or param["size_model_check_reference_doi"]
                != REFERENCE_DOIS[reference_key]
            ):
                raise ValueError("J0437: model-check citation does not match its DOI")
            lower = optional_float(param["size_model_check_lower_uas"])
            upper = optional_float(param["size_model_check_upper_uas"])
            if lower is None or upper is None or not 0.0 < lower < upper:
                raise ValueError("J0437: invalid one-sided model-sensitivity range")
            if not math.isclose(
                lower,
                J0437_MODEL_CHECK_LOWER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ) or not math.isclose(
                upper,
                J0437_MODEL_CHECK_UPPER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ):
                raise ValueError(
                    "J0437: model-check bounds do not match reviewed endpoints"
                )
        elif key == "Mrk1029":
            if model_key != "deepshikha_2024_mrk1029_mestel_disk_mass_model":
                raise ValueError("Mrk1029: required Mestel-disk diagnostic is missing")
            if param["size_model_check_type"] != "disk_mass_model_choice":
                raise ValueError("Mrk1029: unsupported disk-model diagnostic type")
            if (
                param["size_model_check_data_relation"]
                != "reanalysis_of_published_gao2017_vlbi_data_data_selection_unreported"
            ):
                raise ValueError("Mrk1029: unsupported model-check data relation")
            if (
                param["size_model_check_range_basis"]
                != "reported_symmetric_mass_error_span"
            ):
                raise ValueError("Mrk1029: unsupported model-check range basis")
            if param["size_model_check_accepted_configuration_count"] != "1":
                raise ValueError("Mrk1029: model count must remain one")
            if param["size_model_check_excluded_configuration_ids"] != "NA":
                raise ValueError("Mrk1029: excluded configurations must remain NA")
            if param["size_model_check_probability_status"] != "non_probabilistic":
                raise ValueError("Mrk1029: model diagnostic must be non-probabilistic")
            if (
                param["size_model_check_adoption_status"]
                != "reported_separately_not_adopted"
            ):
                raise ValueError("Mrk1029: model diagnostic must remain separate")
            reference_key = param["size_model_check_reference_key"]
            if (
                reference_key != "Deepshikha2024"
                or param["size_model_check_reference_doi"]
                != REFERENCE_DOIS[reference_key]
            ):
                raise ValueError("Mrk1029: model-check citation does not match its DOI")
            lower = optional_float(param["size_model_check_lower_uas"])
            upper = optional_float(param["size_model_check_upper_uas"])
            if lower is None or upper is None or not 0.0 < lower < upper:
                raise ValueError("Mrk1029: invalid standalone model-diagnostic range")
            if not math.isclose(
                lower,
                MRK1029_MODEL_CHECK_LOWER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ) or not math.isclose(
                upper,
                MRK1029_MODEL_CHECK_UPPER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ):
                raise ValueError(
                    "Mrk1029: model-check bounds do not match reviewed endpoints"
                )
            if optional_float(param["size_uas"]) is not None:
                raise ValueError("Mrk1029: standalone diagnostic cannot become a primary")
        elif key == "NGC1320":
            if model_key != "deepshikha_2024_ngc1320_mestel_disk_mass_model":
                raise ValueError("NGC1320: required Mestel-disk diagnostic is missing")
            if param["size_model_check_type"] != "disk_mass_model_choice":
                raise ValueError("NGC1320: unsupported disk-model diagnostic type")
            if (
                param["size_model_check_data_relation"]
                != "reanalysis_of_published_gao2017_vlbi_data_data_selection_unreported"
            ):
                raise ValueError("NGC1320: unsupported model-check data relation")
            if (
                param["size_model_check_range_basis"]
                != "reported_symmetric_mass_error_span"
            ):
                raise ValueError("NGC1320: unsupported model-check range basis")
            if param["size_model_check_accepted_configuration_count"] != "1":
                raise ValueError("NGC1320: model count must remain one")
            if param["size_model_check_excluded_configuration_ids"] != "NA":
                raise ValueError("NGC1320: excluded configurations must remain NA")
            if param["size_model_check_probability_status"] != "non_probabilistic":
                raise ValueError("NGC1320: model diagnostic must be non-probabilistic")
            if (
                param["size_model_check_adoption_status"]
                != "reported_separately_not_adopted"
            ):
                raise ValueError("NGC1320: model diagnostic must remain separate")
            reference_key = param["size_model_check_reference_key"]
            if (
                reference_key != "Deepshikha2024"
                or param["size_model_check_reference_doi"]
                != REFERENCE_DOIS[reference_key]
            ):
                raise ValueError("NGC1320: model-check citation does not match its DOI")
            lower = optional_float(param["size_model_check_lower_uas"])
            upper = optional_float(param["size_model_check_upper_uas"])
            if lower is None or upper is None or not 0.0 < lower < upper:
                raise ValueError("NGC1320: invalid standalone model-diagnostic range")
            if not math.isclose(
                lower,
                NGC1320_MODEL_CHECK_LOWER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ) or not math.isclose(
                upper,
                NGC1320_MODEL_CHECK_UPPER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ):
                raise ValueError(
                    "NGC1320: model-check bounds do not match reviewed endpoints"
                )
            if optional_float(param["size_uas"]) is not None:
                raise ValueError("NGC1320: standalone diagnostic cannot become a primary")
        elif key == "NGC2960":
            if model_key != "kuo_2018_ngc2960_physically_allowed_disk_mass_models":
                raise ValueError("NGC2960: required disk-model check is missing")
            if param["size_model_check_type"] != "disk_mass_model_choice":
                raise ValueError("NGC2960: unsupported disk-model check type")
            if (
                param["size_model_check_data_relation"]
                != "same_fit_data_alternative_accretion_disk_models"
            ):
                raise ValueError("NGC2960: unsupported model-check data relation")
            if (
                param["size_model_check_range_basis"]
                != "published_physically_allowed_model_envelope"
            ):
                raise ValueError("NGC2960: unsupported model-check range basis")
            if param["size_model_check_accepted_configuration_count"] != "3":
                raise ValueError("NGC2960: accepted model count must remain three")
            if param["size_model_check_excluded_configuration_ids"] != "2;4":
                raise ValueError("NGC2960: excluded unphysical models must remain 2;4")
            if param["size_model_check_probability_status"] != "non_probabilistic":
                raise ValueError("NGC2960: model envelope must be non-probabilistic")
            if (
                param["size_model_check_adoption_status"]
                != "reported_separately_not_adopted"
            ):
                raise ValueError("NGC2960: model envelope must remain separate")
            reference_key = param["size_model_check_reference_key"]
            if (
                reference_key != "Kuo2018Disk"
                or param["size_model_check_reference_doi"]
                != REFERENCE_DOIS[reference_key]
            ):
                raise ValueError("NGC2960: model-check citation does not match its DOI")
            lower = optional_float(param["size_model_check_lower_uas"])
            upper = optional_float(param["size_model_check_upper_uas"])
            if lower is None or upper is None or not 0.0 < lower < upper:
                raise ValueError("NGC2960: invalid one-sided model-sensitivity range")
            if not math.isclose(
                lower,
                NGC2960_MODEL_CHECK_LOWER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ) or not math.isclose(
                upper,
                NGC2960_MODEL_CHECK_UPPER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ):
                raise ValueError(
                    "NGC2960: model-check bounds do not match reviewed endpoints"
                )
            if param["size_systematic_status"] != "unquantified_astrometric":
                raise ValueError(
                    "NGC2960: astrometric completeness must remain unquantified"
                )
            if param["size_total_status"] != "incomplete_unquantified_systematic":
                raise ValueError("NGC2960: complete shadow total must remain unavailable")
            reviewed_metadata = {
                "mass_reference_key": "Kuo2018Disk",
                "distance_reference_key": "Kuo2010",
                "inclination_reference_key": "Kuo2010",
                "inclination_tex": r"\mathrm{e.o.}",
            }
            for field, expected in reviewed_metadata.items():
                if meta[field] != expected:
                    raise ValueError(
                        f"NGC2960: {field} does not match the reviewed evidence"
                    )
        elif key == "Circinus":
            if model_key != "kuo_2020_circinus_turbulent_support_bias":
                raise ValueError("Circinus: required turbulent-support check is missing")
            if param["size_model_check_type"] != "turbulent_support_mass_bias":
                raise ValueError("Circinus: unsupported turbulent-support check type")
            if (
                param["size_model_check_data_relation"]
                != "published_analytic_model_applied_to_circinus_disk_parameters"
            ):
                raise ValueError("Circinus: unsupported model-check data relation")
            if (
                param["size_model_check_range_basis"]
                != "one_sided_fractional_mass_underestimate"
            ):
                raise ValueError("Circinus: unsupported model-check range basis")
            if param["size_model_check_accepted_configuration_count"] != "NA":
                raise ValueError("Circinus: configuration count must remain NA")
            if param["size_model_check_excluded_configuration_ids"] != "NA":
                raise ValueError("Circinus: excluded configurations must remain NA")
            if param["size_model_check_probability_status"] != "non_probabilistic":
                raise ValueError("Circinus: model-sensitivity range must be non-probabilistic")
            if (
                param["size_model_check_adoption_status"]
                != "reported_separately_not_adopted"
            ):
                raise ValueError("Circinus: model-sensitivity range must remain separate")
            reference_key = param["size_model_check_reference_key"]
            if (
                reference_key != "Kuo2020"
                or param["size_model_check_reference_doi"]
                != REFERENCE_DOIS[reference_key]
            ):
                raise ValueError("Circinus: model-check citation does not match its DOI")
            lower = optional_float(param["size_model_check_lower_uas"])
            upper = optional_float(param["size_model_check_upper_uas"])
            if lower is None or upper is None or not 0.0 < lower < upper:
                raise ValueError("Circinus: invalid one-sided model-sensitivity range")
            if not math.isclose(
                lower,
                CIRCINUS_MODEL_CHECK_LOWER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ) or not math.isclose(
                upper,
                CIRCINUS_MODEL_CHECK_UPPER_UAS,
                rel_tol=0.0,
                abs_tol=1e-17,
            ):
                raise ValueError("Circinus: model-check bounds do not match reviewed endpoints")
        elif model_key not in {"", "NA"} or any(
            param[field] not in {"", "NA"} for field in model_fields
        ):
            raise ValueError(f"{key}: unexpected model-sensitivity check")

        nonadopted_fields = (
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
        )
        reviewed_fit = REVIEWED_NONADOPTED_FITS.get(key)
        if reviewed_fit is not None:
            expected_treatment = reviewed_fit["primary_treatment"]
            if param["treatment"] != expected_treatment:
                raise ValueError(
                    f"{key}: non-adopted fit requires {expected_treatment}"
                )
            primary_result_status = reviewed_fit["primary_result_status"]
            if primary_result_status == "absent":
                if param["mass_msun"] not in {"", "NA"} or param[
                    "size_uas"
                ] not in {"", "NA"}:
                    raise ValueError(
                        f"{key}: non-adopted fit cannot populate primary cells"
                    )
            elif primary_result_status == "retained":
                primary_mass = optional_float(param["mass_msun"])
                primary_size = optional_float(param["size_uas"])
                if (
                    primary_mass is None
                    or primary_size is None
                    or primary_mass <= 0.0
                    or primary_size <= 0.0
                ):
                    raise ValueError(
                        f"{key}: non-adopted fit requires retained primary cells"
                    )
            else:
                raise ValueError(f"{key}: unsupported non-adopted primary policy")
            if param["nonadopted_fit_key"] != reviewed_fit["fit_key"]:
                raise ValueError(f"{key}: required non-adopted fit is missing")
            reviewed_numeric = {
                "nonadopted_fit_ratio_msun_per_mpc": reviewed_fit["ratio"],
                "nonadopted_fit_ratio_unc_minus_msun_per_mpc": reviewed_fit[
                    "ratio_minus"
                ],
                "nonadopted_fit_ratio_unc_plus_msun_per_mpc": reviewed_fit[
                    "ratio_plus"
                ],
                "nonadopted_fit_size_uas": reviewed_fit["size"],
                "nonadopted_fit_size_unc_minus_uas": reviewed_fit["size_minus"],
                "nonadopted_fit_size_unc_plus_uas": reviewed_fit["size_plus"],
            }
            for field, expected in reviewed_numeric.items():
                actual = optional_float(param[field])
                if actual is None or not math.isclose(
                    actual, expected, rel_tol=0.0, abs_tol=1e-15
                ):
                    raise ValueError(
                        f"{key}: {field} does not match the reviewed non-adopted fit"
                    )
            if (
                param["nonadopted_fit_interval_status"]
                != "bayesian_formal_conditional_on_model"
            ):
                raise ValueError(f"{key}: non-adopted interval status is invalid")
            if (
                param["nonadopted_fit_adoption_status"]
                != "reported_separately_not_adopted"
            ):
                raise ValueError(f"{key}: conditional fit must remain not adopted")
            if param["nonadopted_fit_basis"] != reviewed_fit["basis"]:
                raise ValueError(f"{key}: non-adoption basis is invalid")
            reference_key = param["nonadopted_fit_reference_key"]
            if (
                reference_key != reviewed_fit["reference_key"]
                or param["nonadopted_fit_reference_doi"]
                != REFERENCE_DOIS[reference_key]
            ):
                raise ValueError(f"{key}: non-adopted fit citation does not match")
        elif param["nonadopted_fit_key"] not in {"", "NA"} or any(
            param[field] not in {"", "NA"} for field in nonadopted_fields
        ):
            raise ValueError(f"{key}: unexpected non-adopted ratio fit")


def reference_order(metadata: list[dict[str, str]]) -> list[str]:
    ordered: list[str] = []

    def add(key: str) -> None:
        if key not in ordered:
            ordered.append(key)

    for row in COMPARISON_ROWS:
        for column in ("size", "mass", "distance", "inclination"):
            add(row[f"{column}_reference_key"])
    for row in metadata:
        for column in ("mass", "distance", "inclination"):
            add(row[f"{column}_reference_key"])
    return ordered


def refmark(key: str, reference_numbers: dict[str, int]) -> str:
    return f"$^{{{reference_numbers[key]}}}$" + r"\,"


def values_close(lower: float, upper: float) -> bool:
    return math.isclose(lower, upper, rel_tol=1e-12, abs_tol=1e-15)


def symmetric_display_interval_bounds(
    lower: float, upper: float
) -> tuple[float, float]:
    """Return the arithmetic-mean half-width used in the manuscript table.

    Exact lower and upper bounds remain unchanged in ``BHS_params.csv``.  The
    table first rounds each source half-width to its displayed two-significant-
    figure value and then prints their arithmetic mean symmetrically.  Averaging
    the displayed values makes the transformation explicit and reproducible;
    for example, NGC 3393's displayed 0.010 and 0.0082 become 0.0091.
    """

    if not math.isfinite(lower) or not math.isfinite(upper):
        raise ValueError("display interval bounds must be finite")
    if lower <= 0 or upper <= 0:
        raise ValueError("display interval half-widths must be positive")
    rounded_lower = Decimal(format_significant_figures(lower, 2))
    rounded_upper = Decimal(format_significant_figures(upper, 2))
    mean_half_width = float((rounded_lower + rounded_upper) / Decimal("2"))
    return mean_half_width, mean_half_width


def quadrature_interval_bounds(
    lower_components: list[float], upper_components: list[float]
) -> tuple[float, float]:
    """Combine compatible components, then print their mean half-width."""

    if not lower_components or len(lower_components) != len(upper_components):
        raise ValueError("quadrature combination requires paired components")
    if any(value <= 0 or not math.isfinite(value) for value in lower_components):
        raise ValueError("lower quadrature components must be finite and positive")
    if any(value <= 0 or not math.isfinite(value) for value in upper_components):
        raise ValueError("upper quadrature components must be finite and positive")
    lower = math.sqrt(sum(value * value for value in lower_components))
    upper = math.sqrt(sum(value * value for value in upper_components))
    return symmetric_display_interval_bounds(lower, upper)


def decimal_places_for_significant_figures(
    value: float, significant_figures: int
) -> int:
    """Return the fixed-point decimal place needed for ``significant_figures``.

    Publication cells are rendered in fixed-point notation, including
    meaningful trailing zeros.  ``Decimal(str(value))`` avoids making the
    display precision depend on binary floating-point representation noise.
    """

    if significant_figures < 1:
        raise ValueError("significant_figures must be positive")
    decimal_value = Decimal(str(value))
    if not decimal_value.is_finite() or decimal_value <= 0:
        raise ValueError(f"publication precision requires a positive value: {value!r}")
    return significant_figures - 1 - decimal_value.adjusted()


def format_at_decimal_place(value: float, decimal_places: int) -> str:
    """Round half up and render a fixed-point value at ``decimal_places``."""

    decimal_value = Decimal(str(value))
    if not decimal_value.is_finite():
        raise ValueError(f"cannot format a non-finite value: {value!r}")
    quantum = Decimal("1").scaleb(-decimal_places)
    rounded = decimal_value.quantize(quantum, rounding=ROUND_HALF_UP)
    if rounded == 0:
        rounded = abs(rounded)
    return f"{rounded:.{max(decimal_places, 0)}f}"


def format_significant_figures(value: float, significant_figures: int) -> str:
    """Render ``value`` with an explicit number of significant figures."""

    return format_at_decimal_place(
        value,
        decimal_places_for_significant_figures(value, significant_figures),
    )


def central_decimal_places(uncertainties: list[float]) -> int:
    """Use the finest displayed two-significant-figure uncertainty place.

    Each uncertainty component keeps two significant figures independently.
    The central value follows the finest component so a smaller statistical or
    formal term remains interpretable alongside a larger systematic term.
    """

    if not uncertainties:
        raise ValueError("at least one uncertainty is required")
    return max(
        decimal_places_for_significant_figures(value, 2)
        for value in uncertainties
    )


def uncertainty_tex(
    lower: float,
    upper: float,
    *,
    label: str | None = None,
) -> str:
    suffix = rf"\,\mathrm{{({label})}}" if label else ""
    lower_text = format_significant_figures(lower, 2)
    upper_text = format_significant_figures(upper, 2)
    if values_close(lower, upper):
        return rf"\pm {lower_text}{suffix}"
    return rf"_{{-{lower_text}}}^{{+{upper_text}}}{suffix}"


def mass_log_display_bounds(param: dict[str, str]) -> tuple[float, float, float] | None:
    mass = optional_float(param["mass_msun"])
    if mass is None:
        return None
    central = math.log10(mass)
    lower = optional_float(param["mass_display_unc_minus_msun"])
    upper = optional_float(param["mass_display_unc_plus_msun"])
    if lower is None and upper is None:
        return central, math.nan, math.nan
    if lower is None or upper is None or lower >= mass:
        raise ValueError(f"{param['source_key']}: invalid total mass interval")
    return central, math.log10(mass / (mass - lower)), math.log10((mass + upper) / mass)


def format_mass_cell(
    param: dict[str, str], meta: dict[str, str], reference_numbers: dict[str, int]
) -> str:
    bounds = mass_log_display_bounds(param)
    if bounds is None:
        return r"\nodata"
    central, lower, upper = bounds
    display_uncertainty = math.isfinite(lower) and math.isfinite(upper)
    if display_uncertainty:
        lower_components = [lower]
        upper_components = [upper]
        separate_component_statuses = {
            "separate_statistical_and_systematic_components",
            "separate_decomposed_and_position_systematic_components",
            "separate_decomposed_and_geometry_systematic_components",
        }
        if param["mass_display_interval_status"] in separate_component_statuses:
            systematic_lower = optional_float(param["log_mass_systematic_unc_minus"])
            systematic_upper = optional_float(param["log_mass_systematic_unc_plus"])
            if systematic_lower is None or systematic_upper is None:
                raise ValueError(
                    f"{param['source_key']}: separate mass systematic bounds are required"
                )
            lower_components.append(systematic_lower)
            upper_components.append(systematic_upper)
        elif optional_float(param["log_mass_systematic_unc_minus"]) is not None or optional_float(
            param["log_mass_systematic_unc_plus"]
        ) is not None:
            raise ValueError(
                f"{param['source_key']}: unexpected separate mass systematic bounds"
            )

        lower, upper = quadrature_interval_bounds(
            lower_components, upper_components
        )
        value = format_at_decimal_place(
            central, central_decimal_places([lower, upper])
        )
        value += uncertainty_tex(lower, upper)
    else:
        value = format_significant_figures(central, 3)
    return refmark(meta["mass_reference_key"], reference_numbers) + f"${value}$"


def format_distance_cell(
    param: dict[str, str], meta: dict[str, str], reference_numbers: dict[str, int]
) -> str:
    central = optional_float(param["distance_mpc"])
    if central is None:
        return r"\nodata"
    stat_lower = optional_float(param["distance_stat_unc_minus_mpc"])
    stat_upper = optional_float(param["distance_stat_unc_plus_mpc"])
    if stat_lower is not None or stat_upper is not None:
        if stat_lower is None or stat_upper is None:
            raise ValueError(f"{param['source_key']}: incomplete distance interval")
    sys_lower = optional_float(param["distance_systematic_unc_minus_mpc"])
    sys_upper = optional_float(param["distance_systematic_unc_plus_mpc"])
    if sys_lower is not None or sys_upper is not None:
        if sys_lower is None or sys_upper is None:
            raise ValueError(f"{param['source_key']}: incomplete distance systematic")

    lower_components = [value for value in (stat_lower, sys_lower) if value is not None]
    upper_components = [value for value in (stat_upper, sys_upper) if value is not None]
    if lower_components:
        lower, upper = quadrature_interval_bounds(
            lower_components, upper_components
        )
        value = format_at_decimal_place(
            central, central_decimal_places([lower, upper])
        )
        value += uncertainty_tex(lower, upper)
    else:
        value = format_significant_figures(central, 3)
    return refmark(meta["distance_reference_key"], reference_numbers) + f"${value}$"


def size_stat_label(status: str) -> str:
    labels = {
        "paired_chain_quantile_interval": "stat.",
        "reported_undecomposed_uncertainty": "reported",
        "fit_plus_absolute_position_mixed_component": "fit+pos.",
        "fixed_distance_formal_fit_interval": "formal",
        "angular_radius_only_component": "angular radius",
    }
    try:
        return labels[status]
    except KeyError as exc:
        raise ValueError(f"unsupported size statistical component status {status!r}") from exc


def display_size_stat_bounds(
    param: dict[str, str], lower: float | None, upper: float | None
) -> tuple[float | None, float | None]:
    """Symmetrize a supported primary size interval for table display."""

    if lower is None and upper is None:
        return None, None
    if lower is None or upper is None:
        raise ValueError(f"{param['source_key']}: incomplete size interval")
    return symmetric_display_interval_bounds(lower, upper)


def size_uncertainty_marker_tex(marker: str) -> str:
    """Render an accepted primary-uncertainty qualifier inside math mode."""

    if marker == "":
        return ""
    if marker == "dagger":
        return r"^{\dagger}"
    raise ValueError(f"unsupported size uncertainty marker {marker!r}")


def format_size_cell(param: dict[str, str], marker: str = "") -> str:
    """Render the adopted size with its primary uncertainty in one cell."""

    central = optional_float(param["size_uas"])
    if central is None:
        return r"\nodata"
    stat_lower, stat_upper = display_size_stat_bounds(
        param,
        optional_float(param["size_stat_unc_minus_uas"]),
        optional_float(param["size_stat_unc_plus_uas"]),
    )

    sys_lower = optional_float(param["size_systematic_unc_minus_uas"])
    sys_upper = optional_float(param["size_systematic_unc_plus_uas"])
    sys_status = param["size_systematic_status"]
    if sys_lower is not None or sys_upper is not None:
        if sys_lower is None or sys_upper is None:
            raise ValueError(f"{param['source_key']}: incomplete size systematic")
        sys_lower, sys_upper = symmetric_display_interval_bounds(
            sys_lower, sys_upper
        )

    uncertainties = [
        value
        for value in (stat_lower, stat_upper, sys_lower, sys_upper)
        if value is not None
    ]
    if uncertainties:
        value = format_at_decimal_place(
            central, central_decimal_places(uncertainties)
        )
    else:
        value = format_significant_figures(central, 3)
    if (
        stat_lower is not None
        and stat_upper is not None
        and param["source_key"] not in INCOMPLETE_SIZE_TOTAL_KEYS
    ):
        value += uncertainty_tex(stat_lower, stat_upper)
    elif marker:
        raise ValueError(
            f"{param['source_key']}: uncertainty marker requires a displayed "
            "primary uncertainty"
        )
    value += size_uncertainty_marker_tex(marker)
    return f"${value}$"


def shortstack(lines: list[str]) -> str:
    """Wrap a compact publication cell without introducing fixed-width columns."""

    if not lines:
        raise ValueError("shortstack requires at least one line")
    if len(lines) == 1:
        return lines[0]
    return r"\shortstack[l]{" + r"\\".join(lines) + "}"


INCOMPLETE_SIZE_TOTAL_KEYS = {
    "NGC1068",
    "NGC4388",
    "Circinus",
    "NGC2960",
    "IC2560",
    "J0437",
}


def format_additional_size_term_cell(param: dict[str, str]) -> str:
    """Render the separate additional size term, if one is adopted."""

    key = param["source_key"]
    central = optional_float(param["size_uas"])
    if central is None:
        if key not in {"Mrk1029", "NGC1320"}:
            raise ValueError(f"{key}: unexpected missing adopted shadow size")
        return "---"

    one_sided_model_checks = {"Circinus", "NGC4945", "NGC3079"}
    if key in one_sided_model_checks:
        lower = optional_float(param["size_model_check_lower_uas"])
        upper = optional_float(param["size_model_check_upper_uas"])
        if lower is None or upper is None or upper <= lower:
            raise ValueError(f"{key}: one-sided model-check bounds are required")
        if not math.isclose(central, lower, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"{key}: model-check lower bound must equal the nominal size")
        shift = format_significant_figures(upper - lower, 2)
        return rf"$+{shift}$\tablenotemark{{b}}"

    stat_lower, stat_upper = display_size_stat_bounds(
        param,
        optional_float(param["size_stat_unc_minus_uas"]),
        optional_float(param["size_stat_unc_plus_uas"]),
    )
    sys_lower = optional_float(param["size_systematic_unc_minus_uas"])
    sys_upper = optional_float(param["size_systematic_unc_plus_uas"])
    stat_status = param["size_stat_component_status"]
    sys_status = param["size_systematic_status"]

    if stat_status == "no_formal_uncertainty":
        if key not in {"NGC4945", "NGC3079"}:
            raise ValueError(f"{key}: unexpected central-only shadow size")
        if any(value is not None for value in (stat_lower, stat_upper, sys_lower, sys_upper)):
            raise ValueError(f"{key}: central-only shadow size has a numerical component")
        return "---"

    if stat_lower is None or stat_upper is None:
        raise ValueError(f"{key}: displayed central size requires a supported component")

    supported_stat_statuses = {
        "paired_chain_quantile_interval",
        "reported_undecomposed_uncertainty",
        "fit_plus_absolute_position_mixed_component",
        "fixed_distance_formal_fit_interval",
        "angular_radius_only_component",
    }
    if stat_status not in supported_stat_statuses:
        raise ValueError(f"unsupported size statistical component status {stat_status!r}")

    if sys_lower is not None or sys_upper is not None:
        if sys_lower is None or sys_upper is None:
            raise ValueError(f"{key}: incomplete size systematic")
        if sys_status == "model_sensitivity_envelope":
            if key != "NGC3393":
                raise ValueError(f"{key}: unexpected displayed model component")
        elif sys_status == "numeric_one_percent_non_keplerian_allowance":
            if stat_status != "paired_chain_quantile_interval":
                raise ValueError(f"{key}: 1% term requires paired-chain statistics")
        elif sys_status == "numeric_systematic" and key == "UGC6093":
            pass
        elif sys_status == "numeric_systematic" and key == "ESO558-G009":
            pass
        else:
            raise ValueError(f"{key}: unsupported numerical size systematic")
        sys_lower, sys_upper = symmetric_display_interval_bounds(sys_lower, sys_upper)
        marker = {
            "numeric_one_percent_non_keplerian_allowance": "a",
            "model_sensitivity_envelope": "b",
        }.get(sys_status)
        if key == "UGC6093":
            marker = "c"
        elif key == "ESO558-G009":
            marker = "d"
        if marker is None:
            raise ValueError(f"{key}: displayed term lacks a table-note marker")
        return rf"${uncertainty_tex(sys_lower, sys_upper)}$\tablenotemark{{{marker}}}"

    incomplete = key in INCOMPLETE_SIZE_TOTAL_KEYS
    if incomplete:
        expected_incomplete_statuses = {
            "NGC1068": "unquantified_model",
            "NGC4388": "unquantified_model",
            "Circinus": "unquantified_model",
            "NGC2960": "unquantified_astrometric",
            "IC2560": "unquantified_model",
            "J0437": "unquantified_model",
        }
        if sys_status != expected_incomplete_statuses[key]:
            raise ValueError(
                f"{key}: incomplete-total marker does not match systematic status"
            )
    elif sys_lower is None:
        permitted_non_numeric = {
            "reported_undecomposed_uncertainty"
        } if stat_status == "reported_undecomposed_uncertainty" else {"none", "NA"}
        if sys_status not in permitted_non_numeric:
            raise ValueError(f"{key}: unsupported non-numeric size systematic status")

    return "---"


def ngc4258_model_check_note(params: list[dict[str, str]]) -> str:
    row = next(param for param in params if param["source_key"] == "NGC4258")
    lower = optional_float(row["size_model_check_lower_uas"])
    upper = optional_float(row["size_model_check_upper_uas"])
    if lower is None or upper is None:
        raise ValueError("NGC4258: model-check bounds are required")
    lower_text = format_significant_figures(lower, 4)
    upper_text = format_significant_figures(upper, 4)
    accepted_count = row["size_model_check_accepted_configuration_count"]
    citation = row["size_model_check_reference_key"]
    return (
        r"\tablenotetext{e}{The individual-component reanalysis of a subset of the "
        r"public 1997--2000 VLBI data by \citet{" + citation + r"} found "
        + accepted_count
        + r" statistically "
        r"acceptable data-selection variants whose published marginal medians imply "
        f"${lower_text}\\text{{--}}{upper_text}~\\mu\\mathrm{{as}}$.  This non-probabilistic range is "
        r"reported separately and is not combined with the paired-chain statistical "
        r"interval or 1\% non-Keplerian allowance.}"
    )


def ngc2960_model_check_note(params: list[dict[str, str]]) -> str:
    row = next(param for param in params if param["source_key"] == "NGC2960")
    lower = optional_float(row["size_model_check_lower_uas"])
    upper = optional_float(row["size_model_check_upper_uas"])
    if lower is None or upper is None:
        raise ValueError("NGC2960: reviewed disk-model check is required")
    lower_text = format_significant_figures(lower, 5)
    upper_text = format_significant_figures(upper, 5)
    citation = row["size_model_check_reference_key"]
    return (
        r"\tablenotetext{d}{For NGC~2960, the displayed interval is the "
        r"fixed-distance 68\% marginalized formal term from \citet{" + citation
        + r"}.  That paper attributes its input positions and velocities to "
        r"unpublished work and does not state whether the earlier 10~mas "
        r"absolute-coordinate contribution discussed by \citet{Kuo2010} was "
        r"resolved.  The legacy 3--5\% estimate is therefore not transplanted, "
        r"but a complete astrometric systematic remains unavailable.  The "
        r"physical-mass interval likewise combines only the available fixed-distance "
        r"fit and adopted-distance components and is not a complete mass-error "
        r"total.  The "
        r"density-constrained linear and Herrnstein disk models give "
        r"$1.232$ and $1.235\times10^7~M_\odot$, respectively, versus the "
        r"fiducial $1.237\times10^7~M_\odot$.  Their physically allowed "
        r"one-sided envelope is $"
        + lower_text
        + r"\text{--}"
        + upper_text
        + r"~\mu\mathrm{as}$ (a maximum 0.404\% downward shift).  It is "
        r"non-probabilistic, reported separately, and neither adopted nor "
        r"combined in quadrature; the unconstrained disk variants are excluded "
        r"because their inferred densities would quench masing.}"
    )


def ngc3079_model_check_note(params: list[dict[str, str]]) -> str:
    row = next(param for param in params if param["source_key"] == "NGC3079")
    lower = optional_float(row["size_model_check_lower_uas"])
    upper = optional_float(row["size_model_check_upper_uas"])
    if lower is None or upper is None:
        raise ValueError("NGC3079: model-check bounds are required")
    lower_text = format_significant_figures(lower, 3)
    upper_text = format_significant_figures(upper, 3)
    citation = row["size_model_check_reference_key"]
    return (
        r"\tablenotetext{f}{For NGC~3079, Appendix~A of \citet{" + citation
        + r"} gives $\Delta M_{\rm turb}/M_{\rm enc}=0.34$, so that "
        r"$M_{\rm BH}=1.34M_{\rm enc}$ under the stated model "
        r"assumptions $r_{\rm in}=0.40$~pc, $v_{\rm turb}=60$~km~s$^{-1}$, "
        r"$\Sigma\propto r^{-1}$, a nearly flat rotation curve "
        r"($\Omega\propto r^{-1}$, giving $\eta=2$ in the mass-bias expression), "
        r"and constant turbulent velocity, with turbulent pressure dominating "
        r"thermal pressure.  "
        r"Scaling the central shadow size by 1.34 gives a one-sided check of $"
        + lower_text
        + r"\text{--}"
        + upper_text
        + r"~\mu\mathrm{as}$.  This model-dependent range is non-probabilistic, "
        r"is reported separately, and is neither adopted as an uncertainty nor "
        r"combined in quadrature.  It does not bound the broader order-unity "
        r"geometry and disk-mass caveat in the original enclosed-mass analysis.}"
    )


def ngc6323_model_check_note(params: list[dict[str, str]]) -> str:
    row = next(param for param in params if param["source_key"] == "NGC6323")
    lower = optional_float(row["size_model_check_lower_uas"])
    upper = optional_float(row["size_model_check_upper_uas"])
    if lower is None or upper is None:
        raise ValueError("NGC6323: reviewed disk-model check is required")
    lower_text = format_significant_figures(lower, 5)
    upper_text = format_significant_figures(upper, 5)
    citation = row["size_model_check_reference_key"]
    return (
        r"\tablenotetext{j}{For NGC~6323, the adopted paired-chain result and "
        r"separate 1\% non-Keplerian allowance remain unchanged.  In a predecessor "
        r"full disk fit, \citet{" + citation
        + r"} obtain $0.984$ and $0.990\times10^7~M_\odot$ from the "
        r"density-constrained linear and Herrnstein disk models, respectively, "
        r"versus the $0.992\times10^7~M_\odot$ fiducial at 106.7~Mpc.  "
        r"Their physically allowed one-sided envelope is $"
        + lower_text
        + r"\text{--}"
        + upper_text
        + r"~\mu\mathrm{as}$, entirely within the adopted 1\% allowance.  It is "
        r"non-probabilistic, reported separately, and neither substituted for the "
        r"paired posterior nor combined in quadrature; the unconstrained disk "
        r"variants are excluded because their inferred densities are outside the "
        r"water-maser-supporting range.}"
    )


def j0437_model_check_note(params: list[dict[str, str]]) -> str:
    row = next(param for param in params if param["source_key"] == "J0437")
    lower = optional_float(row["size_model_check_lower_uas"])
    upper = optional_float(row["size_model_check_upper_uas"])
    if lower is None or upper is None:
        raise ValueError("J0437: reviewed Mestel-disk check is required")
    lower_text = format_significant_figures(lower, 4)
    upper_text = format_significant_figures(upper, 4)
    citation = row["size_model_check_reference_key"]
    return (
        r"\tablenotetext{k}{For J0437+2456, \citet{" + citation
        + r"} fit a one-dimensional, edge-on, coplanar black-hole plus Mestel-disk "
        r"model to literature maser data.  Their supplement reports "
        r"$M_{\rm BH}=(1.85\pm0.23)\times10^6~M_\odot$ and "
        r"$M_{\rm disk}=(0.82\pm0.26)\times10^6~M_\odot$ at an imported "
        r"$D=68.27\pm0.45$~Mpc.  Holding that reported distance fixed, as "
        r"required for the shared-distance mass normalization, gives the separate "
        r"$"
        + lower_text
        + r"\text{--}"
        + upper_text
        + r"~\mu\mathrm{as}$ black-hole-only span.  The exact retained spot set "
        r"and confidence semantics are not reported, so this is a "
        r"non-probabilistic model-sensitivity check rather than a complete "
        r"model envelope; it is neither adopted nor combined with the primary "
        r"formal interval.}"
    )


def mrk1029_model_check_note(params: list[dict[str, str]]) -> str:
    row = next(param for param in params if param["source_key"] == "Mrk1029")
    lower = optional_float(row["size_model_check_lower_uas"])
    upper = optional_float(row["size_model_check_upper_uas"])
    ratio_size = optional_float(row["nonadopted_fit_size_uas"])
    ratio_unc = optional_float(row["nonadopted_fit_size_unc_minus_uas"])
    if lower is None or upper is None or ratio_size is None or ratio_unc is None:
        raise ValueError("Mrk1029: both reviewed standalone diagnostics are required")
    lower_text = format_significant_figures(lower, 4)
    upper_text = format_significant_figures(upper, 4)
    ratio_text = format_at_decimal_place(
        ratio_size, central_decimal_places([ratio_unc])
    )
    ratio_unc_text = format_significant_figures(ratio_unc, 2)
    model_citation = row["size_model_check_reference_key"]
    ratio_citation = row["nonadopted_fit_reference_key"]
    return (
        r"\tablenotetext{l}{Mrk~1029 has no adopted black-hole mass or shadow size: "
        r"the source spectrum and mapped morphology do not support a robust disk "
        r"assignment.  The displayed $79\pm2\arcdeg$ inclination is Gao's "
        r"thin-disk geometry estimate, fixed rather than inferred in the dynamical "
        r"fit, and is therefore marked as an estimate.  The restricted high-velocity "
        r"fit of \citet{" + ratio_citation
        + r"} gives a conditional formal shadow size of $"
        + ratio_text
        + r"\pm"
        + ratio_unc_text
        + r"~\mu\mathrm{as}$ after fixing $i=79\arcdeg$, $y_0=0$, and a selected "
        r"scatter angle.  Separately, \citet{" + model_citation
        + r"} report a one-dimensional edge-on black-hole plus Mestel-disk "
        r"diagnostic; holding their imported distance fixed gives a black-hole-only "
        r"span of $"
        + lower_text
        + r"\text{--}"
        + upper_text
        + r"~\mu\mathrm{as}$.  The exact spot selection, distance provenance, and "
        r"confidence construction for that diagnostic are not reported.  Neither "
        r"result is adopted, combined, or plotted as a primary value.}"
    )


def ngc1320_model_check_note(params: list[dict[str, str]]) -> str:
    row = next(param for param in params if param["source_key"] == "NGC1320")
    lower = optional_float(row["size_model_check_lower_uas"])
    upper = optional_float(row["size_model_check_upper_uas"])
    ratio_size = optional_float(row["nonadopted_fit_size_uas"])
    ratio_minus = optional_float(row["nonadopted_fit_size_unc_minus_uas"])
    ratio_plus = optional_float(row["nonadopted_fit_size_unc_plus_uas"])
    if any(
        value is None
        for value in (lower, upper, ratio_size, ratio_minus, ratio_plus)
    ):
        raise ValueError("NGC1320: both reviewed standalone diagnostics are required")
    lower_text = format_significant_figures(float(lower), 4)
    upper_text = format_significant_figures(float(upper), 4)
    ratio_text = format_at_decimal_place(
        float(ratio_size),
        central_decimal_places([float(ratio_minus), float(ratio_plus)]),
    )
    ratio_minus_text = format_significant_figures(float(ratio_minus), 2)
    ratio_plus_text = format_significant_figures(float(ratio_plus), 2)
    model_citation = row["size_model_check_reference_key"]
    ratio_citation = row["nonadopted_fit_reference_key"]
    return (
        r"\tablenotetext{m}{NGC~1320 has no adopted black-hole mass or shadow "
        r"size: Gao's VLBI map lacks systemic masers, and its off-axis red- and "
        r"blueshifted spots allow a strong warp or non-disk/outflow origin.  The "
        r"inclination is therefore unconstrained; $i=90\arcdeg$ was fixed only as "
        r"a model assumption.  The restricted same-data fit of \citet{" + ratio_citation
        + r"} fits $M/D$, $z_{\rm rec}$, $x_0$, and $y_0$ while fixing "
        r"$i=90\arcdeg$ and selecting a $17\arcdeg$ scatter angle "
        r"($\chi^2_\nu=1.127$).  It gives the conditional formal shadow size $"
        + ratio_text
        + r"_{-"
        + ratio_minus_text
        + r"}^{+"
        + ratio_plus_text
        + r"}~\mu\mathrm{as}$; the reported $\pm10\arcdeg$ position-angle "
        r"variation changes $M/D$ by 2.64\% and is not added to that interval.  "
        r"Separately, \citet{" + model_citation
        + r"} report a one-dimensional edge-on black-hole plus Mestel-disk "
        r"diagnostic; holding their imported distance fixed gives a black-hole-only "
        r"span of $"
        + lower_text
        + r"\text{--}"
        + upper_text
        + r"~\mu\mathrm{as}$.  Its exact spot selection and confidence construction "
        r"are not reported.  Appendix~A of \citet{Kuo2020} also gives a conditional "
        r"1\% turbulent-support correction under its stated analytic assumptions; "
        r"because no robust primary mass exists, that result remains an audit-only "
        r"unattached check.  None of these diagnostics is adopted, combined, or "
        r"plotted as a primary value.}"
    )


def ngc4945_model_check_note(params: list[dict[str, str]]) -> str:
    row = next(param for param in params if param["source_key"] == "NGC4945")
    lower = optional_float(row["size_model_check_lower_uas"])
    upper = optional_float(row["size_model_check_upper_uas"])
    if lower is None or upper is None:
        raise ValueError("NGC4945: model-check bounds are required")
    lower_text = format_significant_figures(lower, 3)
    upper_text = format_significant_figures(upper, 3)
    citation = row["size_model_check_reference_key"]
    return (
        r"\tablenotetext{g}{For NGC~4945, Appendix~A of \citet{" + citation
        + r"} estimates that turbulent support could make the enclosed mass "
        r"inferred without that term 24\% too low under the stated model "
        r"assumptions $r_{\rm in}=0.16$~pc, $v_{\rm turb}=60$~km~s$^{-1}$, "
        r"$\Sigma\propto r^{-1}$, Keplerian rotation "
        r"($\Omega\propto r^{-3/2}$, giving $\eta=5/2$ in the mass-bias expression), "
        r"and constant turbulent velocity, with turbulent pressure dominating "
        r"thermal pressure.  "
        r"Scaling the central shadow size by 1.24 gives a one-sided check of $"
        + lower_text
        + r"\text{--}"
        + upper_text
        + r"~\mu\mathrm{as}$.  This model-dependent range is non-probabilistic, "
        r"is reported separately, and is neither adopted as an uncertainty nor "
        r"combined in quadrature.}"
    )


def circinus_model_check_note(params: list[dict[str, str]]) -> str:
    row = next(param for param in params if param["source_key"] == "Circinus")
    lower = optional_float(row["size_model_check_lower_uas"])
    upper = optional_float(row["size_model_check_upper_uas"])
    diagnostic_mass = optional_float(row["mass_total_unc_plus_msun"])
    if lower is None or upper is None or diagnostic_mass is None:
        raise ValueError("Circinus: reviewed check and mass diagnostic are required")
    lower_text = format_significant_figures(lower, 3)
    upper_text = format_significant_figures(upper, 3)
    diagnostic_mass_text = format_significant_figures(diagnostic_mass / 1e6, 3)
    citation = row["size_model_check_reference_key"]
    return (
        r"\tablenotetext{h}{For Circinus, the shadow-size component propagates the "
        r"reported inner angular radius $5.3\pm0.1$~mas after the shared distance "
        r"cancels; it is an angular-radius-only component, not a formal fit "
        r"interval.  The mass cell uses the source-reported "
        r"$(1.7\pm0.3)\times10^6~M_\odot$ interval.  The machine table also "
        r"retains a diagnostic decomposition using $0.1/5.3$ for the angular "
        r"radius and $0.8/4.2$ for distance, whose quadrature mass term is $"
        + diagnostic_mass_text
        + r"\times10^6~M_\odot$ and rounds to the reported $0.3\times10^6~M_\odot$; "
        r"that reconstructed value is not used for publication display.  "
        r"Appendix~A of \citet{" + citation
        + r"} estimates a 3\% turbulent-support mass bias for $r_{\rm in}=0.11$~pc, "
        r"$v_{\rm turb}=30$~km~s$^{-1}$, $\Sigma\propto r^{-1}$, Keplerian "
        r"rotation ($\eta=5/2$), and constant turbulent velocity with turbulent "
        r"pressure dominating thermal pressure.  Scaling the central shadow size "
        r"by 1.03 gives $"
        + lower_text
        + r"\text{--}"
        + upper_text
        + r"~\mu\mathrm{as}$.  This one-sided check is non-probabilistic, reported "
        r"separately, and neither adopted nor combined in quadrature.}"
    )


def ngc4388_model_check_note(params: list[dict[str, str]]) -> str:
    row = next(param for param in params if param["source_key"] == "NGC4388")
    lower = optional_float(row["size_model_check_lower_uas"])
    upper = optional_float(row["size_model_check_upper_uas"])
    if lower is None or upper is None:
        raise ValueError("NGC4388: reviewed model check is required")
    lower_text = format_significant_figures(lower, 4)
    upper_text = format_significant_figures(upper, 4)
    citation = row["size_model_check_reference_key"]
    return (
        r"\tablenotetext{i}{For NGC~4388, the $(0.84\pm0.02)\times10^7~M_\odot$ value in Table~4 "
        r"of \citet{Kuo2010} is adopted; Section~3.3.5 and Table~6 of that paper "
        r"instead give $0.85\times10^7~M_\odot$.  Its fixed-distance interval combines "
        r"Keplerian fitting and absolute source-position contributions.  The "
        r"adopted 19~Mpc Tully--Fisher distance has an 11\% rms uncertainty tied "
        r"to the Cepheid calibration and is therefore labeled systematic.  Only "
        r"five maser spots and no systemic features were detected: an edge-on "
        r"disk was assumed rather than measured, and the disk geometry and "
        r"Keplerian interpretation remain unquantified caveats.  Appendix~A of "
        r"\citet{" + citation
        + r"} estimates a 0.1\% turbulent-support mass bias for "
        r"$r_{\rm in}=0.24$~pc, $v_{\rm turb}=10$~km~s$^{-1}$, "
        r"$\Sigma\propto r^{-1}$, Keplerian rotation ($\eta=5/2$), and constant "
        r"turbulent velocity with turbulent pressure dominating thermal pressure.  "
        r"Scaling the central shadow size by 1.001 gives $"
        + lower_text
        + r"\text{--}"
        + upper_text
        + r"~\mu\mathrm{as}$.  This one-sided check is non-probabilistic, reported "
        r"separately, and neither adopted nor combined in quadrature.}"
    )


def format_inclination_cell(
    meta: dict[str, str], reference_numbers: dict[str, int]
) -> str:
    ref = refmark(meta["inclination_reference_key"], reference_numbers)
    key = meta["source_key"]
    qualitative = {
        "NGC1194": "nearly edge-on",
        "NGC3393": "nearly edge-on",
        "NGC4388": "unconstrained",
        "Circinus": "warped disk + outflow",
        "NGC4945": "nearly edge-on",
        "NGC2273": "nearly edge-on",
        "NGC2960": "nearly edge-on",
        "IC2560": "nearly edge-on",
        "NGC3079": "thick, self-gravitating disk",
        "J0437": "unconstrained",
        "NGC1320": "unconstrained",
    }
    if key in qualitative:
        return ref + qualitative[key]

    cgcg_combined = format_significant_figures(math.hypot(0.6, 1.1), 2)
    ugc6093_symmetric = format_significant_figures(3.95, 2)
    displays = {
        "NGC4258": rf"{ref}$87.05\pm0.09$ (at $r=6.1\,\mathrm{{mas}}$)",
        "NGC5765b": rf"{ref}$72.4\pm0.5$ (at center)",
        "CGCG074-064": rf"{ref}$90.8\pm{cgcg_combined}$ (at center)",
        "NGC6264": rf"{ref}$91.3\pm2.3$ (at center)",
        "UGC3789": rf"{ref}$84.9\pm0.6$ (at center)",
        "UGC6093": rf"{ref}$93.8\pm{ugc6093_symmetric}$",
        "ESO558-G009": rf"{ref}$98\pm1$ (geometry estimate)",
        "NGC6323": rf"{ref}$91.5\pm0.3$ (at center)",
        "Mrk1029": rf"{ref}$79\pm2$ (geometry estimate)",
        "NGC1068": rf"{ref}$\sim90$",
    }
    try:
        return displays[key]
    except KeyError as exc:
        raise ValueError(f"{key}: missing publication inclination display") from exc


def render_comparison_row(
    row: dict[str, str], reference_numbers: dict[str, int]
) -> str:
    size_tex = row["size_tex"] + size_uncertainty_marker_tex(
        row["size_uncertainty_marker"]
    )
    cells = [
        row["source_tex"],
        refmark(row["size_reference_key"], reference_numbers)
        + f"${size_tex}$",
        "---",
    ]
    for column in ("mass", "distance"):
        cells.append(
            refmark(row[f"{column}_reference_key"], reference_numbers)
            + f"${row[f'{column}_tex']}$"
        )
    cells.append(
        refmark(row["inclination_reference_key"], reference_numbers)
        + f"${row['inclination_tex']}"
        + "$"
    )
    return " & ".join(cells) + r" \\"


def render_sample_row(
    param: dict[str, str], meta: dict[str, str], reference_numbers: dict[str, int]
) -> str:
    source = meta["source_tex"]
    if meta["source_key"] in INDIVIDUALLY_DISCUSSED_SOURCE_KEYS:
        source = rf"\textbf{{{source}}}"
    cells = [
        source,
        format_size_cell(param, meta["size_uncertainty_marker"]),
        format_additional_size_term_cell(param),
        format_mass_cell(param, meta, reference_numbers),
        format_distance_cell(param, meta, reference_numbers),
        format_inclination_cell(meta, reference_numbers),
    ]
    return " & ".join(cells) + r" \\"


def render_table(
    params: list[dict[str, str]], metadata: list[dict[str, str]]
) -> str:
    validate_inputs(params, metadata)
    params_by_source = {row["source_key"]: row for row in params}
    metadata_by_source = {row["source_key"]: row for row in metadata}
    if set(params_by_source) != set(DISPLAY_SOURCE_KEYS):
        raise ValueError("parameter inventory differs from the publication display order")
    if set(metadata_by_source) != set(DISPLAY_SOURCE_KEYS):
        raise ValueError("metadata inventory differs from the publication display order")
    display_params = [params_by_source[key] for key in DISPLAY_SOURCE_KEYS]
    display_metadata = [metadata_by_source[key] for key in DISPLAY_SOURCE_KEYS]
    ordered_keys = reference_order(display_metadata)
    reference_numbers = {key: index for index, key in enumerate(ordered_keys, start=1)}

    lines = [
        r"\begin{longrotatetable}",
        r"\setlength{\tabcolsep}{2pt}",
        r"\begin{deluxetable*}{llllll}",
        r"\tablecaption{Predicted Schwarzschild Shadow Diameter, Black Hole Mass, Distance, and Maser Disk Geometry}",
        r"\label{tab:md}",
        r"\tabletypesize{\footnotesize}",
        r"\tablehead{",
        r"\colhead{Source} & \colhead{$\theta_{\rm BHS}$ [$\mu\mathrm{as}$]} & \colhead{\shortstack{Systematic/model term\\{}[$\mu\mathrm{as}$]}} & \colhead{$\log_{10}(M_\bullet/M_\odot)$} & \colhead{$D$ [Mpc]} & \colhead{$i$ [deg] / geometry}",
        r"}",
        r"\startdata",
    ]
    lines.extend(render_comparison_row(row, reference_numbers) for row in COMPARISON_ROWS)
    lines.extend(
        render_sample_row(param, meta, reference_numbers)
        for param, meta in zip(display_params, display_metadata, strict=True)
    )
    lines.extend(
        [
            r"\enddata",
            r"\end{deluxetable*}",
            r"\end{longrotatetable}",
            r"\tablecomments{Boldface identifies the six sources discussed individually in Section~\ref{sec:individual-sources}.  Unmarked uncertainties accompanying $\theta_{\rm BHS}$ are statistical.  A dagger marks an uncertainty that includes other effects or was published without separating them.  Lettered entries in the systematic/model column denote (a) a separate allowance for the fitted $M_\bullet/D$, and therefore the BHS size, to be 1\% higher or lower because of small imperfections in the Keplerian disk model, (b) a model sensitivity, (c) a position systematic, and (d) a geometry systematic; these terms remain separate from the uncertainty accompanying $\theta_{\rm BHS}$.  A leading $+$ denotes a one-sided increase in the predicted BHS size.  Asymmetric intervals are displayed using the average of their upper and lower half-widths; exact bounds and incomplete components remain machine-readable.}",
            r"\tablerefs{\footnotesize{" + "; ".join(
                f"({reference_numbers[key]})~\\citet{{{key}}}" for key in ordered_keys
            ) + ".}}",
            r"\tabletypesize{\scriptsize}",
        ]
    )
    return "\n".join(lines) + "\n"


def build_table(
    params_path: Path = DEFAULT_PARAMS,
    metadata_path: Path = DEFAULT_METADATA,
) -> str:
    return render_table(load_params(params_path), load_metadata(metadata_path))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--params", type=Path, default=DEFAULT_PARAMS)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="explicit manuscript Table 1 output path",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rendered = build_table(args.params, args.metadata)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
