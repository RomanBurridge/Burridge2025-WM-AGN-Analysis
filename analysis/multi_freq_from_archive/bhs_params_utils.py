from __future__ import annotations

import csv
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BHS_PATH = ROOT / "analysis" / "CountingBHS_and_Spectra" / "BHS_params.csv"


def _normalize_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", str(value)).upper()


def _to_float(value: str):
    if value in ("", "NA", None):
        return None
    return float(value)


_ALIASES = {
    "WISEAJ043703": "J0437",
    "J0347": "J0437",
    "J0437": "J0437",
    "CIRCINUSGALAXY": "Circinus",
    "CIRCINUS": "Circinus",
    "CGCG07406": "CGCG074-064",
    "CGCG074064": "CGCG074-064",
    "07406": "CGCG074-064",
    "074064": "CGCG074-064",
    "ESO558": "ESO558-G009",
    "ESO558G009": "ESO558-G009",
    "MRK1029": "Mrk1029",
    "NGC5765": "NGC5765b",
}


_FLOAT_FIELDS = {
    "size_uas",
    "size_unc_independent_uas",
    "size_unc_correlated_uas",
    "size_unc_adopted_uas",
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
    "distance_unc_minus_mpc",
    "distance_unc_plus_mpc",
    "size_posterior_q16_uas",
    "size_posterior_q50_uas",
    "size_posterior_q84_uas",
    "size_unc_minus_uas",
    "size_unc_plus_uas",
    "mass_reference_distance_mpc",
    "mass_fit_unc_minus_msun",
    "mass_fit_unc_plus_msun",
    "mass_shared_distance_unc_minus_msun",
    "mass_shared_distance_unc_plus_msun",
    "mass_total_unc_minus_msun",
    "mass_total_unc_plus_msun",
    "mass_display_unc_minus_msun",
    "mass_display_unc_plus_msun",
    "distance_stat_unc_minus_mpc",
    "distance_stat_unc_plus_mpc",
    "distance_systematic_unc_minus_mpc",
    "distance_systematic_unc_plus_mpc",
    "size_stat_unc_minus_uas",
    "size_stat_unc_plus_uas",
    "size_systematic_unc_minus_uas",
    "size_systematic_unc_plus_uas",
    "size_total_unc_minus_uas",
    "size_total_unc_plus_uas",
    "size_model_check_accepted_configuration_count",
    "size_model_check_lower_uas",
    "size_model_check_upper_uas",
    "nonadopted_fit_ratio_msun_per_mpc",
    "nonadopted_fit_ratio_unc_minus_msun_per_mpc",
    "nonadopted_fit_ratio_unc_plus_msun_per_mpc",
    "nonadopted_fit_size_uas",
    "nonadopted_fit_size_unc_minus_uas",
    "nonadopted_fit_size_unc_plus_uas",
    "source_order",
}

_SIZE_STAT_COMPONENT_STATUSES = {
    "paired_chain_quantile_interval",
    "reported_undecomposed_uncertainty",
    "fit_plus_absolute_position_mixed_component",
    "fixed_distance_formal_fit_interval",
    "angular_radius_only_component",
    "no_formal_uncertainty",
    "not_available",
}


def _validate_size_stat_component(row: dict[str, object]) -> None:
    status = row.get("size_stat_component_status")
    if status not in _SIZE_STAT_COMPONENT_STATUSES:
        raise ValueError(
            f"unsupported size statistical component status {status!r} "
            f"for {row.get('source_label')}"
        )
    if row.get("source_key") == "NGC4388":
        if status != "fit_plus_absolute_position_mixed_component":
            raise ValueError(
                "NGC4388 reviewed size component must remain fit plus position"
            )
        if (
            row.get("distance_stat_component_status") != "not_available"
            or row.get("distance_stat_unc_minus_mpc") is not None
            or row.get("distance_stat_unc_plus_mpc") is not None
            or row.get("distance_systematic_unc_minus_mpc") != 2.09
            or row.get("distance_systematic_unc_plus_mpc") != 2.09
        ):
            raise ValueError(
                "NGC4388 reviewed Tully--Fisher distance term must remain a "
                "2.09 Mpc systematic"
            )


def load_bhs_params(path: Path | None = None) -> dict[str, dict[str, object]]:
    csv_path = Path(path) if path is not None else DEFAULT_BHS_PATH
    if not csv_path.exists():
        raise FileNotFoundError(f"BHS parameter table not found: {csv_path}")

    rows: dict[str, dict[str, object]] = {}
    with csv_path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for raw_row in reader:
            row = {
                key: _to_float(value) if key in _FLOAT_FIELDS else value
                for key, value in raw_row.items()
            }
            _validate_size_stat_component(row)
            rows[_normalize_name(row["source_key"])] = row
            rows[_normalize_name(row["source_label"])] = row
    return rows


def require_final_size_uncertainties(
    table: dict[str, dict[str, object]] | None = None,
) -> None:
    """Block a final release while an explicitly requested author term is pending.

    Rows whose source papers do not quantify a model systematic remain usable
    as formal-only results with a visible caveat.  This gate is narrower: it
    catches any row still carrying an explicit pending-author status.
    """

    if table is None:
        table = load_bhs_params()
    pending = sorted(
        {
            str(row["source_label"])
            for row in table.values()
            if row.get("size_systematic_status") == "pending_nakai"
        }
    )
    if pending:
        raise RuntimeError(
            "Final size-dependent output awaits an explicit formal-only or "
            f"systematic decision for: {', '.join(pending)}"
        )


def get_bhs_row(source_name: str, table: dict[str, dict[str, object]] | None = None):
    if table is None:
        table = load_bhs_params()

    normalized = _normalize_name(source_name)
    normalized = _normalize_name(_ALIASES.get(normalized, source_name))

    if normalized in table:
        return table[normalized]

    for key, row in table.items():
        if normalized in key or key in normalized:
            return row

    return None


def get_bhs_size_and_unc(source_name: str, table: dict[str, dict[str, object]] | None = None):
    row = get_bhs_row(source_name, table=table)
    if row is None:
        return None, None, None
    return row["size_uas"], row["size_unc_adopted_uas"], row


def get_bhs_values(source_name: str, table: dict[str, dict[str, object]] | None = None):
    row = get_bhs_row(source_name, table=table)
    if row is None or row["size_uas"] is None:
        return None, "NA", "NA", row

    size_uas = row["size_uas"]
    size_unc = row["size_unc_adopted_uas"]
    if size_unc is None:
        return size_uas, "NA", "NA", row

    frac_unc = (size_unc / size_uas) * 100 if size_uas != 0 else "NA"
    return size_uas, size_unc, frac_unc, row


def get_bhs_plot_uncertainties(
    source_name: str,
    table: dict[str, dict[str, object]] | None = None,
):
    """Return exact asymmetric linear and log-space plotting uncertainties.

    A numerical total is preferred.  If the total error budget is incomplete,
    the available uncertainty component is returned with ``complete`` set to
    ``False`` so callers can use a distinct line style and legend.
    """

    row = get_bhs_row(source_name, table=table)
    if row is None or row.get("size_uas") is None:
        return None
    _validate_size_stat_component(row)

    size = float(row["size_uas"])
    total_status = row.get("size_total_status")
    total_minus = row.get("size_total_unc_minus_uas")
    total_plus = row.get("size_total_unc_plus_uas")
    if total_status in {
        "quadrature_statistical_plus_systematic",
        "quadrature_mixed_fit_position_plus_systematic",
    }:
        if total_minus is None or total_plus is None:
            raise ValueError(f"calibrated total is missing for {row['source_label']}")
        minus = float(total_minus)
        plus = float(total_plus)
        complete = True
        component = "total"
    elif total_status == "quadrature_for_display_model_envelope_not_probability_interval":
        if total_minus is None or total_plus is None:
            raise ValueError(f"model display span is missing for {row['source_label']}")
        minus = float(total_minus)
        plus = float(total_plus)
        complete = False
        component = "model_span"
    elif total_status in {
        "formal_component_only",
        "reported_undecomposed_uncertainty_component_only",
    }:
        stat_minus = row.get("size_stat_unc_minus_uas")
        stat_plus = row.get("size_stat_unc_plus_uas")
        if stat_minus is None or stat_plus is None:
            raise ValueError(
                f"component-only uncertainty is missing for {row['source_label']}"
            )
        minus = float(stat_minus)
        plus = float(stat_plus)
        complete = False
        component = "component_only"
    elif total_status == "incomplete_unquantified_systematic":
        stat_minus = row.get("size_stat_unc_minus_uas")
        stat_plus = row.get("size_stat_unc_plus_uas")
        if stat_minus is None and stat_plus is None:
            if (
                row.get("source_key") != "NGC3079"
                or row.get("size_stat_component_status")
                != "no_formal_uncertainty"
                or row.get("size_systematic_status") != "unquantified_model"
            ):
                raise ValueError(
                    f"incomplete uncertainty has no retained component for "
                    f"{row['source_label']}"
                )
            return {
                "size_uas": size,
                "linear_minus_uas": None,
                "linear_plus_uas": None,
                "log10_minus": None,
                "log10_plus": None,
                "complete": False,
                "component": "none",
                "row": row,
            }
        if stat_minus is None or stat_plus is None:
            raise ValueError(
                f"component-only uncertainty is missing for {row['source_label']}"
            )
        minus = float(stat_minus)
        plus = float(stat_plus)
        complete = False
        component = "component_only"
    elif total_status == "no_formal_uncertainty":
        return {
            "size_uas": size,
            "linear_minus_uas": None,
            "linear_plus_uas": None,
            "log10_minus": None,
            "log10_plus": None,
            "complete": False,
            "component": "none",
            "row": row,
        }
    else:
        raise ValueError(
            f"unsupported size_total_status {total_status!r} for {row['source_label']}"
        )

    if minus < 0 or plus < 0 or minus >= size:
        raise ValueError(f"invalid size uncertainty for {row['source_label']}")

    log_size = math.log10(size)
    log_minus = log_size - math.log10(size - minus)
    log_plus = math.log10(size + plus) - log_size
    return {
        "size_uas": size,
        "linear_minus_uas": minus,
        "linear_plus_uas": plus,
        "log10_minus": log_minus,
        "log10_plus": log_plus,
        "complete": complete,
        "component": component,
        "row": row,
    }
