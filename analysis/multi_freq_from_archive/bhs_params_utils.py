from __future__ import annotations

import csv
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


def load_bhs_params(path: Path | None = None) -> dict[str, dict[str, object]]:
    csv_path = Path(path) if path is not None else DEFAULT_BHS_PATH
    if not csv_path.exists():
        raise FileNotFoundError(f"BHS parameter table not found: {csv_path}")

    rows: dict[str, dict[str, object]] = {}
    with csv_path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for raw_row in reader:
            row = {
                "source_key": raw_row["source_key"],
                "source_label": raw_row["source_label"],
                "size_uas": _to_float(raw_row["size_uas"]),
                "size_unc_independent_uas": _to_float(raw_row["size_unc_independent_uas"]),
                "size_unc_correlated_uas": _to_float(raw_row["size_unc_correlated_uas"]),
                "size_unc_adopted_uas": _to_float(raw_row["size_unc_adopted_uas"]),
                "adopted_uncertainty_model": raw_row["adopted_uncertainty_model"],
                "log_mass": _to_float(raw_row["log_mass"]),
                "log_mass_unc": _to_float(raw_row["log_mass_unc"]),
                "distance_mpc": _to_float(raw_row["distance_mpc"]),
                "distance_unc_mpc": _to_float(raw_row["distance_unc_mpc"]),
                "mass_frac_percent": _to_float(raw_row["mass_frac_percent"]),
            }
            rows[_normalize_name(row["source_key"])] = row
            rows[_normalize_name(row["source_label"])] = row
    return rows


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
