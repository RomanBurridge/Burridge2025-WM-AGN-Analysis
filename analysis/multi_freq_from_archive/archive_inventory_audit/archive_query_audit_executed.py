#!/usr/bin/env python3
"""Read-only ALMA archive snapshot comparison for the continuum inventory."""

from __future__ import annotations

import ast
import csv
import html
import io
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


AUDIT_DIRECTORY = Path(__file__).resolve().parent
BASE = AUDIT_DIRECTORY.parent
CUTOFF = "2025-01-27"
# Keep a live re-query separate from the bundled, frozen census.  This relative
# default also makes the preserved executable portable with the software tree.
OUTPUT = AUDIT_DIRECTORY / "alma_pre20250127_unmatched_mous_disposition_requery.csv"
TAP = "https://almascience.eso.org/tap/sync"
MIRRORS = ("https://almascience.eso.org", "https://almascience.nrao.edu")

CHEMICALS = (
    "SiO", "HeII", "CS", "H13CN", "HCN", "CO6", "13CO", "CH3CN",
    "CH3OH", "HC3N", "C18O", "CO32", "CS76", "HCOp43", "C17O", "CN",
    "H2CO", "NHC", "C3H2", "SO", "CO21", "CO", "HCOp", "H20", "HNC",
    "HC15N", "co32",
)
EARLY_EXCLUDES = ("_ph", "_bp", "README", "fluxcal", "phasecal")
CATEGORY_EXCLUDES = (
    "contsub", "tt0", "tt1", "alpha", "line", "mfs", "cube", "repBW",
    "sd", "concat", "split.cal", "Window_a", "SBring", "manual",
)


def source_coordinates() -> list[tuple[str, float, float]]:
    records = ast.literal_eval((BASE / "namesandcoords.txt").read_text())
    values = []
    for record in records:
        source = record[0][0]
        ra = [float(value) for value in record[1][1].split(":")]
        dec = record[1][2].split(":")
        ra_deg = 15.0 * (ra[0] + ra[1] / 60.0 + ra[2] / 3600.0)
        dec_deg = (
            (-1.0 if dec[0].startswith("-") else 1.0)
            * (abs(float(dec[0])) + float(dec[1]) / 60.0 + float(dec[2]) / 3600.0)
        )
        values.append((source, ra_deg, dec_deg))
    return values


def local_uids() -> dict[str, set[str]]:
    output: dict[str, set[str]] = {}
    table = BASE / "fitsumfiles/sorted/fitsummary_fixtochosen_sorted.txt"
    for line in table.read_text().splitlines():
        if not line.strip():
            continue
        row = ast.literal_eval(line)
        if row[1] != "ALMA":
            continue
        match = re.search(r"member\.uid___([^.]*)", row[7])
        if match:
            output.setdefault(row[0], set()).add(
                "uid://" + match.group(1).replace("_", "/")
            )
    output.setdefault("NGC4945", set()).add("uid://A002/X5a9a13/X528")
    return output


def tap_rows(source: str, ra_deg: float, dec_deg: float) -> list[dict[str, str]]:
    query = (
        "SELECT DISTINCT proposal_id,member_ous_uid,target_name,obs_release_date,"
        "schedblock_name,spatial_resolution,t_min,t_max,frequency,band_list,"
        "antenna_arrays FROM ivoa.obscore WHERE science_observation='T' "
        "AND data_rights='Public' AND INTERSECTS(CIRCLE('ICRS',"
        f"{ra_deg},{dec_deg},0.000555556),s_region)=1"
    )
    url = TAP + "?" + urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}
    )
    with urllib.request.urlopen(url, timeout=60) as handle:
        records = list(csv.DictReader(io.StringIO(handle.read().decode("utf-8"))))
    unique: dict[str, dict[str, str]] = {}
    for record in records:
        if record["obs_release_date"][:10] > CUTOFF:
            continue
        uid = record["member_ous_uid"]
        if uid not in unique:
            unique[uid] = {**record, "query_source": source}
        else:
            # Preserve the band/frequency ranges across spectral-window rows.
            for field in ("frequency", "band_list"):
                old = set(filter(None, unique[uid][field].split(";")))
                old.add(record[field])
                unique[uid][field] = ";".join(sorted(old))
    return list(unique.values())


def curl(url: str) -> tuple[bool, str]:
    command = [
        "curl", "--retry", "3", "--retry-all-errors", "--connect-timeout", "10",
        "--max-time", "60", "-fsSL", url,
    ]
    result = subprocess.run(command, text=True, capture_output=True)
    return result.returncode == 0, result.stdout


def datalink_text(uid: str) -> tuple[bool, str, str]:
    token = uid.replace("uid://", "uid___").replace("/", "_")
    errors = []
    for mirror in MIRRORS:
        ok, top = curl(f"{mirror}/datalink/sync?ID={token}")
        if not ok or not top.strip():
            errors.append(f"{mirror}:top_failed")
            continue
        ids = sorted(set(re.findall(r"/datalink/sync\?ID=([^\" <]+)", html.unescape(top))))
        if not ids:
            return True, top, mirror
        pieces = [top]
        all_ok = True
        for identifier in ids:
            ok, expanded = curl(f"{mirror}/datalink/sync?ID={identifier}")
            if not ok or not expanded.strip():
                all_ok = False
                errors.append(f"{mirror}:expand_failed:{identifier}")
                break
            pieces.append(expanded)
        if all_ok:
            return True, "\n".join(pieces), mirror
    return False, ";".join(errors), ""


def product_urls(text: str) -> list[str]:
    return sorted(
        set(
            html.unescape(value)
            for value in re.findall(r"https://[^<\" ]+", text)
            if "/dataPortal/" in value
        )
    )


def choose_candidate(source: str, candidates: list[str]) -> str:
    pool = list(candidates)
    numbers = re.findall(r"\d+", source)
    if numbers:
        width = max(len(value) for value in numbers)
        preferred = []
        for token in (value for value in numbers if len(value) == width):
            preferred.extend(value for value in pool if token.lower() in value.lower())
        if preferred:
            pool = preferred
    # A delivered QA2 image is the primary archival product when the same MOUS
    # also exposes an ARI-L reprocessing.  ARI-L remains an auditable alternate
    # reduction, not an additional publication row.  This reviewed rule matters
    # for uid://A001/X122/X600, whose original QA2 continuum is the admitted
    # product and chosen-frequency candidate.
    original_qa2 = [value for value in pool if "ari_l" not in value.lower()]
    if original_qa2:
        pool = original_qa2
    for token in ("pbcor", "regcal"):
        preferred = [value for value in pool if token in value.lower()]
        if preferred:
            pool = preferred
    return sorted(set(pool))[0]


def scientific_chemical_hit(url: str) -> bool:
    # The historical substring test treated the "SO" inside the galaxy-name
    # prefix "ESO" as a molecular-line token.  Remove that source-name prefix
    # before applying the otherwise unchanged chemical-name screen.
    without_eso_source_prefix = url.replace("ESO", "")
    return any(token in without_eso_source_prefix for token in CHEMICALS)


def classify(record: dict[str, str]) -> dict[str, str]:
    ok, payload, mirror = datalink_text(record["member_ous_uid"])
    base = {
        "source": record["query_source"],
        "member_ous_uid": record["member_ous_uid"],
        "proposal_id": record["proposal_id"],
        "target_name": record["target_name"],
        "obs_release_date": record["obs_release_date"],
        "public_before_cutoff": "yes",
        "schedblock_name": record["schedblock_name"],
        "spatial_resolution_arcsec": record["spatial_resolution"],
        "t_min_mjd": record["t_min"],
        "t_max_mjd": record["t_max"],
        "frequency_ghz_values": record["frequency"],
        "band_list": record["band_list"],
        "antenna_arrays": record["antenna_arrays"],
        "datalink_mirror": mirror,
    }
    if not ok:
        return {
            **base,
            "datalink_expansion_status": "failed",
            "eligible_pbcor_continuum": "unknown",
            "chosen_product_url": "",
            "alternative_eligible_urls": "",
            "historical_notebook_filename_disposition": "unknown",
            "exclusion_reason": "datalink_expansion_failed:" + payload,
        }
    urls = product_urls(payload)
    fits = [value for value in urls if ".fit" in value.lower()]
    initial = [
        value for value in fits if not any(token in value for token in EARLY_EXCLUDES)
    ]
    continuum_named = [value for value in initial if "cont" in value]
    after_categories = [
        value
        for value in continuum_named
        if not any(token in value for token in CATEGORY_EXCLUDES)
    ]
    candidates = [
        value
        for value in after_categories
        if "pbcor" in value.lower() and not scientific_chemical_hit(value)
    ]
    if candidates:
        chosen = choose_candidate(record["query_source"], candidates)
        alternatives = [value for value in sorted(set(candidates)) if value != chosen]
        literal_candidates = [
            value
            for value in candidates
            if not any(token in value for token in CHEMICALS)
        ]
        historical_disposition = (
            "would_pass_literal_notebook_filter"
            if literal_candidates
            else "historically_false_excluded_by_chemical_token_in_source_name"
        )
        return {
            **base,
            "datalink_expansion_status": "complete",
            "eligible_pbcor_continuum": "yes",
            "chosen_product_url": chosen,
            "alternative_eligible_urls": ";".join(alternatives),
            "historical_notebook_filename_disposition": historical_disposition,
            "exclusion_reason": "",
        }
    if not fits:
        reason = "no_fits_products"
    elif not continuum_named:
        reason = "no_continuum_named_fits"
    elif any("tt0" in value or "tt1" in value for value in continuum_named):
        reason = "only_taylor_term_or_nonintensity_continuum_fits"
    elif after_categories:
        reason = "no_pbcor_continuum_intensity_fits"
    else:
        reason = "continuum_fits_excluded_by_existing_filename_rules"
    return {
        **base,
        "datalink_expansion_status": "complete",
        "eligible_pbcor_continuum": "no",
        "chosen_product_url": "",
        "alternative_eligible_urls": "",
        "historical_notebook_filename_disposition": "not_eligible_after_other_rules",
        "exclusion_reason": reason,
    }


def main() -> int:
    included = local_uids()
    aliases = {
        "CircinusGalaxy": "Circinus",
        "WISEAJ043703.68+245607.0": "WISEAJ043703",
    }
    archive = []
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = [
            executor.submit(tap_rows, source, ra_deg, dec_deg)
            for source, ra_deg, dec_deg in source_coordinates()
        ]
        for future in as_completed(futures):
            archive.extend(future.result())
    unmatched = []
    seen = set()
    for record in archive:
        identity = (record["query_source"], record["member_ous_uid"])
        if identity in seen:
            continue
        seen.add(identity)
        key = aliases.get(record["query_source"], record["query_source"])
        if record["member_ous_uid"] not in included.get(key, set()):
            unmatched.append(record)
    completed: set[tuple[str, str]] = set()
    if OUTPUT.exists():
        with OUTPUT.open(newline="", encoding="utf-8") as handle:
            completed = {
                (row["source"], row["member_ous_uid"]) for row in csv.DictReader(handle)
            }
    pending = [
        record
        for record in unmatched
        if (record["query_source"], record["member_ous_uid"]) not in completed
    ][:12]
    print(
        f"unmatched={len(unmatched)} already_screened={len(completed)} "
        f"screening_batch={len(pending)}",
        flush=True,
    )
    rows = []
    with ThreadPoolExecutor(max_workers=6) as executor:
        futures = {executor.submit(classify, record): record for record in pending}
        for index, future in enumerate(as_completed(futures), 1):
            row = future.result()
            rows.append(row)
            if row["eligible_pbcor_continuum"] != "no":
                print(
                    index,
                    row["source"],
                    row["member_ous_uid"],
                    row["eligible_pbcor_continuum"],
                    row["chosen_product_url"],
                    row["exclusion_reason"],
                    flush=True,
                )
    if not rows:
        return 0
    fields = list(rows[0])
    write_header = not OUTPUT.exists()
    with OUTPUT.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        if write_header:
            writer.writeheader()
        writer.writerows(rows)
    print(
        f"appended {OUTPUT}: rows={len(rows)} eligible="
        f"{sum(row['eligible_pbcor_continuum'] == 'yes' for row in rows)} "
        f"unknown={sum(row['eligible_pbcor_continuum'] == 'unknown' for row in rows)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
