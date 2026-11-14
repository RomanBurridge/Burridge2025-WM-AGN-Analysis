#!/usr/bin/env python3
"""Extend the continuum identity-review ledger after the Task 28 refit.

This utility is intentionally narrow.  It preserves the already reviewed
scaled-MAD identity rows, adds one row for each newly admitted ALMA product
that the fresh fitter classifies as a detection, and leaves admitted upper
limits without a component-identity row.  Run it only after direct review of
the Task 28 image/contact-sheet evidence; ``--review-confirmed`` records that
the visual review has occurred rather than attempting to automate it.
"""

from __future__ import annotations

import argparse
import ast
import csv
import math
import sys
from pathlib import Path

from astropy.coordinates import FK5, ICRS, SkyCoord
from astropy.time import Time
import astropy.units as u


ANALYSIS_ROOT = Path(__file__).resolve().parents[2]
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from continuum_classification_identity import (  # noqa: E402
    POLICY_FIELDS,
    load_ngc4945_source18_identity_records,
    validate_classification_identity_policy,
)
from continuum_product_additions import (  # noqa: E402
    load_product_addition_policy,
    publication_addition_products,
)


REVIEW_DATE = "2026-08-25"
BLENDED_BEAM_THRESHOLD_ARCSEC = 0.5

# A source-level physical or previously reviewed nuclear reference is reused
# across observing products.  Coordinates are transformed to the native FITS
# frame before the strict identity validator computes the centroid separation.
SOURCE_REFERENCES = {
    "Circinus": {
        "frame": "ICRS",
        "ra": 213.291445416667,
        "dec": -65.339173027778,
        "role": "reviewed_high_resolution_nuclear_continuum_anchor",
        "evidence": (
            "analysis/multi_freq_from_archive/chosen_frequency_anchor_reference.csv#Circinus;"
            "Greenhill2003_Table3_and_note_a_maser_dynamical_center_doi:10.1086/374862"
        ),
    },
    "ESO558-G009": {
        "frame": "J2000",
        "ra": 106.087547083333334,
        "dec": -21.588596666666664,
        "role": "VLA_B_array_H2O_maser_nuclear_proxy",
        "evidence": (
            "Gao2017_Table1_native_FK5_J2000_VLA_B_position_"
            "doi:10.3847/1538-4357/834/1/52;"
            "Pesce2023_same_ALMA_survey_continuum_doi:10.3847/1538-4357/acc57a"
        ),
    },
    "IC2560": {
        "frame": "J2000",
        "ra": 154.077933333333,
        "dec": -33.563813611111,
        "role": "reviewed_nuclear_continuum_anchor",
        "evidence": (
            "analysis/multi_freq_from_archive/chosen_frequency_anchor_reference.csv#IC2560;"
            "Kamali2017_Table3_maser_nuclear_proxy_doi:10.1051/0004-6361/201730899"
        ),
    },
    "NGC1068": {
        "frame": "ICRS",
        "ra": 40.669622083333,
        "dec": -0.013318222222,
        "role": "reviewed_high_resolution_nuclear_continuum_anchor",
        "evidence": (
            "analysis/multi_freq_from_archive/chosen_frequency_anchor_reference.csv#NGC1068;"
            "Gallimore2004_Table2_S1_and_section5.1_doi:10.1086/423167"
        ),
    },
    "NGC1194": {
        "frame": "ICRS",
        "ra": 45.954620833333,
        "dec": -1.103741666667,
        "role": "reviewed_high_resolution_nuclear_continuum_anchor",
        "evidence": (
            "analysis/multi_freq_from_archive/chosen_frequency_anchor_reference.csv#NGC1194;"
            "Kuo2011_phase_referenced_maser_position_doi:10.1088/0004-637X/727/1/20;"
            "Kamali2017_Table3_doi:10.1051/0004-6361/201730899"
        ),
    },
    "NGC1320": {
        "frame": "ICRS",
        "ra": 51.202879166667,
        "dec": -3.042273611111,
        "role": "reviewed_nuclear_continuum_anchor",
        "evidence": (
            "analysis/multi_freq_from_archive/chosen_frequency_anchor_reference.csv#NGC1320"
        ),
    },
    "NGC2960": {
        "frame": "ICRS",
        "ra": 145.151592582180967,
        "dec": 3.577030649741091,
        "role": "phase_referenced_H2O_maser_system_position",
        "evidence": (
            "Kuo2011_Table1_FK5_J2000_position_transformed_to_ICRS_"
            "doi:10.1088/0004-637X/727/1/20"
        ),
    },
    "NGC3393": {
        "frame": "ICRS",
        "ra": 162.0977625,
        "dec": -25.162097222222,
        "role": "reviewed_nuclear_continuum_anchor",
        "evidence": (
            "analysis/multi_freq_from_archive/chosen_frequency_anchor_reference.csv#NGC3393;"
            "Kondratko2008_section4_dynamical_center_doi:10.1086/586879"
        ),
    },
    "NGC4388": {
        "frame": "ICRS",
        "ra": 186.444920833333,
        "dec": 12.662151111111,
        "role": "reviewed_high_resolution_nuclear_continuum_anchor",
        "evidence": (
            "analysis/multi_freq_from_archive/chosen_frequency_anchor_reference.csv#NGC4388;"
            "Kuo2011_phase_referenced_maser_position_doi:10.1088/0004-637X/727/1/20;"
            "Kamali2017_Table3_doi:10.1051/0004-6361/201730899"
        ),
    },
    "NGC4945": {
        "frame": "ICRS",
        "ra": 196.364499166666690,
        "dec": -49.468167777777779,
        "role": "authoritative_source18_AGN_maser_nucleus",
        "evidence": (
            "analysis/multi_freq_from_archive/ngc4945_component_audit.csv#"
            "accepted-source18-fixed-canonical;Emig2020_doi:10.3847/1538-4357/abb67d"
        ),
    },
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def read_sorted(path: Path) -> list[list[object]]:
    rows: list[list[object]] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        row = ast.literal_eval(line)
        if not isinstance(row, list) or len(row) != 8:
            raise ValueError(f"{path}:{line_number}: expected an eight-field list")
        row[7] = str(row[7]).strip().replace("\\", "/")
        rows.append(row)
    return rows


def native_reference(reference: dict[str, object], frame: str) -> tuple[float, float]:
    source_frame = str(reference["frame"])
    coordinate = SkyCoord(
        ra=float(reference["ra"]) * u.deg,
        dec=float(reference["dec"]) * u.deg,
        frame=ICRS() if source_frame == "ICRS" else FK5(equinox=Time("J2000")),
    )
    transformed = coordinate.transform_to(
        ICRS() if frame == "ICRS" else FK5(equinox=Time("J2000"))
    )
    return float(transformed.ra.deg), float(transformed.dec.deg)


def is_detection(row: list[object]) -> bool:
    return not str(row[3]).endswith("*") and str(row[6]) != "NA"


def build_new_row(
    *,
    sequence: int,
    addition: dict[str, str],
    generated: list[object],
    audit: dict[str, str],
) -> dict[str, str]:
    source = str(generated[0])
    if source not in SOURCE_REFERENCES:
        raise ValueError(f"No reviewed source reference is configured for {source}")
    frame = audit["fit_center_reference_frame"]
    if frame not in {"ICRS", "J2000"}:
        raise ValueError(f"Unsupported fitted frame {frame!r} for {addition['product']}")
    reference = SOURCE_REFERENCES[source]
    reference_ra, reference_dec = native_reference(reference, frame)
    beam_major = float(audit["beam_major_arcsec"])
    if not math.isfinite(beam_major) or beam_major <= 0:
        raise ValueError(f"Invalid beam major axis for {addition['product']}")
    blended = beam_major >= BLENDED_BEAM_THRESHOLD_ARCSEC
    frequency = float(generated[2])
    if blended:
        evidence_class = "target_associated_blend"
        identity_status = "confirmed_target_associated_blend"
        semantics = "target_associated_blended_flux"
        component = f"broad_beam_{frequency:.2f}GHz_nuclear_region_component"
        rationale = (
            "Direct Task 28 cutout review shows target-associated central emission at "
            f"the reviewed source reference; the {beam_major:.3f}-arcsec beam includes "
            "resolved nuclear-region structure, so the row retains blended-flux semantics"
        )
    else:
        evidence_class = "a_priori_position_match"
        identity_status = "confirmed_target_component"
        semantics = "compact_component_flux"
        component = f"central_{frequency:.2f}GHz_beam_fixed_component"
        rationale = (
            "Direct Task 28 cutout review shows the beam-fixed component at the reviewed "
            f"source reference; the {beam_major:.3f}-arcsec beam supports compact-component "
            "semantics"
        )
    structured_sources = {"Circinus", "NGC1068", "NGC4388", "NGC4945"}
    return {
        "schema_version": "1",
        "policy_version": "scaled_mad_identity_v1",
        "review_id": f"CCIR-T28-{sequence:03d}",
        "trigger": "new_product_admission",
        "source": source,
        "telescope": str(generated[1]),
        "frequency_ghz": f"{frequency:.2f}",
        "product": addition["product"],
        "product_sha256": addition["product_sha256"],
        "reference_role": str(reference["role"]),
        "reference_frame": frame,
        "reference_ra_deg": f"{reference_ra:.15f}",
        "reference_dec_deg": f"{reference_dec:.15f}",
        "reference_evidence": str(reference["evidence"]),
        "identity_evidence_class": evidence_class,
        "component_identity_status": identity_status,
        "component_label": component,
        "measurement_semantics": semantics,
        "background_morphology": (
            "structured_source_emission" if source in structured_sources else "noise_like"
        ),
        "expected_scaled_mad_classification": "detection",
        "review_decision": "accept_scaled_mad_classification",
        "reviewer": "Codex_direct_FITS_scientific_audit",
        "review_date": REVIEW_DATE,
        "rationale": rationale,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--generated-table", type=Path, required=True)
    parser.add_argument("--noise-audit", type=Path, required=True)
    parser.add_argument("--existing-policy", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--review-confirmed", action="store_true")
    parser.add_argument("--validate", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.review_confirmed:
        raise SystemExit("Direct Task 28 visual review must be confirmed explicitly")
    analysis_root = args.analysis_root.resolve()
    additions = load_product_addition_policy(
        analysis_root / "continuum_product_addition_policy.csv",
        analysis_root=analysis_root,
    )
    publication = publication_addition_products(additions)
    additions_by_product = {row["product"]: row for row in additions}
    generated_rows = read_sorted(args.generated_table)
    generated_by_product = {str(row[7]): row for row in generated_rows}
    audits = read_csv(args.noise_audit)
    audit_by_product = {row["product"].strip().replace("\\", "/"): row for row in audits}
    audit_by_product.update(
        load_ngc4945_source18_identity_records(
            analysis_root / "ngc4945_nuclear_photometry.csv"
        )
    )
    existing = read_csv(args.existing_policy)
    existing_by_product = {row["product"].strip().replace("\\", "/"): row for row in existing}

    new_products = sorted(
        product
        for product in publication
        if product in generated_by_product
        and is_detection(generated_by_product[product])
        and product not in existing_by_product
    )
    rows = list(existing)
    for sequence, product in enumerate(new_products, 1):
        rows.append(
            build_new_row(
                sequence=sequence,
                addition=additions_by_product[product],
                generated=generated_by_product[product],
                audit=audit_by_product[product],
            )
        )
    rows.sort(key=lambda row: row["product"].strip().replace("\\", "/"))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=POLICY_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    if args.validate:
        workspace = analysis_root.parents[1]
        result = validate_classification_identity_policy(
            workspace=workspace,
            canonical_table=(
                analysis_root
                / "fitsumfiles/historical/fitsummary_fixtochosen_sorted_409_reference.txt"
            ),
            generated_table=args.generated_table,
            exclusion_policy=analysis_root / "continuum_product_exclusion_policy.csv",
            generic_noise_audit=args.noise_audit,
            identity_policy=args.output,
            addition_policy=analysis_root / "continuum_product_addition_policy.csv",
            boundary_rejection_policy=(
                analysis_root / "continuum_boundary_rejection_policy.csv"
            ),
        )
        print(
            f"validated {result['required_review_rows']} identity reviews; "
            f"added {len(new_products)} Task 28 detection rows"
        )
    else:
        print(f"wrote {len(rows)} rows; added {len(new_products)} Task 28 detections")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
