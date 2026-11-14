#!/usr/bin/env python3
"""Build the 19-row S230 component-identity and astrometry audit.

The builder joins the current S230 selection, the identity disposition audit,
the preserved-image WCS snapshot, and the source-level astrometry policy.  It
fails closed when a selected observation changes, a detection lacks component
evidence, or an upper limit accidentally acquires a fitted component.

No third-party package is required.  The FK5 J2000 to ICRS rotation implements
the USNO Circular 179 frame-bias angles used by Astropy for an FK5 coordinate
at equinox J2000.  This keeps the release audit reproducible outside CASA.
"""

from __future__ import annotations

import argparse
import csv
import io
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ANALYSIS_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = ANALYSIS_DIR.parents[1]
S230_TABLE = ANALYSIS_DIR / "machinetables/S230_extrapolation_comparison.csv"
IDENTITY_AUDIT = ANALYSIS_DIR / "s230_anchor_identity_audit.csv"
WCS_EVIDENCE = ANALYSIS_DIR / "s230_wcs_evidence.csv"
POLICY = ANALYSIS_DIR / "s230_astrometry_policy.csv"
ANCHOR_REFERENCE = ANALYSIS_DIR / "chosen_frequency_anchor_reference.csv"
FINAL_NOISE_AUDIT = (
    ANALYSIS_DIR
    / "fitsumfiles/noise_audits/230316_2_final_publication_noise_audit.csv"
)
OUTPUT = ANALYSIS_DIR / "machinetables/s230_anchor_astrometry_audit.csv"

NA = "NA"
COMMON_FRAME = "ICRS"
EXPECTED_COUNT = 19
EXPECTED_AMBIGUOUS = {"IC2560", "NGC1320"}
EXPECTED_BLENDED = {"NGC3079", "NGC4258"}
EXPECTED_HISTORICAL_ANCHORS = {
    "Circinus",
    "IC2560",
    "NGC1068",
    "NGC1194",
    "NGC1320",
    "NGC3393",
    "NGC4388",
    "NGC4945",
    "NGC5765B",
}
TASK28_REPLACED_ANCHORS = {"NGC1068", "NGC1194"}
UPPER_LIMIT_DISPOSITIONS = {"upper_limit", "conservative_upper_limit"}
DETECTION_DISPOSITIONS = {"nuclear", "ambiguous", "blended"}

WCS_OUTPUT_FIELDS = (
    "recovery_record_id",
    "product_sha256",
    "product_sha256_algorithm",
    "wcs_evidence_status",
    "wcs_direction_reference_code",
    "fits_radesys_keyword",
    "fits_equinox_keyword",
    "fits_ctype1",
    "fits_ctype2",
    "wcs_projection",
    "wcs_axis1_name",
    "wcs_axis2_name",
    "wcs_axis1_unit",
    "wcs_axis2_unit",
    "wcs_crpix1_zero_based",
    "wcs_crpix2_zero_based",
    "wcs_crval1_deg",
    "wcs_crval2_deg",
    "wcs_cdelt1_deg_per_pixel",
    "wcs_cdelt2_deg_per_pixel",
    "wcs_pc11",
    "wcs_pc12",
    "wcs_pc21",
    "wcs_pc22",
    "wcs_transform_convention",
)

FIELDS = (
    "source",
    "selected_flux_status",
    "selected_observation_id",
    "selected_reference",
    "selected_telescope",
    "selected_date",
    "selected_frequency_ghz",
    "selected_geometric_beam_arcsec",
    "identity_disposition",
    "component_identity",
    "identity_evidence_basis",
    "required_handling",
    "disposition_evidence",
    "component_record_status",
    "component_coordinate_method",
    "component_coordinate_evidence",
    "component_native_frame",
    "component_ra_native_deg",
    "component_dec_native_deg",
    "component_common_frame",
    "component_ra_icrs_deg",
    "component_dec_icrs_deg",
    "component_position_uncertainty_status",
    "component_position_uncertainty_detail",
    "reference_coordinate_role",
    "reference_coordinate_method",
    "reference_coordinate_citation",
    "reference_coordinate_evidence",
    "reference_native_frame",
    "reference_ra_native_deg",
    "reference_dec_native_deg",
    "reference_common_frame",
    "reference_ra_icrs_deg",
    "reference_dec_icrs_deg",
    "reference_position_uncertainty_status",
    "reference_position_uncertainty_detail",
    "authoritative_nuclear_coordinate_status",
    "coordinate_transform_method",
    "signed_east_offset_arcsec",
    "signed_north_offset_arcsec",
    "total_separation_arcsec",
    "separation_geometric_beams",
    "offset_interpretation",
    "positional_significance_status",
    "selected_image_phase_calibrator_status",
    "selected_image_phase_calibrator",
    "selected_image_phase_calibrator_evidence",
    "selected_image_astrometric_qa_status",
    "selected_image_astrometric_qa_evidence",
    *WCS_OUTPUT_FIELDS,
    "accepted_component_record",
    "rejected_component_record",
    "rejected_component_native_frame",
    "rejected_component_ra_native_deg",
    "rejected_component_dec_native_deg",
    "rejected_component_common_frame",
    "rejected_component_ra_icrs_deg",
    "rejected_component_dec_icrs_deg",
    "rejected_component_reason",
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"missing or duplicate CSV fields in {path}")
        rows: list[dict[str, str]] = []
        for line_number, raw in enumerate(reader, 2):
            if None in raw:
                raise ValueError(f"unexpected extra CSV values at {path}:{line_number}")
            # Older source tables use empty cells for unavailable values; the
            # generated audit normalizes every such value to the explicit NA
            # sentinel required by its schema.
            row = {key: ((value or "").strip() or NA) for key, value in raw.items()}
            rows.append(row)
        return rows


def _index_unique(rows: Sequence[Mapping[str, str]], key: str, label: str) -> dict[str, Mapping[str, str]]:
    indexed: dict[str, Mapping[str, str]] = {}
    for row in rows:
        value = row[key]
        if value in indexed:
            raise ValueError(f"duplicate {label} key {value!r}")
        indexed[value] = row
    return indexed


def _matmul(left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]) -> list[list[float]]:
    return [
        [sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3)]
        for i in range(3)
    ]


def _transpose(matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    return [list(row) for row in zip(*matrix)]


def _rotation_x(angle: float) -> list[list[float]]:
    cosine, sine = math.cos(angle), math.sin(angle)
    return [[1.0, 0.0, 0.0], [0.0, cosine, sine], [0.0, -sine, cosine]]


def _rotation_y(angle: float) -> list[list[float]]:
    cosine, sine = math.cos(angle), math.sin(angle)
    return [[cosine, 0.0, -sine], [0.0, 1.0, 0.0], [sine, 0.0, cosine]]


def _rotation_z(angle: float) -> list[list[float]]:
    cosine, sine = math.cos(angle), math.sin(angle)
    return [[cosine, sine, 0.0], [-sine, cosine, 0.0], [0.0, 0.0, 1.0]]


def _fk5_j2000_to_icrs_matrix() -> list[list[float]]:
    milliarcsec_to_rad = math.pi / (180.0 * 3600.0 * 1000.0)
    eta0 = -19.9 * milliarcsec_to_rad
    xi0 = 9.1 * milliarcsec_to_rad
    da0 = -22.9 * milliarcsec_to_rad
    icrs_to_fk5 = _matmul(
        _matmul(_rotation_x(-eta0), _rotation_y(xi0)),
        _rotation_z(da0),
    )
    return _transpose(icrs_to_fk5)


FK5_J2000_TO_ICRS = _fk5_j2000_to_icrs_matrix()


def _coordinate_to_icrs(frame: str, ra_deg: str, dec_deg: str) -> tuple[float, float, str]:
    ra = float(ra_deg)
    dec = float(dec_deg)
    if not math.isfinite(ra) or not math.isfinite(dec) or not 0.0 <= ra < 360.0 or not -90.0 <= dec <= 90.0:
        raise ValueError(f"invalid {frame} coordinate ({ra_deg}, {dec_deg})")
    if frame == "ICRS":
        return ra, dec, "ICRS_identity"
    if frame != "FK5_J2000":
        raise ValueError(f"unsupported native coordinate frame {frame!r}")

    ra_rad = math.radians(ra)
    dec_rad = math.radians(dec)
    vector = (
        math.cos(dec_rad) * math.cos(ra_rad),
        math.cos(dec_rad) * math.sin(ra_rad),
        math.sin(dec_rad),
    )
    transformed = [
        sum(FK5_J2000_TO_ICRS[i][j] * vector[j] for j in range(3))
        for i in range(3)
    ]
    out_ra = math.degrees(math.atan2(transformed[1], transformed[0])) % 360.0
    out_dec = math.degrees(
        math.atan2(transformed[2], math.hypot(transformed[0], transformed[1]))
    )
    return (
        out_ra,
        out_dec,
        "FK5_J2000_to_ICRS_USNO_Circular_179_frame_bias_Astropy_equivalent",
    )


def _offsets(component: tuple[float, float], reference: tuple[float, float]) -> tuple[float, float, float]:
    component_ra, component_dec = map(math.radians, component)
    reference_ra, reference_dec = map(math.radians, reference)
    delta_ra = math.atan2(
        math.sin(component_ra - reference_ra),
        math.cos(component_ra - reference_ra),
    )
    arcsec_per_rad = 180.0 * 3600.0 / math.pi
    east = delta_ra * math.cos((component_dec + reference_dec) / 2.0) * arcsec_per_rad
    north = (component_dec - reference_dec) * arcsec_per_rad
    haversine = (
        math.sin((component_dec - reference_dec) / 2.0) ** 2
        + math.cos(reference_dec)
        * math.cos(component_dec)
        * math.sin(delta_ra / 2.0) ** 2
    )
    total = 2.0 * math.asin(min(1.0, math.sqrt(haversine))) * arcsec_per_rad
    return east, north, total


def _number(value: float, places: int) -> str:
    if abs(value) < 0.5 * 10 ** (-places):
        value = 0.0
    return f"{value:.{places}f}"


def _parse_sexagesimal(ra: str, dec: str) -> tuple[float, float]:
    ra_parts = [float(part) for part in ra.split(":")]
    dec_sign = -1.0 if dec.startswith("-") else 1.0
    dec_parts = [float(part) for part in dec.lstrip("+-").split(":")]
    if len(ra_parts) != 3 or len(dec_parts) != 3:
        raise ValueError(f"invalid sexagesimal coordinate ({ra!r}, {dec!r})")
    ra_deg = 15.0 * (ra_parts[0] + ra_parts[1] / 60.0 + ra_parts[2] / 3600.0)
    dec_deg = dec_sign * (dec_parts[0] + dec_parts[1] / 60.0 + dec_parts[2] / 3600.0)
    return ra_deg, dec_deg


def _resolve_record(reference: str) -> Mapping[str, str]:
    relative, separator, record_id = reference.partition("#")
    if not separator or not relative or not record_id:
        raise ValueError(f"invalid record reference {reference!r}")
    path = REPOSITORY_ROOT / relative
    rows = _read_csv(path)
    matches = [row for row in rows if row.get("record_id") == record_id]
    if len(matches) != 1:
        raise ValueError(f"record reference does not resolve uniquely: {reference}")
    return matches[0]


def _validate_source_sets(
    selected: Sequence[Mapping[str, str]],
    identity: Mapping[str, Mapping[str, str]],
    wcs: Mapping[str, Mapping[str, str]],
    policy: Mapping[str, Mapping[str, str]],
) -> None:
    if len(selected) != EXPECTED_COUNT:
        raise ValueError(f"expected {EXPECTED_COUNT} selected S230 rows; found {len(selected)}")
    selected_sources = [row["Source"] for row in selected]
    if len(selected_sources) != len(set(selected_sources)):
        raise ValueError("duplicate sources in current S230 selection")
    expected = set(selected_sources)
    for label, indexed in (("identity", identity), ("WCS", wcs), ("policy", policy)):
        if set(indexed) != expected:
            missing = sorted(expected - set(indexed))
            stale = sorted(set(indexed) - expected)
            raise ValueError(f"{label} source mismatch; missing={missing}; stale={stale}")


def _validate_ngc4945(policy: Mapping[str, str]) -> None:
    accepted = _resolve_record(policy["accepted_component_record"])
    rejected = _resolve_record(policy["rejected_component_record"])
    if accepted["record_id"] != "accepted-source18-fixed-canonical" or accepted["decision"] != "accepted":
        raise ValueError("NGC4945 accepted source18 record is not canonical")
    if rejected["record_id"] != "rejected-source17-archived-anchor" or rejected["decision"] != "rejected":
        raise ValueError("NGC4945 rejected source17 record is not preserved")
    accepted_coord = _parse_sexagesimal(accepted["ra"], accepted["dec"])
    rejected_coord = _parse_sexagesimal(rejected["ra"], rejected["dec"])
    policy_accepted = (float(policy["component_ra_native_deg"]), float(policy["component_dec_native_deg"]))
    policy_rejected = (
        float(policy["rejected_component_ra_native_deg"]),
        float(policy["rejected_component_dec_native_deg"]),
    )
    for label, actual, expected in (
        ("accepted source18", policy_accepted, accepted_coord),
        ("rejected source17", policy_rejected, rejected_coord),
    ):
        if any(abs(left - right) > 5e-13 for left, right in zip(actual, expected)):
            raise ValueError(f"NGC4945 {label} coordinate drifted: {actual} != {expected}")


def _validate_historical_anchor_coordinates(
    selected_by_source: Mapping[str, Mapping[str, str]],
    wcs_by_source: Mapping[str, Mapping[str, str]],
    policy_by_source: Mapping[str, Mapping[str, str]],
) -> None:
    """Lock policy coordinates to the current high-resolution anchor snapshot.

    NGC 4945 is intentionally the exception: its historical anchor is the
    rejected source17 seed, and the canonical policy coordinate comes from the
    accepted source18 sidecar validated above.
    """

    anchors = _index_unique(_read_csv(ANCHOR_REFERENCE), "source", "anchor reference")
    if set(anchors) != EXPECTED_HISTORICAL_ANCHORS:
        raise ValueError(
            "historical anchor source mismatch; "
            f"missing={sorted(EXPECTED_HISTORICAL_ANCHORS - set(anchors))}; "
            f"stale={sorted(set(anchors) - EXPECTED_HISTORICAL_ANCHORS)}"
        )
    final_noise_by_product = _index_unique(
        _read_csv(FINAL_NOISE_AUDIT),
        "product",
        "final-noise product",
    )
    for source, anchor in anchors.items():
        selected = selected_by_source[source]
        wcs = wcs_by_source[source]
        policy = policy_by_source[source]
        replacement = source in TASK28_REPLACED_ANCHORS
        if replacement:
            if anchor["input_product"] == selected["nu0_reference"]:
                raise ValueError(f"{source}: expected the reviewed Task 28 anchor replacement")
            noise = final_noise_by_product.get(selected["nu0_reference"])
            if noise is None or noise["source"] != source:
                raise ValueError(
                    f"{source}: replacement component is absent from the final-noise audit"
                )
            if (
                noise["classification"] != "detection"
                or noise["guard_status"] != "passed"
            ):
                raise ValueError(
                    f"{source}: replacement component lacks a validated final detection"
                )
            expected_evidence = (
                "analysis/multi_freq_from_archive/fitsumfiles/noise_audits/"
                "230316_2_final_publication_noise_audit.csv#"
                f"{selected['nu0_reference']}"
            )
            if (
                policy["component_coordinate_method"]
                != "fresh_final_fixed_beam_gaussian_centroid"
                or policy["component_coordinate_evidence"] != expected_evidence
            ):
                raise ValueError(
                    f"{source}: replacement component provenance is not bound to the final fit"
                )
            expected_frame = (
                "FK5_J2000"
                if noise["fit_center_reference_frame"] == "J2000"
                else noise["fit_center_reference_frame"]
            )
            if policy["component_native_frame"] != expected_frame:
                raise ValueError(
                    f"{source}: replacement component frame drifted from the final fit"
                )
            policy_coord = (
                float(policy["component_ra_native_deg"]),
                float(policy["component_dec_native_deg"]),
            )
            noise_coord = (
                float(noise["fit_center_ra_deg"]),
                float(noise["fit_center_dec_deg"]),
            )
            if any(
                abs(left - right) > 5e-13
                for left, right in zip(policy_coord, noise_coord)
            ):
                raise ValueError(
                    f"{source}: replacement component coordinate drifted from the final fit"
                )
            continue
        else:
            if anchor["input_product"] != selected["nu0_reference"]:
                raise ValueError(f"{source}: historical anchor points to a stale selected product")
            if anchor["input_sha256"] != wcs["product_sha256"]:
                raise ValueError(f"{source}: historical anchor input hash drifted")
        anchor_coord = (float(anchor["anchor_ra_deg"]), float(anchor["anchor_dec_deg"]))
        if source == "NGC4945":
            if anchor["scientific_status"] != "rejected_component_seed":
                raise ValueError("NGC4945 historical source17 anchor lost its rejected status")
            policy_coord = (
                float(policy["rejected_component_ra_native_deg"]),
                float(policy["rejected_component_dec_native_deg"]),
            )
        else:
            expected_frame = "FK5_J2000" if anchor["anchor_frame"] == "J2000" else anchor[
                "anchor_frame"
            ]
            if policy["component_native_frame"] != expected_frame:
                raise ValueError(f"{source}: component frame drifted from anchor snapshot")
            policy_coord = (
                float(policy["component_ra_native_deg"]),
                float(policy["component_dec_native_deg"]),
            )
        if any(abs(left - right) > 5e-13 for left, right in zip(policy_coord, anchor_coord)):
            raise ValueError(f"{source}: policy component coordinate drifted from anchor snapshot")


def build_rows() -> list[dict[str, str]]:
    selected_rows = _read_csv(S230_TABLE)
    selected_by_source = _index_unique(selected_rows, "Source", "selected S230")
    identity_by_source = _index_unique(_read_csv(IDENTITY_AUDIT), "source", "identity")
    wcs_by_source = _index_unique(_read_csv(WCS_EVIDENCE), "source", "WCS")
    policy_by_source = _index_unique(_read_csv(POLICY), "source", "policy")
    _validate_source_sets(selected_rows, identity_by_source, wcs_by_source, policy_by_source)

    ambiguous = {
        source for source, row in identity_by_source.items() if row["identity_disposition"] == "ambiguous"
    }
    blended = {
        source for source, row in identity_by_source.items() if row["identity_disposition"] == "blended"
    }
    if ambiguous != EXPECTED_AMBIGUOUS or blended != EXPECTED_BLENDED:
        raise ValueError(f"material disposition drift: ambiguous={ambiguous}; blended={blended}")
    _validate_ngc4945(policy_by_source["NGC4945"])
    _validate_historical_anchor_coordinates(
        selected_by_source,
        wcs_by_source,
        policy_by_source,
    )

    output: list[dict[str, str]] = []
    for selected in selected_rows:
        source = selected["Source"]
        identity = identity_by_source[source]
        wcs = wcs_by_source[source]
        policy = policy_by_source[source]
        for label, actual, expected in (
            ("flux status", identity["selected_flux_status"], selected["flux_status"]),
            ("identity observation ID", identity["selected_observation_id"], selected["nu0_observation_id"]),
            ("identity reference", identity["selected_reference"], selected["nu0_reference"]),
            ("WCS observation ID", wcs["selected_observation_id"], selected["nu0_observation_id"]),
            ("WCS reference", wcs["selected_reference"], selected["nu0_reference"]),
        ):
            if actual != expected:
                raise ValueError(f"{source}: stale {label}: {actual!r} != {expected!r}")
        literature_only = wcs["wcs_evidence_status"] == "literature_only_no_local_product"
        if literature_only:
            if source != "NGC3079":
                raise ValueError(f"{source}: unexpected literature-only S230 reference")
            if wcs["product_sha256"] != NA or wcs["product_sha256_algorithm"] != "not_applicable":
                raise ValueError(f"{source}: literature-only WCS hash fields are stale")
        else:
            if wcs["wcs_evidence_status"] != "local_product_hash_verified":
                raise ValueError(f"{source}: selected image WCS is not hash-verified")
            if len(wcs["product_sha256"]) != 64:
                raise ValueError(f"{source}: invalid selected-product SHA-256")
        if any("/home/" in value for value in (*identity.values(), *wcs.values(), *policy.values())):
            raise ValueError(f"{source}: workstation path leaked into audit inputs")

        disposition = identity["identity_disposition"]
        is_upper_limit = selected["flux_status"] == "upper_limit"
        if is_upper_limit:
            if disposition not in UPPER_LIMIT_DISPOSITIONS:
                raise ValueError(f"{source}: upper limit has detection disposition {disposition!r}")
            if policy["component_record_status"] != "upper_limit_no_component":
                raise ValueError(f"{source}: upper limit unexpectedly has a component record")
            coordinate_fields = (
                "component_native_frame",
                "component_ra_native_deg",
                "component_dec_native_deg",
                "reference_native_frame",
                "reference_ra_native_deg",
                "reference_dec_native_deg",
            )
            if any(policy[field] != NA for field in coordinate_fields):
                raise ValueError(f"{source}: upper limit contains an adopted coordinate")
            component_common = reference_common = (NA, NA)
            component_frame = reference_frame = NA
            transform_method = NA
            east = north = total = beams = NA
            significance = "NOT_APPLICABLE_NO_COMPONENT"
        elif literature_only:
            if disposition != "blended":
                raise ValueError(
                    f"{source}: literature-only aperture measurement must retain blended semantics"
                )
            if policy["component_record_status"] != "literature_detection_no_astrometric_component":
                raise ValueError(
                    f"{source}: literature-only detection has an unsupported component record"
                )
            coordinate_fields = (
                "component_native_frame",
                "component_ra_native_deg",
                "component_dec_native_deg",
                "reference_native_frame",
                "reference_ra_native_deg",
                "reference_dec_native_deg",
            )
            if any(policy[field] != NA for field in coordinate_fields):
                raise ValueError(
                    f"{source}: literature-only detection contains unreported image astrometry"
                )
            component_common = reference_common = (NA, NA)
            component_frame = reference_frame = NA
            transform_method = NA
            east = north = total = beams = NA
            significance = "NOT_COMPUTABLE_NO_LOCAL_IMAGE_ASTROMETRY"
        else:
            if disposition not in DETECTION_DISPOSITIONS:
                raise ValueError(f"{source}: detection has unsupported disposition {disposition!r}")
            if policy["component_record_status"] != "detection_component":
                raise ValueError(f"{source}: detection lacks a component record")
            required_evidence = (
                "component_coordinate_method",
                "component_coordinate_evidence",
                "disposition_evidence",
                "reference_coordinate_role",
                "reference_coordinate_citation",
            )
            if any(policy[field] == NA for field in required_evidence):
                raise ValueError(f"{source}: detection lacks evidence-backed disposition")
            if policy["component_native_frame"].replace("FK5_", "") != wcs[
                "wcs_direction_reference_code"
            ]:
                raise ValueError(
                    f"{source}: component frame {policy['component_native_frame']} "
                    f"does not match selected-image WCS {wcs['wcs_direction_reference_code']}"
                )
            component_ra, component_dec, component_transform = _coordinate_to_icrs(
                policy["component_native_frame"],
                policy["component_ra_native_deg"],
                policy["component_dec_native_deg"],
            )
            reference_ra, reference_dec, reference_transform = _coordinate_to_icrs(
                policy["reference_native_frame"],
                policy["reference_ra_native_deg"],
                policy["reference_dec_native_deg"],
            )
            component_common = (_number(component_ra, 12), _number(component_dec, 12))
            reference_common = (_number(reference_ra, 12), _number(reference_dec, 12))
            component_frame = reference_frame = COMMON_FRAME
            transform_method = f"component={component_transform};reference={reference_transform}"
            east_value, north_value, total_value = _offsets(
                (component_ra, component_dec), (reference_ra, reference_dec)
            )
            beam = float(selected["nu0_beam_arcsec"])
            if not math.isfinite(beam) or beam <= 0.0:
                raise ValueError(f"{source}: invalid geometric beam {selected['nu0_beam_arcsec']}")
            east = _number(east_value, 6)
            north = _number(north_value, 6)
            total = _number(total_value, 6)
            beams = _number(total_value / beam, 6)
            significance = "NOT_COMPUTABLE_COMPONENT_POSITION_UNCERTAINTY_UNAVAILABLE"

        rejected_frame = rejected_ra = rejected_dec = NA
        if policy["rejected_component_record"] != NA:
            rejected_ra_value, rejected_dec_value, _ = _coordinate_to_icrs(
                policy["rejected_component_native_frame"],
                policy["rejected_component_ra_native_deg"],
                policy["rejected_component_dec_native_deg"],
            )
            rejected_frame = COMMON_FRAME
            rejected_ra = _number(rejected_ra_value, 12)
            rejected_dec = _number(rejected_dec_value, 12)

        row = {
            "source": source,
            "selected_flux_status": selected["flux_status"],
            "selected_observation_id": selected["nu0_observation_id"],
            "selected_reference": selected["nu0_reference"],
            "selected_telescope": selected["nu0_telescope"],
            "selected_date": selected["nu0_date"],
            "selected_frequency_ghz": selected["nu0_GHz"],
            "selected_geometric_beam_arcsec": selected["nu0_beam_arcsec"],
            "identity_disposition": disposition,
            "component_identity": identity["component_identity"],
            "identity_evidence_basis": identity["evidence_basis"],
            "required_handling": identity["required_handling"],
            "disposition_evidence": policy["disposition_evidence"],
            "component_record_status": policy["component_record_status"],
            "component_coordinate_method": policy["component_coordinate_method"],
            "component_coordinate_evidence": policy["component_coordinate_evidence"],
            "component_native_frame": policy["component_native_frame"],
            "component_ra_native_deg": policy["component_ra_native_deg"],
            "component_dec_native_deg": policy["component_dec_native_deg"],
            "component_common_frame": component_frame,
            "component_ra_icrs_deg": component_common[0],
            "component_dec_icrs_deg": component_common[1],
            "component_position_uncertainty_status": policy["component_position_uncertainty_status"],
            "component_position_uncertainty_detail": policy["component_position_uncertainty_detail"],
            "reference_coordinate_role": policy["reference_coordinate_role"],
            "reference_coordinate_method": policy["reference_coordinate_method"],
            "reference_coordinate_citation": policy["reference_coordinate_citation"],
            "reference_coordinate_evidence": policy["reference_coordinate_evidence"],
            "reference_native_frame": policy["reference_native_frame"],
            "reference_ra_native_deg": policy["reference_ra_native_deg"],
            "reference_dec_native_deg": policy["reference_dec_native_deg"],
            "reference_common_frame": reference_frame,
            "reference_ra_icrs_deg": reference_common[0],
            "reference_dec_icrs_deg": reference_common[1],
            "reference_position_uncertainty_status": policy["reference_position_uncertainty_status"],
            "reference_position_uncertainty_detail": policy["reference_position_uncertainty_detail"],
            "authoritative_nuclear_coordinate_status": policy["authoritative_nuclear_coordinate_status"],
            "coordinate_transform_method": transform_method,
            "signed_east_offset_arcsec": east,
            "signed_north_offset_arcsec": north,
            "total_separation_arcsec": total,
            "separation_geometric_beams": beams,
            "offset_interpretation": policy["offset_interpretation"],
            "positional_significance_status": significance,
            "selected_image_phase_calibrator_status": policy["selected_image_phase_calibrator_status"],
            "selected_image_phase_calibrator": policy["selected_image_phase_calibrator"],
            "selected_image_phase_calibrator_evidence": policy["selected_image_phase_calibrator_evidence"],
            "selected_image_astrometric_qa_status": policy["selected_image_astrometric_qa_status"],
            "selected_image_astrometric_qa_evidence": policy["selected_image_astrometric_qa_evidence"],
            **{field: wcs[field] for field in WCS_OUTPUT_FIELDS},
            "accepted_component_record": policy["accepted_component_record"],
            "rejected_component_record": policy["rejected_component_record"],
            "rejected_component_native_frame": policy["rejected_component_native_frame"],
            "rejected_component_ra_native_deg": policy["rejected_component_ra_native_deg"],
            "rejected_component_dec_native_deg": policy["rejected_component_dec_native_deg"],
            "rejected_component_common_frame": rejected_frame,
            "rejected_component_ra_icrs_deg": rejected_ra,
            "rejected_component_dec_icrs_deg": rejected_dec,
            "rejected_component_reason": policy["rejected_component_reason"],
        }
        if tuple(row) != FIELDS:
            raise AssertionError(f"internal output field mismatch for {source}")
        output.append(row)
    return output


def render(rows: Iterable[Mapping[str, str]]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    content = render(build_rows())
    if args.check:
        if not args.output.is_file() or args.output.read_text(encoding="utf-8") != content:
            raise SystemExit(f"stale S230 astrometry audit: {args.output}")
        print(f"S230 astrometry audit is current: {EXPECTED_COUNT}/{EXPECTED_COUNT} rows ({args.output})")
        return 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(content, encoding="utf-8")
    print(f"Wrote {EXPECTED_COUNT}/{EXPECTED_COUNT} S230 astrometry rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
