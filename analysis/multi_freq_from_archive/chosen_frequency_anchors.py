#!/usr/bin/env python3
"""Deterministic chosen-frequency anchor selection and consolidation.

The historical continuum fitter generated one candidate CASA model per source
in its ``getchoosercoords`` pass and then delegated consolidation to a notebook
that was not preserved.  This module replaces that hidden step with a small,
CASA-independent manifest and a deterministic selection rule.

Candidate model tables are generated data products.  They are copied into
``images/chosenfreq`` for the final fitter, but they are intentionally not part
of the lightweight source release.  The CSV manifest is the durable numerical
and provenance record.  Generic minimum-beam selection is not a physical-
identity decision: the historical NGC 4945 seed is explicitly rejected and
must be replaced downstream by the guarded source-18 fixed-nucleus rerun.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import shutil
import tempfile
from collections import defaultdict
from pathlib import Path, PurePosixPath
from typing import Iterable, Mapping, Sequence

from continuum_product_additions import (
    load_product_addition_policy,
    publication_addition_products,
)


FIELDS = (
    "source",
    "telescope",
    "frequency_ghz",
    "beam_major_arcsec",
    "beam_minor_arcsec",
    "beam_pa_deg",
    "beam_geometric_arcsec",
    "anchor_frame",
    "anchor_ra_deg",
    "anchor_dec_deg",
    "reviewed_frame",
    "reviewed_ra_deg",
    "reviewed_dec_deg",
    "reviewed_coordinate_basis",
    "component_separation_arcsec",
    "component_tolerance_arcsec",
    "generic_seed_action",
    "icrs_ra",
    "icrs_dec",
    "j2000_ra",
    "j2000_dec",
    "input_product",
    "input_sha256",
    "model_name",
    "model_path",
    "model_tree_sha256",
    "run_context",
    "selection_rule",
    "identity_disposition",
    "component_identity",
    "identity_evidence",
    "required_handling",
    "scientific_disposition",
)

REFERENCE_FIELDS = (
    "source",
    "telescope",
    "frequency_ghz",
    "beam_major_arcsec",
    "beam_minor_arcsec",
    "beam_pa_deg",
    "beam_geometric_arcsec",
    "anchor_frame",
    "anchor_ra_deg",
    "anchor_dec_deg",
    "input_product",
    "input_sha256",
    "model_name",
    "model_tree_sha256",
    "reference_409_row",
    "scientific_status",
    "scientific_note",
)

SELECTION_RULE = (
    "detected_200_400GHz_geometric_beam_lt_1arcsec_then_minimum_geometric_beam"
)

COMPONENT_POLICY_FIELDS = (
    "source",
    "selected_reference",
    "identity_disposition",
    "component_identity",
    "evidence_basis",
    "required_handling",
    "reviewed_frame",
    "reviewed_ra_deg",
    "reviewed_dec_deg",
    "reviewed_coordinate_basis",
    "max_separation_arcsec",
    "generic_seed_action",
)

PENDING_COMPONENT_REVIEW = "pending_component_policy_review"
GENERIC_SEED_ACTIONS = {
    "require_coordinate_match",
    "require_coordinate_match_with_caveat",
    "reject_generic_seed_and_rerun_fixed_reviewed_nucleus",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def tree_sha256(path: Path) -> str:
    """Hash the relative names and bytes of a directory in lexical order."""

    if not path.is_dir():
        raise NotADirectoryError(path)
    digest = hashlib.sha256()
    files = sorted(
        (candidate for candidate in path.rglob("*") if candidate.is_file()),
        key=lambda candidate: candidate.relative_to(path).as_posix(),
    )
    if not files:
        raise ValueError(f"CASA model tree is empty: {path}")
    for item in files:
        relative = item.relative_to(path).as_posix().encode("utf-8")
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        with item.open("rb") as handle:
            while chunk := handle.read(16 * 1024 * 1024):
                digest.update(chunk)
    return digest.hexdigest()


def path_sha256(path: Path) -> str:
    if path.is_file():
        return sha256_file(path)
    return tree_sha256(path)


def _finite_float(row: Mapping[str, object], field: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {field} for anchor candidate: {row!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"Non-finite {field} for anchor candidate: {row!r}")
    return value


def _portable_relative_path(value: object, field: str) -> str:
    text = str(value).strip().replace("\\", "/")
    path = PurePosixPath(text)
    windows_absolute = len(text) >= 3 and text[1:3] == ":/"
    if not text or path.is_absolute() or windows_absolute or ".." in path.parts:
        raise ValueError(f"{field} must be a portable relative path, got {text!r}")
    return path.as_posix()


def _sha256(value: object, field: str) -> str:
    digest = str(value).strip().lower()
    if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
        raise ValueError(f"{field} is not a SHA-256 digest: {digest!r}")
    return digest


def _format_number(value: float) -> str:
    return f"{value:.15g}"


def _normalize_frame(value: object, *, source: str, field: str) -> str:
    frame = str(value).strip().upper()
    if frame in {"FK5", "J2000.0"}:
        frame = "J2000"
    if frame not in {"ICRS", "J2000"}:
        raise ValueError(f"Unsupported {field} for {source}: {frame!r}")
    return frame


def angular_separation_arcsec(
    ra1_deg: float,
    dec1_deg: float,
    ra2_deg: float,
    dec2_deg: float,
) -> float:
    """Return the great-circle separation for two coordinates in degrees."""

    ra1 = math.radians(float(ra1_deg))
    dec1 = math.radians(float(dec1_deg))
    ra2 = math.radians(float(ra2_deg))
    dec2 = math.radians(float(dec2_deg))
    haversine = (
        math.sin((dec2 - dec1) / 2.0) ** 2
        + math.cos(dec1) * math.cos(dec2) * math.sin((ra2 - ra1) / 2.0) ** 2
    )
    haversine = min(1.0, max(0.0, haversine))
    return math.degrees(2.0 * math.asin(math.sqrt(haversine))) * 3600.0


def sexagesimal_ra_to_deg(value: object) -> float:
    text = str(value).strip().lower()
    for character in "hms":
        text = text.replace(character, ":")
    parts = [part for part in text.rstrip(":").split(":") if part]
    if len(parts) != 3:
        raise ValueError(f"Invalid sexagesimal right ascension: {value!r}")
    hours, minutes, seconds = map(float, parts)
    degrees = 15.0 * (hours + minutes / 60.0 + seconds / 3600.0)
    if not 0 <= degrees < 360:
        raise ValueError(f"Right ascension is outside 0--360 degrees: {value!r}")
    return degrees


def sexagesimal_dec_to_deg(value: object) -> float:
    text = str(value).strip().lower()
    for character in "dms":
        text = text.replace(character, ":")
    parts = [part for part in text.rstrip(":").split(":") if part]
    if len(parts) != 3:
        raise ValueError(f"Invalid sexagesimal declination: {value!r}")
    sign = -1.0 if parts[0].startswith("-") else 1.0
    degrees = abs(float(parts[0])) + float(parts[1]) / 60.0 + float(parts[2]) / 3600.0
    degrees *= sign
    if not -90 <= degrees <= 90:
        raise ValueError(f"Declination is outside -90--90 degrees: {value!r}")
    return degrees


def normalize_candidate(row: Mapping[str, object]) -> dict[str, str]:
    """Validate and canonicalize one eligible detected candidate."""

    source = str(row.get("source", "")).strip()
    telescope = str(row.get("telescope", "")).strip().upper()
    if not source:
        raise ValueError("Chosen-frequency candidate has no source")
    if telescope not in {"ALMA", "VLA", "SMA"}:
        raise ValueError(f"Unsupported telescope for {source}: {telescope!r}")

    frequency = _finite_float(row, "frequency_ghz")
    major = _finite_float(row, "beam_major_arcsec")
    minor = _finite_float(row, "beam_minor_arcsec")
    pa = _finite_float(row, "beam_pa_deg")
    geometric = _finite_float(row, "beam_geometric_arcsec")
    expected_geometric = math.sqrt(major * minor) if major > 0 and minor > 0 else -1
    if not 200 <= frequency <= 400:
        raise ValueError(f"Anchor frequency is outside 200--400 GHz for {source}")
    if major <= 0 or minor <= 0 or not geometric < 1:
        raise ValueError(f"Anchor beam is not eligible for {source}")
    if not math.isclose(geometric, expected_geometric, rel_tol=2e-12, abs_tol=1e-15):
        raise ValueError(
            f"Geometric beam mismatch for {source}: {geometric} != {expected_geometric}"
        )

    frame = _normalize_frame(row.get("anchor_frame", ""), source=source, field="anchor frame")
    ra = _finite_float(row, "anchor_ra_deg")
    dec = _finite_float(row, "anchor_dec_deg")
    if not 0 <= ra < 360 or not -90 <= dec <= 90:
        raise ValueError(f"Invalid anchor sky coordinate for {source}: {ra}, {dec}")

    model_name = PurePosixPath(str(row.get("model_name", "")).strip()).name
    if not model_name.endswith(".fit") or model_name != str(row.get("model_name", "")).strip():
        raise ValueError(f"Invalid CASA model name for {source}: {model_name!r}")

    reviewed_values = {
        "reviewed_frame": str(
            row.get("reviewed_frame", PENDING_COMPONENT_REVIEW)
        ).strip(),
        "reviewed_ra_deg": str(
            row.get("reviewed_ra_deg", PENDING_COMPONENT_REVIEW)
        ).strip(),
        "reviewed_dec_deg": str(
            row.get("reviewed_dec_deg", PENDING_COMPONENT_REVIEW)
        ).strip(),
        "reviewed_coordinate_basis": str(
            row.get("reviewed_coordinate_basis", PENDING_COMPONENT_REVIEW)
        ).strip(),
        "component_separation_arcsec": str(
            row.get("component_separation_arcsec", PENDING_COMPONENT_REVIEW)
        ).strip(),
        "component_tolerance_arcsec": str(
            row.get("component_tolerance_arcsec", PENDING_COMPONENT_REVIEW)
        ).strip(),
        "generic_seed_action": str(
            row.get("generic_seed_action", PENDING_COMPONENT_REVIEW)
        ).strip(),
    }
    pending_fields = [
        field
        for field, value in reviewed_values.items()
        if value == PENDING_COMPONENT_REVIEW
    ]
    if pending_fields and len(pending_fields) != len(reviewed_values):
        raise ValueError(
            f"Incomplete component-review metadata for {source}: {pending_fields} pending"
        )
    if not pending_fields:
        reviewed_values["reviewed_frame"] = _normalize_frame(
            reviewed_values["reviewed_frame"], source=source, field="reviewed frame"
        )
        for field in (
            "reviewed_ra_deg",
            "reviewed_dec_deg",
            "component_separation_arcsec",
            "component_tolerance_arcsec",
        ):
            value = _finite_float(reviewed_values, field)
            if field in {"component_separation_arcsec", "component_tolerance_arcsec"}:
                if value < 0 or (field == "component_tolerance_arcsec" and value == 0):
                    raise ValueError(f"Invalid {field} for {source}: {value}")
            reviewed_values[field] = _format_number(value)
        reviewed_ra = float(reviewed_values["reviewed_ra_deg"])
        reviewed_dec = float(reviewed_values["reviewed_dec_deg"])
        if not 0 <= reviewed_ra < 360 or not -90 <= reviewed_dec <= 90:
            raise ValueError(
                f"Invalid reviewed sky coordinate for {source}: "
                f"{reviewed_ra}, {reviewed_dec}"
            )
        if reviewed_values["generic_seed_action"] not in GENERIC_SEED_ACTIONS:
            raise ValueError(
                f"Unsupported generic-seed action for {source}: "
                f"{reviewed_values['generic_seed_action']!r}"
            )

    normalized = {
        "source": source,
        "telescope": telescope,
        "frequency_ghz": _format_number(frequency),
        "beam_major_arcsec": _format_number(major),
        "beam_minor_arcsec": _format_number(minor),
        "beam_pa_deg": _format_number(pa),
        "beam_geometric_arcsec": _format_number(geometric),
        "anchor_frame": frame,
        "anchor_ra_deg": _format_number(ra),
        "anchor_dec_deg": _format_number(dec),
        **reviewed_values,
        "icrs_ra": str(row.get("icrs_ra", "")).strip(),
        "icrs_dec": str(row.get("icrs_dec", "")).strip(),
        "j2000_ra": str(row.get("j2000_ra", "")).strip(),
        "j2000_dec": str(row.get("j2000_dec", "")).strip(),
        "input_product": _portable_relative_path(row.get("input_product"), "input_product"),
        "input_sha256": _sha256(row.get("input_sha256"), "input_sha256"),
        "model_name": model_name,
        "model_path": _portable_relative_path(row.get("model_path"), "model_path"),
        "model_tree_sha256": _sha256(row.get("model_tree_sha256"), "model_tree_sha256"),
        "run_context": str(row.get("run_context", "")).strip(),
        "selection_rule": str(row.get("selection_rule", SELECTION_RULE)).strip(),
        "identity_disposition": str(
            row.get("identity_disposition", PENDING_COMPONENT_REVIEW)
        ).strip(),
        "component_identity": str(
            row.get("component_identity", PENDING_COMPONENT_REVIEW)
        ).strip(),
        "identity_evidence": str(
            row.get("identity_evidence", PENDING_COMPONENT_REVIEW)
        ).strip(),
        "required_handling": str(
            row.get("required_handling", PENDING_COMPONENT_REVIEW)
        ).strip(),
        "scientific_disposition": str(
            row.get("scientific_disposition", PENDING_COMPONENT_REVIEW)
        ).strip(),
    }
    if normalized["selection_rule"] != SELECTION_RULE:
        raise ValueError(
            f"Unrecognized chosen-frequency selection rule for {source}: "
            f"{normalized['selection_rule']!r}"
        )
    for field in (
        "identity_disposition",
        "component_identity",
        "identity_evidence",
        "required_handling",
        "scientific_disposition",
    ):
        if not normalized[field]:
            raise ValueError(f"Missing {field} for {source}")
    return normalized


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        missing = set(FIELDS).difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} is missing fields: {sorted(missing)}")
        rows = [normalize_candidate(row) for row in reader]
    if not rows:
        raise ValueError(f"Chosen-frequency manifest is empty: {path}")
    return rows


def write_manifest(path: Path, rows: Iterable[Mapping[str, object]]) -> None:
    normalized = [normalize_candidate(row) for row in rows]
    if not normalized:
        raise ValueError("Refusing to write an empty chosen-frequency manifest")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            newline="",
            prefix=f".{path.name}.",
            suffix=".tmp",
            dir=path.parent,
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
            writer.writeheader()
            writer.writerows(sorted(normalized, key=lambda row: row["source"]))
            handle.flush()
            os.fsync(handle.fileno())
        temporary_path.chmod(0o644)
        os.replace(temporary_path, path)
        temporary_path = None
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def select_anchors(rows: Iterable[Mapping[str, object]]) -> list[dict[str, str]]:
    """Select one candidate per source with a deterministic tie-break."""

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        normalized = normalize_candidate(row)
        grouped[normalized["source"]].append(normalized)
    if not grouped:
        raise ValueError("No chosen-frequency candidates were supplied")

    selected = []
    for source in sorted(grouped):
        winner = min(
            grouped[source],
            key=lambda row: (
                float(row["beam_geometric_arcsec"]),
                row["input_product"],
                row["telescope"],
                row["model_name"],
                row["run_context"],
            ),
        )
        selected.append(winner)
    return selected


def read_component_policy(path: Path) -> dict[str, dict[str, str]]:
    """Read the reviewed component-identity decisions for selected products."""

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        missing = set(COMPONENT_POLICY_FIELDS).difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} is missing component-policy fields: {sorted(missing)}")
        rows = [dict(row) for row in reader]
    policy: dict[str, dict[str, str]] = {}
    for row in rows:
        source = row["source"].strip()
        if not source or source in policy:
            raise ValueError(f"Duplicate or empty component-policy source: {source!r}")
        row["selected_reference"] = _portable_relative_path(
            row["selected_reference"], "selected_reference"
        )
        for field in (
            "identity_disposition",
            "component_identity",
            "evidence_basis",
            "required_handling",
            "reviewed_coordinate_basis",
            "generic_seed_action",
        ):
            row[field] = row[field].strip()
            if not row[field]:
                raise ValueError(f"Empty {field} in component policy for {source}")
        row["reviewed_frame"] = _normalize_frame(
            row["reviewed_frame"], source=source, field="reviewed frame"
        )
        for field in ("reviewed_ra_deg", "reviewed_dec_deg", "max_separation_arcsec"):
            row[field] = _format_number(_finite_float(row, field))
        reviewed_ra = float(row["reviewed_ra_deg"])
        reviewed_dec = float(row["reviewed_dec_deg"])
        tolerance = float(row["max_separation_arcsec"])
        if not 0 <= reviewed_ra < 360 or not -90 <= reviewed_dec <= 90:
            raise ValueError(
                f"Invalid reviewed sky coordinate for {source}: "
                f"{reviewed_ra}, {reviewed_dec}"
            )
        if tolerance <= 0:
            raise ValueError(
                f"Component-policy tolerance must be positive for {source}: {tolerance}"
            )
        if row["generic_seed_action"] not in GENERIC_SEED_ACTIONS:
            raise ValueError(
                f"Unsupported generic-seed action for {source}: "
                f"{row['generic_seed_action']!r}"
            )
        policy[source] = row
    if not policy:
        raise ValueError(f"Component policy is empty: {path}")
    return policy


def review_selected_anchors(
    rows: Iterable[Mapping[str, object]],
    component_policy: Mapping[str, Mapping[str, str]],
) -> list[dict[str, str]]:
    """Fail closed unless each generic winner matches a reviewed component.

    This check is deliberately performed *after* the generic minimum-beam
    selection.  A smaller-beam image of a brighter nearby component therefore
    cannot silently become the nuclear anchor: a product/component mismatch is
    a hard error requiring a reviewed policy change.
    """

    reviewed = []
    reviewed_sources = set()
    for raw_row in rows:
        row = normalize_candidate(raw_row)
        source = row["source"]
        if source in reviewed_sources:
            raise ValueError(f"Duplicate selected anchor source: {source}")
        reviewed_sources.add(source)
        policy = component_policy.get(source)
        if policy is None:
            raise ValueError(f"No reviewed nuclear/component disposition for {source}")
        expected_product = str(policy["selected_reference"]).strip().replace("\\", "/")
        if row["input_product"] != expected_product:
            raise ValueError(
                f"Reviewed nuclear product mismatch for {source}: generic minimum-beam "
                f"winner {row['input_product']!r}, reviewed product {expected_product!r}"
            )
        disposition = str(policy["identity_disposition"]).strip()
        handling = str(policy["required_handling"]).strip()
        action = str(policy["generic_seed_action"]).strip()
        reviewed_frame = str(policy["reviewed_frame"]).strip()
        reviewed_ra = float(policy["reviewed_ra_deg"])
        reviewed_dec = float(policy["reviewed_dec_deg"])
        tolerance = float(policy["max_separation_arcsec"])
        if row["anchor_frame"] != reviewed_frame:
            raise ValueError(
                f"Reviewed component frame mismatch for {source}: fitted anchor uses "
                f"{row['anchor_frame']}, policy uses {reviewed_frame}"
            )
        separation = angular_separation_arcsec(
            float(row["anchor_ra_deg"]),
            float(row["anchor_dec_deg"]),
            reviewed_ra,
            reviewed_dec,
        )
        if disposition not in {"nuclear", "ambiguous"}:
            raise ValueError(
                f"Chosen-frequency detection for {source} has non-usable reviewed "
                f"disposition {disposition!r}"
            )
        if source == "NGC4945":
            if handling != "use_source18_fixed_photometry_exclude_source17":
                raise ValueError("NGC 4945 policy does not require the source-18 correction")
            if action != "reject_generic_seed_and_rerun_fixed_reviewed_nucleus":
                raise ValueError(
                    "NGC 4945 policy does not reject the generic source-17 seed"
                )
            scientific_disposition = (
                "generic_source17_seed_rejected_source18_fixed_nucleus_rerun_required"
            )
        elif disposition == "ambiguous":
            if handling != "retain_with_ambiguity_caveat":
                raise ValueError(f"Ambiguous anchor lacks caveat handling for {source}")
            if action != "require_coordinate_match_with_caveat":
                raise ValueError(
                    f"Ambiguous anchor has unexpected generic-seed action for {source}"
                )
            if separation > tolerance:
                raise ValueError(
                    f"Reviewed component coordinate mismatch for {source}: fitted anchor "
                    f"is {separation:.6g} arcsec from the reviewed position "
                    f"(maximum {tolerance:.6g} arcsec)"
                )
            scientific_disposition = "reviewed_ambiguous_anchor_retain_with_caveat"
        else:
            if handling != "use_selected_observation":
                raise ValueError(f"Nuclear anchor has unexpected handling for {source}")
            if action != "require_coordinate_match":
                raise ValueError(
                    f"Nuclear anchor has unexpected generic-seed action for {source}"
                )
            if separation > tolerance:
                raise ValueError(
                    f"Reviewed component coordinate mismatch for {source}: fitted anchor "
                    f"is {separation:.6g} arcsec from the reviewed position "
                    f"(maximum {tolerance:.6g} arcsec)"
                )
            scientific_disposition = "reviewed_nuclear_anchor"
        row.update(
            {
                "reviewed_frame": reviewed_frame,
                "reviewed_ra_deg": _format_number(reviewed_ra),
                "reviewed_dec_deg": _format_number(reviewed_dec),
                "reviewed_coordinate_basis": str(
                    policy["reviewed_coordinate_basis"]
                ).strip(),
                "component_separation_arcsec": _format_number(separation),
                "component_tolerance_arcsec": _format_number(tolerance),
                "generic_seed_action": action,
                "identity_disposition": disposition,
                "component_identity": str(policy["component_identity"]).strip(),
                "identity_evidence": str(policy["evidence_basis"]).strip(),
                "required_handling": handling,
                "scientific_disposition": scientific_disposition,
            }
        )
        reviewed.append(row)
    unused_policy = sorted(set(component_policy).difference(reviewed_sources))
    if unused_policy:
        raise ValueError(
            "No eligible chosen-frequency candidate for reviewed source(s): "
            + ", ".join(unused_policy)
        )
    return reviewed


def _inside(root: Path, relative: str) -> Path:
    candidate = (root / relative).resolve()
    try:
        candidate.relative_to(root.resolve())
    except ValueError as exc:
        raise ValueError(f"Path escapes the analysis workspace: {relative!r}") from exc
    return candidate


def materialize_models(
    rows: Sequence[Mapping[str, object]],
    *,
    analysis_root: Path,
    output_dir: Path,
    replace: bool = False,
) -> list[dict[str, str]]:
    """Copy selected CASA models atomically into the final fitter directory."""

    normalized = [normalize_candidate(row) for row in rows]
    analysis_root = analysis_root.resolve()
    output_dir = output_dir.resolve()
    try:
        output_dir.relative_to(analysis_root)
    except ValueError as exc:
        raise ValueError("Chosen-model output must be inside the analysis workspace") from exc
    if output_dir.is_symlink():
        raise ValueError(f"Refusing to replace symlinked chosen-model directory: {output_dir}")
    if output_dir.exists() and not replace:
        raise FileExistsError(output_dir)

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".chosenfreq.", dir=output_dir.parent))
    completed: list[dict[str, str]] = []
    try:
        for row in normalized:
            source_model = _inside(analysis_root, row["model_path"])
            if not source_model.is_dir() or source_model.is_symlink():
                raise ValueError(f"Candidate CASA model is not an isolated directory: {source_model}")
            actual_hash = tree_sha256(source_model)
            if actual_hash != row["model_tree_sha256"]:
                raise ValueError(
                    f"CASA model hash mismatch for {row['source']}: "
                    f"{actual_hash} != {row['model_tree_sha256']}"
                )
            destination = temporary / row["source"] / row["model_name"]
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(source_model, destination)
            if tree_sha256(destination) != actual_hash:
                raise ValueError(f"Copied CASA model hash mismatch for {row['source']}")
            updated = dict(row)
            updated["model_path"] = destination.relative_to(temporary).as_posix()
            completed.append(updated)

        backup = None
        if output_dir.exists():
            backup = Path(tempfile.mkdtemp(prefix=".chosenfreq-old.", dir=output_dir.parent))
            backup.rmdir()
            output_dir.rename(backup)
        try:
            temporary.rename(output_dir)
        except Exception:
            if backup is not None and backup.exists() and not output_dir.exists():
                backup.rename(output_dir)
            raise
        if backup is not None:
            shutil.rmtree(backup)

        canonical = []
        for row in completed:
            updated = dict(row)
            updated["model_path"] = (
                output_dir.relative_to(analysis_root) / row["model_path"]
            ).as_posix()
            canonical.append(updated)
        return canonical
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def consolidate(
    candidate_manifests: Sequence[Path],
    *,
    analysis_root: Path,
    output_manifest: Path,
    models_output: Path,
    component_policy: Path,
    replace_models: bool = False,
) -> list[dict[str, str]]:
    candidates: list[dict[str, str]] = []
    for manifest in candidate_manifests:
        candidates.extend(read_manifest(manifest))
    winners = review_selected_anchors(
        select_anchors(candidates), read_component_policy(component_policy)
    )
    materialized = materialize_models(
        winners,
        analysis_root=analysis_root,
        output_dir=models_output,
        replace=replace_models,
    )
    write_manifest(output_manifest, materialized)
    return materialized


def _reference_identity(row: Mapping[str, object]) -> tuple[str, str, str]:
    return (
        str(row["telescope"]).strip().upper(),
        str(row["input_product"]).strip().replace("\\", "/"),
        str(row["model_name"]).strip(),
    )


def compare_reference_identities(
    actual_rows: Iterable[Mapping[str, object]],
    reference_rows: Iterable[Mapping[str, object]],
    *,
    reviewed_addition_products: set[str] | None = None,
    component_policy: Mapping[str, Mapping[str, str]] | None = None,
) -> list[str]:
    """Preserve nine historical identities and allow only reviewed new anchors."""

    actual = {str(row["source"]).strip(): _reference_identity(row) for row in actual_rows}
    reference = {
        str(row["source"]).strip(): _reference_identity(row) for row in reference_rows
    }
    problems = []
    admitted = reviewed_addition_products or set()
    for source in sorted(reference):
        if source not in actual:
            problems.append(f"missing actual anchor: {source}")
        elif actual[source] != reference[source]:
            product = actual[source][1]
            policy = component_policy.get(source) if component_policy is not None else None
            if product not in admitted or policy is None or policy.get("selected_reference") != product:
                problems.append(
                    f"identity mismatch for {source}: actual={actual[source]!r}, "
                    f"reference={reference[source]!r}"
                )
    for source in sorted(set(actual) - set(reference)):
        product = actual[source][1]
        policy = component_policy.get(source) if component_policy is not None else None
        if product not in admitted:
            problems.append(f"unexpected actual anchor: {source}")
        elif policy is None or policy.get("selected_reference") != product:
            problems.append(
                f"reviewed addition anchor lacks matching component policy: {source}"
            )
    return problems


def validate_reference_artifacts(
    reference_rows: Iterable[Mapping[str, object]],
    *,
    archive_analysis_root: Path,
    archive_chosen_root: Path,
) -> None:
    """Validate the preserved input and CASA-model bytes behind the reference."""

    rows = list(reference_rows)
    if len(rows) != 9:
        raise ValueError(f"Expected nine preserved historical anchors, found {len(rows)}")
    for row in rows:
        source = str(row["source"]).strip()
        input_path = archive_analysis_root / str(row["input_product"]).strip()
        model_path = archive_chosen_root / source / str(row["model_name"]).strip()
        if not input_path.is_file():
            raise FileNotFoundError(input_path)
        if sha256_file(input_path) != str(row["input_sha256"]).strip():
            raise ValueError(f"Historical input hash mismatch for {source}")
        if not model_path.is_dir():
            raise FileNotFoundError(model_path)
        if tree_sha256(model_path) != str(row["model_tree_sha256"]).strip():
            raise ValueError(f"Historical model-tree hash mismatch for {source}")


def read_reference(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        missing = set(REFERENCE_FIELDS).difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} is missing reference fields: {sorted(missing)}")
        rows = [dict(row) for row in reader]
    if len(rows) != 9:
        raise ValueError(f"Expected nine historical anchor identities in {path}")
    sources = [row.get("source", "").strip() for row in rows]
    if any(not source for source in sources) or len(set(sources)) != len(sources):
        raise ValueError(f"Historical anchor sources are empty or duplicated in {path}")
    for row in rows:
        source = row["source"].strip()
        row["input_product"] = _portable_relative_path(
            row["input_product"], "input_product"
        )
        row["input_sha256"] = _sha256(row["input_sha256"], "input_sha256")
        row["model_tree_sha256"] = _sha256(
            row["model_tree_sha256"], "model_tree_sha256"
        )
        row["anchor_frame"] = _normalize_frame(
            row["anchor_frame"], source=source, field="anchor frame"
        )
        numbers = {}
        for field in (
            "frequency_ghz",
            "beam_major_arcsec",
            "beam_minor_arcsec",
            "beam_pa_deg",
            "beam_geometric_arcsec",
            "anchor_ra_deg",
            "anchor_dec_deg",
        ):
            numbers[field] = _finite_float(row, field)
        if not 200 <= numbers["frequency_ghz"] <= 400:
            raise ValueError(f"Historical anchor frequency is invalid for {source}")
        if (
            numbers["beam_major_arcsec"] <= 0
            or numbers["beam_minor_arcsec"] <= 0
            or not 0 < numbers["beam_geometric_arcsec"] < 1
            or not math.isclose(
                numbers["beam_geometric_arcsec"],
                math.sqrt(
                    numbers["beam_major_arcsec"] * numbers["beam_minor_arcsec"]
                ),
                rel_tol=2e-12,
                abs_tol=1e-15,
            )
        ):
            raise ValueError(f"Historical anchor beam is invalid for {source}")
        if not 0 <= numbers["anchor_ra_deg"] < 360 or not -90 <= numbers[
            "anchor_dec_deg"
        ] <= 90:
            raise ValueError(f"Historical anchor coordinate is invalid for {source}")
        model_name = str(row["model_name"]).strip()
        if PurePosixPath(model_name).name != model_name or not model_name.endswith(".fit"):
            raise ValueError(f"Historical CASA model name is invalid for {source}")
        try:
            reference_row = int(row["reference_409_row"])
        except ValueError as exc:
            raise ValueError(f"Historical 409-row link is invalid for {source}") from exc
        if not 1 <= reference_row <= 409 or not row["scientific_note"].strip():
            raise ValueError(f"Historical reference metadata is invalid for {source}")
    if len({row["reference_409_row"] for row in rows}) != len(rows):
        raise ValueError(f"Historical 409-row links are duplicated in {path}")
    ngc4945 = [row for row in rows if row.get("source", "").strip() == "NGC4945"]
    if len(ngc4945) != 1 or ngc4945[0].get("scientific_status") != "rejected_component_seed":
        raise ValueError(
            "The historical NGC 4945 source-17 anchor must be explicitly marked "
            "as a rejected component seed"
        )
    for row in rows:
        if row["source"].strip() != "NGC4945" and row["scientific_status"] != (
            "historical_identity_reference"
        ):
            raise ValueError(
                f"Unexpected historical scientific status for {row['source']}: "
                f"{row['scientific_status']!r}"
            )
    return rows


def check_materialized_manifest(path: Path, analysis_root: Path) -> list[dict[str, str]]:
    rows = read_manifest(path)
    seen = set()
    for row in rows:
        if row["source"] in seen:
            raise ValueError(f"Duplicate selected anchor source: {row['source']}")
        seen.add(row["source"])
        model = _inside(analysis_root.resolve(), row["model_path"])
        if not model.is_dir() or model.is_symlink():
            raise ValueError(f"Missing isolated chosen model for {row['source']}: {model}")
        if tree_sha256(model) != row["model_tree_sha256"]:
            raise ValueError(f"Chosen model hash mismatch for {row['source']}")
        for field in (
            "reviewed_frame",
            "reviewed_ra_deg",
            "reviewed_dec_deg",
            "reviewed_coordinate_basis",
            "component_separation_arcsec",
            "component_tolerance_arcsec",
            "generic_seed_action",
            "identity_disposition",
            "component_identity",
            "identity_evidence",
            "required_handling",
            "scientific_disposition",
        ):
            if PENDING_COMPONENT_REVIEW in row[field]:
                raise ValueError(
                    f"Chosen anchor lacks reviewed component policy for {row['source']}"
                )
        if row["anchor_frame"] != row["reviewed_frame"]:
            raise ValueError(
                f"Chosen anchor/reviewed frame mismatch for {row['source']}"
            )
        recomputed_separation = angular_separation_arcsec(
            float(row["anchor_ra_deg"]),
            float(row["anchor_dec_deg"]),
            float(row["reviewed_ra_deg"]),
            float(row["reviewed_dec_deg"]),
        )
        recorded_separation = float(row["component_separation_arcsec"])
        if not math.isclose(
            recomputed_separation,
            recorded_separation,
            rel_tol=1e-10,
            abs_tol=1e-9,
        ):
            raise ValueError(
                f"Stale component separation for {row['source']}: "
                f"{recorded_separation} != {recomputed_separation} arcsec"
            )
        tolerance = float(row["component_tolerance_arcsec"])
        if row["source"] == "NGC4945" and (
            row["required_handling"]
            != "use_source18_fixed_photometry_exclude_source17"
            or row["generic_seed_action"]
            != "reject_generic_seed_and_rerun_fixed_reviewed_nucleus"
            or "source18_fixed_nucleus_rerun_required"
            not in row["scientific_disposition"]
        ):
            raise ValueError("NGC 4945 chosen seed is not guarded by source-18 handling")
        if row["source"] != "NGC4945":
            if recorded_separation > tolerance:
                raise ValueError(
                    f"Chosen anchor fails reviewed component tolerance for {row['source']}"
                )
            expected_action = (
                "require_coordinate_match_with_caveat"
                if row["identity_disposition"] == "ambiguous"
                else "require_coordinate_match"
            )
            if row["generic_seed_action"] != expected_action:
                raise ValueError(
                    f"Chosen anchor has inconsistent generic-seed action for {row['source']}"
                )
    return rows


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    consolidate_parser = subparsers.add_parser("consolidate")
    consolidate_parser.add_argument("--candidate", action="append", type=Path, required=True)
    consolidate_parser.add_argument("--analysis-root", type=Path, default=Path("."))
    consolidate_parser.add_argument(
        "--output-manifest",
        type=Path,
        default=Path("chosen_frequency/chosen_frequency_anchors.csv"),
    )
    consolidate_parser.add_argument(
        "--models-output", type=Path, default=Path("images/chosenfreq")
    )
    consolidate_parser.add_argument("--replace-models", action="store_true")
    consolidate_parser.add_argument("--component-policy", type=Path, required=True)

    check_parser = subparsers.add_parser("check")
    check_parser.add_argument(
        "--manifest", type=Path, default=Path("chosen_frequency/chosen_frequency_anchors.csv")
    )
    check_parser.add_argument("--analysis-root", type=Path, default=Path("."))
    check_parser.add_argument("--reference", type=Path)
    check_parser.add_argument("--addition-policy", type=Path)
    check_parser.add_argument("--component-policy", type=Path)
    reference_parser = subparsers.add_parser("check-reference-artifacts")
    reference_parser.add_argument("--reference", type=Path, required=True)
    reference_parser.add_argument("--archive-analysis-root", type=Path, required=True)
    reference_parser.add_argument("--archive-chosen-root", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.command == "consolidate":
        rows = consolidate(
            args.candidate,
            analysis_root=args.analysis_root,
            output_manifest=args.output_manifest,
            models_output=args.models_output,
            component_policy=args.component_policy,
            replace_models=args.replace_models,
        )
        print(f"selected {len(rows)} chosen-frequency anchors")
        return 0

    if args.command == "check-reference-artifacts":
        rows = read_reference(args.reference)
        validate_reference_artifacts(
            rows,
            archive_analysis_root=args.archive_analysis_root,
            archive_chosen_root=args.archive_chosen_root,
        )
        print("validated nine preserved historical anchor input/model artifacts")
        return 0

    rows = check_materialized_manifest(args.manifest, args.analysis_root)
    if args.reference:
        admitted: set[str] = set()
        component_policy = None
        if args.addition_policy is not None:
            admitted = publication_addition_products(
                load_product_addition_policy(args.addition_policy)
            )
            if args.component_policy is None:
                raise ValueError("--addition-policy requires --component-policy")
            component_policy = read_component_policy(args.component_policy)
        problems = compare_reference_identities(
            rows,
            read_reference(args.reference),
            reviewed_addition_products=admitted,
            component_policy=component_policy,
        )
        if problems:
            raise ValueError("; ".join(problems))
    print(f"validated {len(rows)} chosen-frequency anchors")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
