#!/usr/bin/env python3
"""Run the four historical continuum-fit modes in an isolated workspace.

This orchestrator never runs against a Git worktree.  Its ``prepare`` command
creates a new isolated copy whose archive FITS inputs are independent files by
default and whose SMA CASA tables are real directory copies.  It refuses an
existing destination and never invokes CASA.  The ``run`` command preserves
one log and one output-hash snapshot per stage.  After an intentional
``review_required`` stop, ``adjudicate`` checks a hash-frozen, phase-specific
delta ledger without invoking CASA or modifying either continuum table.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

from continuum_classification_identity import (
    validate_classification_identity_policy,
)
from continuum_measurements import portable_release_product
from continuum_noise_inventory import validate_final_noise_inventory
from frozen_fitter_inputs import (
    load_frozen_fitter_inputs,
    validate_and_select_frozen_inputs,
)
from ngc4945_bootstrap import (
    derive_bootstrap as derive_source18_bootstrap,
    initialize_bootstrap as initialize_source18_bootstrap,
    validate_bootstrap_snapshot as validate_source18_bootstrap_snapshot,
    validate_preapplication_lineage,
)
from continuum_boundary_rejections import (
    boundary_rejection_products,
    load_boundary_rejection_policy,
)
from continuum_product_additions import (
    ADDITION_POLICY_FIELDS,
    addition_products,
    fitter_output_source,
    load_product_addition_policy,
    publication_addition_products,
    publication_addition_products_in_fitter_order,
)
from rerun_ngc4945_nuclear import (
    COMPONENT_FIELDS,
    CORRECTION_FIELDS,
    CORRECTION_SCHEMA_VERSION,
    LEGACY_CORRECTION_FIELDS,
    PHOTOMETRY_FIELDS,
    EXPECTED_ROWS as EXPECTED_NGC4945_ROWS,
    read_exact_csv,
    validate_component_records,
    validate_photometry_records,
)


MARKER_NAME = ".continuum-fit-isolated-workspace.json"
PREPARED_INPUT_SCHEMA_VERSION = 5
STAGES = (
    ("prelimchecker", "prelim"),
    ("fit_sizenoloc", "noloc"),
    ("getchoosercoords", "getchooser"),
    ("fixtochosenfreq", "fixedtochosen"),
)
WORKFLOW_INPUT_NAMES = (
    "run_continuum_four_pass.py",
    "get_continuum_params.py",
    "chosen_frequency_anchors.py",
    "chosen_frequency_component_policy.csv",
    "chosen_frequency_anchor_reference.csv",
    "continuum_product_exclusion_policy.csv",
    "continuum_boundary_rejections.py",
    "continuum_boundary_rejection_policy.csv",
    "continuum_boundary_rejection_evidence.csv",
    "continuum_product_additions.py",
    "continuum_product_addition_policy.csv",
    "continuum_measurements.py",
    "continuum_output_paths.py",
    "continuum_units.py",
    "continuum_noise_inventory.py",
    "continuum_classification_identity.py",
    "continuum_classification_identity_policy.csv",
    "celestial_frames.py",
    "selected_clean_manifest.py",
    "namesandcoords.txt",
    "rerun_ngc4945_nuclear.py",
    "ngc4945_bootstrap.py",
    "source18_initial_history/ngc4945_nuclear_photometry.csv",
    "source18_initial_history/validated_continuum_corrections.csv",
    "source18_initial_history/public_source_provenance.json",
    "frozen_fitter_inputs.py",
    "frozen_fitter_inputs.csv",
)
EXPECTED_ARCHIVE_INPUTS = {"ALMA": 133, "VLA": 307}
EXPECTED_SMA_DISCOVERY_INPUTS = 15
EXPECTED_CANONICAL_ROWS = 409
EXPECTED_REVIEWED_EXCLUSIONS = 3
EXPECTED_CORRECTED_ROWS = 406
HISTORICAL_CANONICAL_RELATIVE = Path(
    "fitsumfiles/historical/fitsummary_fixtochosen_sorted_409_reference.txt"
)

SORTED_REPORTED_FIELD_NAMES = (
    "source",
    "telescope",
    "frequency_ghz",
    "flux_field",
    "date_label",
    "beam_arcsec",
    "snr",
)
REVIEW_SCOPE = "final_local_source_excluded_scaled_mad"
REVIEW_STATUS = "frozen_after_exact_inventory_delta_and_identity_review"
REVIEW_DELTA_FIELDS = (
    "corrected_row_number",
    "source",
    "telescope",
    "frequency_ghz",
    "date_label",
    "product",
    "current_flux_field",
    "generated_flux_field",
    "current_snr",
    "generated_snr",
    "changed_fields",
    "review_reason",
)
ADDITION_REVIEW_SCOPE = "reviewed_new_product_admission"
REVIEW_ADDITION_FIELDS = (
    "admission_id",
    "source",
    "telescope",
    "product",
    "product_sha256",
    "expected_run_disposition",
    "observed_run_disposition",
    "generated_row_number",
    "frequency_ghz",
    "flux_field",
    "date_label",
    "beam_arcsec",
    "snr",
    "review_reason",
)
REVIEW_FREEZE_FIELDS = (
    "schema_version",
    "review_scope",
    "review_status",
    "run_id",
    "expected_counts",
    "artifacts",
)
REVIEW_FREEZE_COUNT_FIELDS = (
    "canonical_rows",
    "reviewed_exclusion_rows",
    "generated_rows",
    "generic_noise_audit_rows",
    "generic_publication_noise_audit_rows",
    "boundary_rejection_policy_rows",
    "boundary_rejected_noise_audit_rows",
    "ngc4945_photometry_rows",
    "ngc4945_component_rows",
    "preapplication_correction_rows",
    "postrun_correction_rows",
    "classification_identity_rows",
    "reviewed_delta_rows",
    "addition_policy_rows",
    "publication_addition_rows",
    "reviewed_addition_rows",
)
REVIEW_FREEZE_ARTIFACT_NAMES = (
    "workflow_plan",
    "prepared_archive_inventory",
    "prepared_input_integrity",
    "initial_comparison",
    "final_stage_snapshot",
    "canonical_table",
    "generated_table",
    "exclusion_policy",
    "addition_policy",
    "boundary_rejection_policy",
    "generic_final_noise_audit",
    "ngc4945_nuclear_photometry",
    "ngc4945_component_audit",
    "final_noise_inventory_validation",
    "correction_manifest_transition_validation",
    "classification_identity_policy",
    "classification_identity_validation",
    "analysis_utils_snapshot",
    "preapplication_cumulative_corrections",
    "postrun_cumulative_corrections",
    "reviewed_delta_csv",
    "reviewed_addition_csv",
)
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def tree_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    entries = sorted(
        path.rglob("*"), key=lambda item: item.relative_to(path).as_posix()
    )
    for item in entries:
        relative = item.relative_to(path).as_posix().encode("utf-8")
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        if item.is_symlink():
            target = os.readlink(item).encode("utf-8")
            digest.update(b"L" + len(target).to_bytes(8, "big") + target)
        elif item.is_file():
            digest.update(b"F")
            with item.open("rb") as handle:
                while chunk := handle.read(16 * 1024 * 1024):
                    digest.update(chunk)
        elif item.is_dir():
            digest.update(b"D")
    return digest.hexdigest()


def analysis_utils_bundle_inventory(path: Path) -> dict[str, object]:
    """Hash the exact importable analysisUtils dependency bundle.

    Python bytecode caches and CASA log files are execution by-products, not
    inputs.  Every other regular file is frozen; symlinks are rejected so the
    evidence directory is self-contained.
    """

    if path.is_symlink():
        raise ValueError(f"analysisUtils bundle must not be a symlink: {path}")
    path = path.resolve()
    if not path.is_dir():
        raise ValueError(f"analysisUtils bundle must be a real directory: {path}")
    files: dict[str, dict[str, object]] = {}
    for item in sorted(path.rglob("*"), key=lambda value: value.relative_to(path).as_posix()):
        relative = item.relative_to(path)
        if "__pycache__" in relative.parts or item.name.startswith("casa-"):
            continue
        if item.suffix == ".pyc":
            continue
        if item.is_symlink():
            raise ValueError(f"analysisUtils bundle contains a symlink: {item}")
        if item.is_dir():
            continue
        if not item.is_file():
            raise ValueError(f"analysisUtils bundle contains a non-file: {item}")
        files[relative.as_posix()] = {
            "size": item.stat().st_size,
            "sha256": sha256_file(item),
        }
    if "analysisUtils.py" not in files or "mpfit.py" not in files:
        raise ValueError("analysisUtils bundle lacks required analysisUtils.py/mpfit.py")
    serialized = json.dumps(files, sort_keys=True, separators=(",", ":")).encode(
        "utf-8"
    )
    return {
        "schema_version": 1,
        "file_count": len(files),
        "file_inventory_sha256": hashlib.sha256(serialized).hexdigest(),
        "files": files,
    }


def freeze_analysis_utils_bundle(
    source: Path,
    destination: Path,
    manifest_path: Path,
) -> dict[str, object]:
    """Copy and verify the complete import bundle used by every CASA stage."""

    if source.is_symlink():
        raise ValueError(f"analysisUtils source must not be a symlink: {source}")
    source = source.resolve()
    if not source.is_dir():
        raise ValueError(f"analysisUtils source must be a real directory: {source}")
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    if manifest_path.exists() or manifest_path.is_symlink():
        raise FileExistsError(manifest_path)
    source_inventory = analysis_utils_bundle_inventory(source)
    shutil.copytree(
        source,
        destination,
        symlinks=False,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "casa-*.log"),
    )
    copied_inventory = analysis_utils_bundle_inventory(destination)
    if copied_inventory != source_inventory:
        raise RuntimeError("analysisUtils dependency bundle changed while being copied")
    manifest = {
        **copied_inventory,
        "source_directory": str(source),
        "execution_directory": str(destination.resolve()),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def _analysis_utils_inventory_from_manifest(
    manifest: dict[str, object],
) -> dict[str, object]:
    """Return and schema-check the immutable inventory portion of a manifest."""

    expected_fields = {
        "schema_version",
        "file_count",
        "file_inventory_sha256",
        "files",
        "source_directory",
        "execution_directory",
    }
    if set(manifest) != expected_fields:
        raise ValueError("analysisUtils bundle manifest has an unexpected schema")
    inventory = {
        key: manifest[key]
        for key in (
            "schema_version",
            "file_count",
            "file_inventory_sha256",
            "files",
        )
    }
    if inventory["schema_version"] != 1:
        raise ValueError("analysisUtils bundle manifest has a stale schema version")
    if not isinstance(inventory["files"], dict):
        raise ValueError("analysisUtils bundle manifest lacks a file inventory")
    if inventory["file_count"] != len(inventory["files"]):
        raise ValueError("analysisUtils bundle manifest has a stale file count")
    digest = inventory["file_inventory_sha256"]
    if not isinstance(digest, str) or not SHA256_PATTERN.fullmatch(digest):
        raise ValueError("analysisUtils bundle manifest has an invalid inventory hash")
    return inventory


def validate_ngc4945_component_artifacts(
    photometry_path: Path,
    component_path: Path,
) -> dict[str, object]:
    """Recompute the semantic source-18/component contract for frozen evidence."""

    photometry = read_exact_csv(
        photometry_path,
        PHOTOMETRY_FIELDS,
        label="NGC 4945 photometry sidecar",
    )
    validate_photometry_records(photometry)
    components = read_exact_csv(
        component_path,
        COMPONENT_FIELDS,
        label="NGC 4945 component audit",
    )
    validate_component_records(components, photometry)
    return {
        "ngc4945_photometry_semantic_status": "exact_validated",
        "ngc4945_photometry_semantic_rows": len(photometry),
        "ngc4945_component_semantic_status": "exact_validated",
        "ngc4945_component_semantic_rows": len(components),
        "ngc4945_component_record_ids": sorted(
            str(record["record_id"]) for record in components
        ),
    }


def repository_root_for_script() -> Path:
    return Path(__file__).resolve().parents[2]


def analysis_root(workspace: Path) -> Path:
    return workspace / "analysis/multi_freq_from_archive"


def assert_isolated_layout(workspace: Path, *, require_marker: bool) -> Path:
    workspace = workspace.resolve()
    if not workspace.is_dir():
        raise FileNotFoundError(workspace)
    if (workspace / ".git").exists():
        raise ValueError(
            f"Refusing to run continuum fits in a Git worktree: {workspace}"
        )
    fitter = analysis_root(workspace) / "get_continuum_params.py"
    if not fitter.is_file():
        raise FileNotFoundError(fitter)
    marker = workspace / MARKER_NAME
    if require_marker:
        if not marker.is_file():
            raise ValueError(
                f"Workspace has not been explicitly initialized: missing {marker}"
            )
        payload = json.loads(marker.read_text())
        if Path(payload.get("workspace", "")).resolve() != workspace:
            raise ValueError(f"Workspace marker points elsewhere: {marker}")
        actual_fitter_hash = sha256_file(fitter)
        if payload.get("fitter_sha256") != actual_fitter_hash:
            raise ValueError(
                "Workspace fitter changed after explicit initialization; "
                f"rebuild and reinitialize the isolated workspace: {fitter}"
            )
    return workspace


def initialize_workspace(workspace: Path) -> Path:
    workspace = assert_isolated_layout(workspace, require_marker=False)
    marker = workspace / MARKER_NAME
    if marker.exists() or marker.is_symlink():
        raise FileExistsError(marker)
    fitter = analysis_root(workspace) / "get_continuum_params.py"
    payload = {
        "schema_version": 1,
        "workspace": str(workspace),
        "initialized_utc": datetime.now(timezone.utc).isoformat(),
        "fitter_sha256": sha256_file(fitter),
        "policy": (
            "ALMA/VLA inputs must match the recorded preparation mode; "
            "SMA CASA tables and all writable parent directories must be copies"
        ),
    }
    marker.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return marker


def _copy_release_ignore(_directory: str, names: list[str]) -> set[str]:
    ignored = {".git", "__pycache__"}
    ignored.update(name for name in names if name.endswith((".pyc", ".pyo")))
    ignored.update(name for name in names if name.startswith("casa-") and name.endswith(".log"))
    return ignored.intersection(names)


def _require_outside(path: Path, protected: Path, label: str) -> None:
    try:
        path.resolve().relative_to(protected.resolve())
    except ValueError:
        return
    raise ValueError(f"{label} must be outside {protected}: {path}")


def _portable_workspace_path(value: str, field: str) -> Path:
    text = value.strip().replace("\\", "/").removeprefix("./")
    path = Path(text)
    windows_absolute = len(text) >= 3 and text[1:3] == ":/"
    if not text or path.is_absolute() or windows_absolute or ".." in path.parts:
        raise ValueError(f"{field} must be a release-relative path: {value!r}")
    return path


def prepare_workspace(
    *,
    release: Path,
    archive_project_root: Path,
    workspace: Path,
    date: str,
    flagnum: int,
    sma_discovery_relative: Path,
    archive_input_mode: str = "copy",
) -> Path:
    """Build an isolated current-code workspace without running CASA."""

    release = release.resolve()
    archive_project_root = archive_project_root.resolve()
    workspace = workspace.absolute()
    if archive_input_mode not in {"copy", "symlink"}:
        raise ValueError(
            "Archive input mode must be 'copy' (strong isolation) or 'symlink'"
        )
    archive_analysis = archive_project_root / "analysis/multi_freq_from_archive"
    selected_manifest = release / "analysis/sma_calibration/selected_clean_iterations.csv"
    addition_policy_source = (
        release
        / "analysis/multi_freq_from_archive/continuum_product_addition_policy.csv"
    )
    boundary_policy_source = (
        release
        / "analysis/multi_freq_from_archive/continuum_boundary_rejection_policy.csv"
    )
    canonical_source = (
        release / "analysis/multi_freq_from_archive" / HISTORICAL_CANONICAL_RELATIVE
    )
    for required in (
        release,
        archive_project_root,
        archive_analysis,
        selected_manifest,
        addition_policy_source,
        boundary_policy_source,
        canonical_source,
    ):
        if not required.exists():
            raise FileNotFoundError(required)
    canonical_rows = read_sorted_table(canonical_source)
    if len(canonical_rows) != EXPECTED_CANONICAL_ROWS:
        raise ValueError(
            f"Immutable canonical baseline must contain {EXPECTED_CANONICAL_ROWS} rows"
        )
    # The current final45/84 sidecars are NOT pre-application lineage. Load the
    # independent authenticated history; bind current release bytes before copy.
    bootstrap_payloads, bootstrap_evidence = derive_source18_bootstrap(
        analysis_root(release)
    )
    addition_rows = load_reviewed_product_additions(
        addition_policy_source,
        analysis_directory=archive_analysis,
        canonical_rows=canonical_rows,
    )
    boundary_rows = load_boundary_rejection_policy(
        boundary_policy_source,
        analysis_root=archive_analysis,
    )
    boundary_products = boundary_rejection_products(boundary_rows)
    canonical_products = {
        portable_release_product(row[7]) for row in canonical_rows
    }
    if boundary_products & canonical_products:
        raise ValueError(
            "Boundary-rejection products must be absent from the immutable "
            "publication baseline"
        )
    if boundary_products & addition_products(addition_rows):
        raise ValueError(
            "Boundary-rejection and newly admitted product inventories overlap"
        )
    additions_by_telescope = {
        telescope: sum(row["telescope"] == telescope for row in addition_rows)
        for telescope in EXPECTED_ARCHIVE_INPUTS
    }
    if workspace.exists() or workspace.is_symlink():
        raise FileExistsError(workspace)
    _require_outside(workspace, release, "Workspace")
    _require_outside(workspace, archive_project_root, "Workspace")
    if release == Path("/") or archive_project_root == Path("/") or workspace == Path("/"):
        raise ValueError("Release, archive, and workspace roots must be specific directories")

    frozen_manifest = analysis_root(release) / "frozen_fitter_inputs.csv"
    frozen_selection = validate_and_select_frozen_inputs(
        archive_project_root, manifest=frozen_manifest
    )
    frozen_rows = load_frozen_fitter_inputs(frozen_manifest)
    frozen_expected = {r['relative_path']: r['expected_sha256'] for r in frozen_rows}

    shutil.copytree(release, workspace, ignore=_copy_release_ignore)
    root = analysis_root(workspace)
    for generated_name in (
        "ALMA",
        "VLA",
        "images",
        "chosen_frequency",
        "bad",
        "separations",
        "guessregions",
        "regrids",
        "workflow_runs",
    ):
        generated = root / generated_name
        if generated.is_symlink():
            raise ValueError(f"Copied release contains generated symlink: {generated}")
        if generated.is_dir():
            shutil.rmtree(generated)
        elif generated.exists():
            generated.unlink()

    # A release can legitimately preserve evidence from an earlier run that
    # used the same date/flag label.  Those files are outputs, not inputs, and
    # must not make a newly prepared no-clobber workspace look dirty.
    prefix = f"{date}_{flagnum}"
    for generated in expected_output_paths(root, prefix, flagnum):
        if generated.is_dir() and not generated.is_symlink():
            shutil.rmtree(generated)
        elif generated.exists() or generated.is_symlink():
            generated.unlink()

    source18_bootstrap = initialize_source18_bootstrap(
        root, bootstrap_payloads, bootstrap_evidence
    )

    archive_counts: dict[str, int] = {}
    archive_input_sha256: dict[str, str] = {}
    for telescope, baseline_count in EXPECTED_ARCHIVE_INPUTS.items():
        expected_count = baseline_count + additions_by_telescope[telescope]
        source_root = archive_analysis / telescope
        destination_root = root / telescope
        destination_root.mkdir()
        count = 0
        for image_path in sorted(frozen_selection.discovery_paths):
            if not image_path.is_relative_to(source_root):
                continue
            if not image_path.is_file() or image_path.is_symlink():
                raise ValueError(f"Frozen archive input is not a real file: {image_path}")
            for ancestor in image_path.parents:
                if ancestor.is_symlink():
                    raise ValueError(f"Frozen archive input has a symlink parent: {ancestor}")
                if ancestor == source_root:
                    break
            destination = destination_root / image_path.relative_to(source_root)
            destination.parent.mkdir(parents=True, exist_ok=True)
            relative_product = destination.relative_to(root).as_posix()
            source_sha256 = sha256_file(image_path)
            if source_sha256 != frozen_expected[relative_product]:
                raise ValueError(f"Frozen archive input changed before copy: {image_path}")
            if archive_input_mode == "copy":
                shutil.copy2(image_path, destination)
                if sha256_file(destination) != source_sha256:
                    raise RuntimeError(f"Copied archive input failed SHA-256 verification: {destination}")
            else:
                destination.symlink_to(image_path)
            archive_input_sha256[relative_product] = source_sha256
            count += 1
        if count != expected_count:
            raise ValueError(
                f"Expected {expected_count} {telescope} archive inputs, found {count}"
            )
        archive_counts[telescope] = count

    sma_relative = _portable_workspace_path(
        sma_discovery_relative.as_posix(), "SMA discovery path"
    )
    archive_sma = archive_project_root / sma_relative
    workspace_sma = workspace / sma_relative
    if not archive_sma.is_dir() or archive_sma.is_symlink():
        raise ValueError(f"SMA discovery path must be a real archive directory: {archive_sma}")
    frozen_sma = [p for p in frozen_selection.discovery_paths
                  if not p.is_relative_to(archive_analysis)]
    if any(p.parent != archive_sma for p in frozen_sma):
        raise ValueError("SMA discovery path differs from the frozen historical selection")
    workspace_sma.mkdir(parents=True)
    for image in frozen_sma:
        shutil.copytree(image, workspace_sma / image.name)
    sma_images = sorted(
        path
        for path in workspace_sma.iterdir()
        if path.is_dir() and path.name.endswith("_test.image")
    )
    if len(sma_images) != EXPECTED_SMA_DISCOVERY_INPUTS:
        raise ValueError(
            f"Expected {EXPECTED_SMA_DISCOVERY_INPUTS} SMA discovery images, "
            f"found {len(sma_images)}"
        )
    sma_discovery_sha256 = {
        path.relative_to(workspace).as_posix(): tree_sha256(path)
        for path in sma_images
    }
    expected_sma = {r['relative_path']: r['expected_sha256'] for r in frozen_rows
                    if r['input_role'] == 'sma_discovery_input'}
    if sma_discovery_sha256 != expected_sma:
        raise ValueError("Copied SMA discovery hashes differ from frozen membership")

    workspace_manifest = workspace / "analysis/sma_calibration/selected_clean_iterations.csv"
    with selected_manifest.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or ())
        selected_rows = [
            dict(row)
            for row in reader
            if row.get("date_code", "").strip() == date
            and row.get("flagnum", "").strip() == str(flagnum)
        ]
    if not selected_rows:
        raise ValueError(f"No selected-clean rows for date={date}, flag={flagnum}")
    selected_requested = {_portable_workspace_path(r[f], f).as_posix()
                          for r in selected_rows for f in ('image_path','residual_path')}
    selected_frozen = {p.relative_to(archive_project_root).as_posix()
                       for p in frozen_selection.selected_clean_paths}
    if selected_requested != selected_frozen:
        raise ValueError("Selected-clean manifest differs from frozen support membership")
    for row in selected_rows:
        for field in ("image_path", "residual_path"):
            relative = _portable_workspace_path(row[field], field)
            source = archive_project_root / relative
            destination = workspace / relative
            if not source.is_dir() or source.is_symlink():
                raise ValueError(f"Selected SMA {field} is not a real archive directory: {source}")
            if destination.exists() or destination.is_symlink():
                raise FileExistsError(destination)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(source, destination)
    with workspace_manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected_rows)
    selected_clean_sha256: dict[str, str] = {}
    for row in selected_rows:
        for field in ("image_path", "residual_path"):
            relative = _portable_workspace_path(row[field], field)
            selected_clean_sha256[relative.as_posix()] = tree_sha256(
                workspace / relative
            )
    expected_selected = {r['relative_path']: r['expected_sha256'] for r in frozen_rows
                         if r['input_role'] == 'sma_selected_clean_support'}
    if selected_clean_sha256 != expected_selected:
        raise ValueError("Copied selected-clean hashes differ from frozen membership")

    marker = initialize_workspace(workspace)
    preparation = {
        "schema_version": PREPARED_INPUT_SCHEMA_VERSION,
        "release": str(release),
        "archive_project_root": str(archive_project_root),
        "workspace": str(workspace),
        "date": date,
        "flagnum": flagnum,
        "source18_bootstrap": source18_bootstrap,
        "frozen_fitter_manifest_sha256": frozen_selection.manifest_sha256,
        "archive_input_counts": archive_counts,
        "archive_input_sha256": archive_input_sha256,
        "baseline_archive_input_counts": EXPECTED_ARCHIVE_INPUTS,
        "addition_policy_sha256": sha256_file(addition_policy_source),
        "addition_policy_rows": len(addition_rows),
        "addition_products": sorted(addition_products(addition_rows)),
        "publication_addition_products": sorted(
            publication_addition_products(addition_rows)
        ),
        "boundary_rejection_policy_sha256": sha256_file(boundary_policy_source),
        "boundary_rejection_policy_rows": len(boundary_rows),
        "boundary_rejection_products": sorted(boundary_products),
        "total_discovery_inputs": sum(archive_counts.values()) + len(sma_images),
        "sma_discovery_relative": sma_relative.as_posix(),
        "sma_discovery_images": len(sma_images),
        "sma_discovery_sha256": sma_discovery_sha256,
        "selected_clean_rows": len(selected_rows),
        "selected_clean_manifest_sha256": sha256_file(workspace_manifest),
        "selected_clean_sha256": selected_clean_sha256,
        "materialization": {
            "release": "copy_without_git_or_generated_fit_outputs",
            "alma_vla": (
                "independent_file_copies_in_real_source_directories"
                if archive_input_mode == "copy"
                else "individual_file_symlinks_in_real_source_directories"
            ),
            "archive_input_mode": archive_input_mode,
            "sma": "isolated_casa_directory_copies",
        },
        "marker": marker.name,
    }
    (workspace / ".continuum-fit-inputs.json").write_text(
        json.dumps(preparation, indent=2, sort_keys=True) + "\n"
    )
    return workspace_sma


def _inside(workspace: Path, path: Path, label: str) -> Path:
    resolved = path.resolve()
    try:
        resolved.relative_to(workspace.resolve())
    except ValueError as exc:
        raise ValueError(f"{label} must be inside the isolated workspace: {path}") from exc
    return resolved


def assert_no_directory_symlinks(root: Path, label: str) -> None:
    if root.is_symlink():
        raise ValueError(f"{label} is a directory symlink: {root}")
    for item in root.rglob("*"):
        if item.is_symlink() and item.resolve().is_dir():
            raise ValueError(f"{label} contains a directory symlink: {item}")


def validate_selected_clean_manifest(workspace: Path, manifest: Path) -> int:
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle)
        missing = {"image_path", "residual_path"}.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(
                f"Selected-clean manifest is missing fields: {sorted(missing)}"
            )
        rows = [dict(row) for row in reader]
    if not rows:
        raise ValueError(f"Selected-clean manifest is empty: {manifest}")
    for row in rows:
        for field in ("image_path", "residual_path"):
            relative = _portable_workspace_path(row[field], field)
            path = workspace / relative
            resolved = _inside(workspace, path, f"selected-clean {field}")
            if path.is_symlink() or not resolved.is_dir():
                raise ValueError(
                    f"Selected-clean {field} must be an isolated CASA directory: {path}"
                )
            assert_no_directory_symlinks(resolved, f"selected-clean {field}")
    return len(rows)


def validate_inputs(
    workspace: Path,
    smapaths: Sequence[Path],
    selected_manifest: Path | None,
) -> tuple[list[Path], Path | None]:
    root = analysis_root(workspace)
    for telescope in ("ALMA", "VLA"):
        directory = root / telescope
        if not directory.is_dir() or directory.is_symlink():
            raise ValueError(
                f"{telescope} must be a real writable directory in the workspace: {directory}"
            )
        assert_no_directory_symlinks(directory, f"{telescope} input tree")

    checked_sma = []
    for supplied in smapaths:
        candidate = supplied if supplied.is_absolute() else workspace / supplied
        if candidate.is_symlink():
            raise ValueError(f"SMA input is a directory symlink: {candidate}")
        resolved = _inside(workspace, candidate, "SMA path")
        if not resolved.is_dir():
            raise FileNotFoundError(resolved)
        assert_no_directory_symlinks(resolved, "SMA input")
        checked_sma.append(resolved)

    checked_manifest = None
    if selected_manifest is not None:
        candidate = (
            selected_manifest
            if selected_manifest.is_absolute()
            else workspace / selected_manifest
        )
        checked_manifest = _inside(workspace, candidate, "selected-clean manifest")
        if not checked_manifest.is_file():
            raise FileNotFoundError(checked_manifest)
        validate_selected_clean_manifest(workspace, checked_manifest)
    return checked_sma, checked_manifest


def _json_sha256(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def validate_prepared_archive_inputs(workspace: Path) -> dict[str, object]:
    """Verify the deterministic 455-plus-admissions input layout and metadata."""

    preparation_path = workspace / ".continuum-fit-inputs.json"
    if not preparation_path.is_file():
        raise ValueError(
            "The workspace was not made by the deterministic prepare command: "
            f"missing {preparation_path}"
        )
    preparation = json.loads(preparation_path.read_text())
    if preparation.get("schema_version") != PREPARED_INPUT_SCHEMA_VERSION:
        raise ValueError(
            "Prepared input manifest has a stale schema version; rebuild the "
            "isolated workspace"
        )
    archive_project_root = Path(preparation.get("archive_project_root", "")).resolve()
    if not archive_project_root.is_dir():
        raise FileNotFoundError(archive_project_root)
    mode = preparation.get("materialization", {}).get("archive_input_mode")
    if mode not in {"copy", "symlink"}:
        raise ValueError(f"Unrecognized prepared archive input mode: {mode!r}")

    root = analysis_root(workspace)
    addition_policy = root / "continuum_product_addition_policy.csv"
    source18_bootstrap = validate_source18_bootstrap_snapshot(
        root, preparation.get("source18_bootstrap")
    )
    frozen_manifest = root / "frozen_fitter_inputs.csv"
    frozen_rows = load_frozen_fitter_inputs(frozen_manifest)
    frozen_manifest_sha256 = sha256_file(frozen_manifest)
    if preparation.get("frozen_fitter_manifest_sha256") != frozen_manifest_sha256:
        raise ValueError("Prepared frozen fitter manifest changed")
    boundary_policy = root / "continuum_boundary_rejection_policy.csv"
    canonical_table = root / HISTORICAL_CANONICAL_RELATIVE
    canonical_rows = read_sorted_table(canonical_table)
    if len(canonical_rows) != EXPECTED_CANONICAL_ROWS:
        raise ValueError(
            f"Immutable canonical baseline must contain {EXPECTED_CANONICAL_ROWS} rows"
        )
    addition_rows = load_reviewed_product_additions(
        addition_policy,
        analysis_directory=root,
        canonical_rows=canonical_rows,
    )
    if preparation.get("addition_policy_sha256") != sha256_file(addition_policy):
        raise ValueError("Prepared addition-policy hash changed")
    expected_addition_products = sorted(addition_products(addition_rows))
    expected_publication_additions = sorted(
        publication_addition_products(addition_rows)
    )
    if (
        preparation.get("addition_policy_rows") != len(addition_rows)
        or preparation.get("addition_products") != expected_addition_products
        or preparation.get("publication_addition_products")
        != expected_publication_additions
    ):
        raise ValueError("Prepared addition-policy inventory changed")
    boundary_rows = load_boundary_rejection_policy(
        boundary_policy,
        analysis_root=root,
    )
    expected_boundary_products = sorted(
        boundary_rejection_products(boundary_rows)
    )
    canonical_products = {
        portable_release_product(row[7]) for row in canonical_rows
    }
    if set(expected_boundary_products) & canonical_products:
        raise ValueError(
            "Boundary-rejection products overlap the immutable publication baseline"
        )
    if set(expected_boundary_products) & set(expected_addition_products):
        raise ValueError(
            "Boundary-rejection and addition-policy products overlap"
        )
    if (
        preparation.get("boundary_rejection_policy_sha256")
        != sha256_file(boundary_policy)
        or preparation.get("boundary_rejection_policy_rows") != len(boundary_rows)
        or preparation.get("boundary_rejection_products")
        != expected_boundary_products
    ):
        raise ValueError("Prepared boundary-rejection policy inventory changed")

    recorded_sha256 = preparation.get("archive_input_sha256")
    if not isinstance(recorded_sha256, dict):
        raise ValueError("Prepared input manifest lacks the archive SHA-256 mapping")
    records: list[dict[str, object]] = []
    observed_sha256: dict[str, str] = {}
    counts: dict[str, int] = {}
    for telescope, baseline_count in EXPECTED_ARCHIVE_INPUTS.items():
        expected_count = baseline_count + sum(
            row["telescope"] == telescope for row in addition_rows
        )
        telescope_root = root / telescope
        archive_telescope = (
            archive_project_root / "analysis/multi_freq_from_archive" / telescope
        ).resolve()
        products = sorted(
            (
                path
                for path in telescope_root.rglob("*")
                if path.name.lower().endswith((".fits", ".imfits"))
            ),
            key=lambda path: path.relative_to(analysis_root(workspace)).as_posix(),
        )
        if len(products) != expected_count:
            raise ValueError(
                f"Prepared {telescope} tree must contain exactly {expected_count} "
                f"FITS inputs, found {len(products)}"
            )
        for product in products:
            if mode == "copy" and (product.is_symlink() or not product.is_file()):
                raise ValueError(f"Prepared archive input is not an independent file: {product}")
            if mode == "symlink":
                if not product.is_symlink() or not product.is_file():
                    raise ValueError(f"Prepared archive input is not a valid file link: {product}")
                target = product.resolve(strict=True)
                try:
                    target.relative_to(archive_telescope)
                except ValueError as exc:
                    raise ValueError(
                        f"Prepared archive link points outside the declared archive: {product}"
                    ) from exc
            else:
                target = product
            stat = target.stat()
            relative_product = product.relative_to(root).as_posix()
            digest = sha256_file(target)
            observed_sha256[relative_product] = digest
            if recorded_sha256.get(relative_product) != digest:
                raise ValueError(
                    f"Prepared archive input failed recorded SHA-256 check: {product}"
                )
            records.append(
                {
                    "product": relative_product,
                    "target": str(target),
                    "size": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                }
            )
        counts[telescope] = len(products)

    if set(observed_sha256) != set(recorded_sha256):
        raise ValueError(
            "Prepared archive SHA-256 inventory does not match the exact base-plus-"
            "addition products"
        )

    sma_relative = _portable_workspace_path(
        str(preparation.get("sma_discovery_relative", "")),
        "prepared SMA discovery path",
    )
    sma_root = workspace / sma_relative
    if sma_root.is_symlink() or not sma_root.is_dir():
        raise ValueError(f"Prepared SMA discovery tree is not an isolated copy: {sma_root}")
    assert_no_directory_symlinks(sma_root, "prepared SMA discovery tree")
    sma_count = sum(
        1
        for path in sma_root.iterdir()
        if path.is_dir() and path.name.endswith("_test.image")
    )
    if sma_count != EXPECTED_SMA_DISCOVERY_INPUTS:
        raise ValueError(
            f"Prepared SMA tree must contain {EXPECTED_SMA_DISCOVERY_INPUTS} images, "
            f"found {sma_count}"
        )
    recorded_sma_sha256 = preparation.get("sma_discovery_sha256")
    if not isinstance(recorded_sma_sha256, dict):
        raise ValueError("Prepared input manifest lacks SMA discovery tree hashes")
    sma_images = sorted(
        (
            path
            for path in sma_root.iterdir()
            if path.is_dir() and path.name.endswith("_test.image")
        ),
        key=lambda path: path.relative_to(workspace).as_posix(),
    )
    observed_sma_sha256 = {
        path.relative_to(workspace).as_posix(): tree_sha256(path)
        for path in sma_images
    }
    if observed_sma_sha256 != recorded_sma_sha256:
        raise ValueError("Prepared SMA discovery tree hashes changed")

    selected_manifest = (
        workspace / "analysis/sma_calibration/selected_clean_iterations.csv"
    )
    if selected_manifest.is_symlink() or not selected_manifest.is_file():
        raise FileNotFoundError(selected_manifest)
    if preparation.get("selected_clean_manifest_sha256") != sha256_file(
        selected_manifest
    ):
        raise ValueError("Prepared selected-clean manifest hash changed")
    with selected_manifest.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if {"image_path", "residual_path"}.difference(reader.fieldnames or ()):
            raise ValueError("Prepared selected-clean manifest has an invalid schema")
        selected_rows = [dict(row) for row in reader]
    if len(selected_rows) != preparation.get("selected_clean_rows"):
        raise ValueError("Prepared selected-clean row count changed")
    selected_paths: list[Path] = []
    for row in selected_rows:
        for field in ("image_path", "residual_path"):
            relative = _portable_workspace_path(row[field], field)
            selected_path = workspace / relative
            if selected_path.is_symlink() or not selected_path.is_dir():
                raise ValueError(
                    f"Prepared selected-clean {field} is not an isolated directory: "
                    f"{selected_path}"
                )
            assert_no_directory_symlinks(selected_path, f"selected-clean {field}")
            selected_paths.append(selected_path)
    selected_relative_paths = [
        path.relative_to(workspace).as_posix() for path in selected_paths
    ]
    if len(selected_relative_paths) != len(set(selected_relative_paths)):
        raise ValueError("Prepared selected-clean manifest repeats a CASA tree")
    observed_selected_sha256 = {
        path.relative_to(workspace).as_posix(): tree_sha256(path)
        for path in sorted(
            selected_paths, key=lambda item: item.relative_to(workspace).as_posix()
        )
    }
    if observed_selected_sha256 != preparation.get("selected_clean_sha256"):
        raise ValueError("Prepared selected-clean CASA tree hashes changed")
    for role, observed in (
        ('archive_fitter_input', observed_sha256),
        ('sma_discovery_input', observed_sma_sha256),
        ('sma_selected_clean_support', observed_selected_sha256),
    ):
        expected = {r['relative_path']: r['expected_sha256'] for r in frozen_rows
                    if r['input_role'] == role}
        if observed != expected:
            raise ValueError(f"Prepared {role} inventory/hash differs from frozen membership")
    total = sum(counts.values()) + sma_count
    expected_total = sum(EXPECTED_ARCHIVE_INPUTS.values()) + len(addition_rows) + sma_count
    if total != expected_total:
        raise ValueError(
            f"Prepared discovery inventory must total {expected_total} inputs, found {total}"
        )
    return {
        "archive_input_mode": mode,
        "source18_bootstrap": source18_bootstrap,
        "frozen_fitter_manifest_sha256": frozen_manifest_sha256,
        "counts": counts,
        "baseline_archive_input_counts": dict(EXPECTED_ARCHIVE_INPUTS),
        "addition_policy_rows": len(addition_rows),
        "addition_policy_sha256": sha256_file(addition_policy),
        "addition_products": expected_addition_products,
        "publication_addition_products": expected_publication_additions,
        "boundary_rejection_policy_rows": len(boundary_rows),
        "boundary_rejection_policy_sha256": sha256_file(boundary_policy),
        "boundary_rejection_products": expected_boundary_products,
        "sma_discovery_inputs": sma_count,
        "sma_discovery_content_sha256": _json_sha256(observed_sma_sha256),
        "selected_clean_rows": len(selected_rows),
        "selected_clean_manifest_sha256": sha256_file(selected_manifest),
        "selected_clean_content_sha256": _json_sha256(
            observed_selected_sha256
        ),
        "total_discovery_inputs": total,
        "file_metadata_sha256": _json_sha256(records),
        "file_content_sha256": _json_sha256(observed_sha256),
        "file_metadata": records,
    }


def read_sorted_table(path: Path) -> list[list[object]]:
    rows = [
        ast.literal_eval(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if any(not isinstance(row, list) or len(row) != 8 for row in rows):
        raise ValueError(f"Unexpected continuum sorted-table schema: {path}")
    return rows


def _rows_by_product(rows: Sequence[list[object]], label: str) -> dict[str, list[object]]:
    by_product: dict[str, list[object]] = {}
    for row in rows:
        product = portable_release_product(row[7])
        if product in by_product:
            raise ValueError(f"Duplicate {label} product: {product}")
        by_product[product] = row
    return by_product


def load_reviewed_exclusion_products(exclusion_policy: Path) -> set[str]:
    """Load the reviewed product set and reject count/schema/identity drift."""

    with exclusion_policy.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if "product" not in (reader.fieldnames or ()):
            raise ValueError(f"{exclusion_policy}: missing product column")
        policy_rows = list(reader)
    policy_products = [
        portable_release_product(row.get("product", "")) for row in policy_rows
    ]
    excluded = set(policy_products)
    if not policy_rows:
        raise ValueError(f"{exclusion_policy}: reviewed exclusion policy is empty")
    if "" in excluded or len(excluded) != len(policy_rows):
        raise ValueError(
            f"{exclusion_policy}: exclusion products must be nonblank and unique"
        )
    return excluded


def load_reviewed_product_additions(
    addition_policy: Path,
    *,
    analysis_directory: Path | None = None,
    canonical_rows: Sequence[list[object]] | None = None,
) -> list[dict[str, str]]:
    """Load additions and require them to be absent from the 409-row baseline."""

    rows = load_product_addition_policy(
        addition_policy,
        analysis_root=analysis_directory,
    )
    if canonical_rows is not None:
        baseline_products = {
            portable_release_product(row[7]) for row in canonical_rows
        }
        overlap = sorted(addition_products(rows) & baseline_products)
        if overlap:
            raise ValueError(
                "Product additions must be absent from the immutable canonical "
                f"baseline: {overlap}"
            )
        baseline_date_keys = {
            (
                str(row[0]),
                str(row[1]),
                f"{float(row[2]):.2f}",
                str(row[4]).split("_", maxsplit=1)[0],
            )
            for row in canonical_rows
        }
        addition_date_keys: set[tuple[str, str, str, str]] = set()
        for row in rows:
            key = (
                fitter_output_source(row),
                row["telescope"],
                row["expected_frequency_ghz"],
                row["expected_date_label"].split("_", maxsplit=1)[0],
            )
            if key in baseline_date_keys or key in addition_date_keys:
                raise ValueError(
                    "Addition shares source/telescope/rounded-frequency/UT-date "
                    "with another product; the current fitter cannot prove that "
                    f"historical date-label suffixes remain stable: {row['product']}"
                )
            addition_date_keys.add(key)
    return rows


def compare_current_corrected_table(
    generated: Path,
    canonical: Path,
    exclusion_policy: Path,
    addition_policy: Path | None = None,
) -> dict[str, object]:
    """Compare a rebuild with strict common lineage plus reviewed additions."""

    generated_rows = read_sorted_table(generated)
    canonical_rows = read_sorted_table(canonical)
    if len(canonical_rows) != EXPECTED_CANONICAL_ROWS:
        raise ValueError(
            f"Canonical baseline must contain {EXPECTED_CANONICAL_ROWS} rows, "
            f"found {len(canonical_rows)}"
        )
    excluded = load_reviewed_exclusion_products(exclusion_policy)
    additions = (
        []
        if addition_policy is None
        else load_reviewed_product_additions(
            addition_policy,
            canonical_rows=canonical_rows,
        )
    )
    admitted = addition_products(additions)
    publication_additions = publication_addition_products(additions)

    canonical_by_product = _rows_by_product(canonical_rows, "canonical")
    generated_by_product = _rows_by_product(generated_rows, "generated")
    if not excluded.issubset(canonical_by_product):
        raise ValueError(
            "Exclusion policy contains products absent from the canonical baseline: "
            f"{sorted(excluded - canonical_by_product)}"
        )
    expected_common_by_product = {
        product: row
        for product, row in canonical_by_product.items()
        if product not in excluded
    }
    expected_common_order = [
        portable_release_product(row[7])
        for row in canonical_rows
        if portable_release_product(row[7]) not in excluded
    ]
    generated_product_order = [
        portable_release_product(row[7]) for row in generated_rows
    ]
    generated_common_order = [
        product for product in generated_product_order if product not in admitted
    ]
    generated_addition_order = [
        product for product in generated_product_order if product in admitted
    ]
    expected_common_rows = len(canonical_rows) - len(excluded)
    expected_corrected_rows = expected_common_rows + len(publication_additions)
    if len(expected_common_by_product) != expected_common_rows:
        raise ValueError(
            f"Corrected canonical expectation must contain {expected_common_rows} "
            f"common rows, found {len(expected_common_by_product)}"
        )
    expected_products = set(expected_common_by_product) | publication_additions
    missing = sorted(expected_products - set(generated_by_product))
    extra = sorted(set(generated_by_product) - expected_products)
    mismatches = []
    for product in sorted(set(expected_common_by_product) & set(generated_by_product)):
        expected_fields = [
            str(value) for value in expected_common_by_product[product][:7]
        ]
        generated_fields = [str(value) for value in generated_by_product[product][:7]]
        if expected_fields != generated_fields:
            mismatches.append(
                {
                    "product": product,
                    "canonical_fields": expected_fields,
                    "generated_fields": generated_fields,
                }
            )
    common_product_order_matches = generated_common_order == expected_common_order
    expected_addition_order = publication_addition_products_in_fitter_order(
        additions
    )
    addition_product_order_matches = (
        generated_addition_order == expected_addition_order
    )
    product_order_matches = (
        common_product_order_matches and addition_product_order_matches
    )
    first_order_mismatch = None
    if not common_product_order_matches:
        for index in range(max(len(expected_common_order), len(generated_common_order))):
            expected_product = (
                expected_common_order[index]
                if index < len(expected_common_order)
                else None
            )
            generated_product = (
                generated_common_order[index]
                if index < len(generated_common_order)
                else None
            )
            if expected_product != generated_product:
                first_order_mismatch = {
                    "corrected_row_number": index + 1,
                    "canonical_product": expected_product,
                    "generated_product": generated_product,
                }
                break
    identities_match = not missing and not extra and product_order_matches
    fields_match = not mismatches
    matched = identities_match and fields_match and not additions
    return {
        "schema_version": 3,
        "status": "exact_match" if matched else "review_required",
        "comparison_basis": (
            f"immutable {len(canonical_rows)}-row canonical baseline after removing "
            f"{len(excluded)} reviewed exclusions plus "
            f"{len(publication_additions)} policy-required publication additions"
        ),
        "canonical": str(canonical),
        "canonical_sha256": sha256_file(canonical),
        "generated": str(generated),
        "generated_sha256": sha256_file(generated),
        "exclusion_policy": str(exclusion_policy),
        "exclusion_policy_sha256": sha256_file(exclusion_policy),
        "addition_policy": str(addition_policy) if addition_policy else None,
        "addition_policy_sha256": (
            sha256_file(addition_policy) if addition_policy else None
        ),
        "canonical_rows": len(canonical_rows),
        "excluded_products": sorted(excluded),
        "addition_policy_rows": len(additions),
        "publication_addition_products": sorted(publication_additions),
        "input_only_addition_products": sorted(admitted - publication_additions),
        "expected_common_rows": expected_common_rows,
        "expected_corrected_rows": expected_corrected_rows,
        "generated_rows": len(generated_rows),
        "missing_products": missing,
        "extra_products": extra,
        "product_order_matches": product_order_matches,
        "common_product_order_matches": common_product_order_matches,
        "addition_product_order_matches": addition_product_order_matches,
        "first_order_mismatch": first_order_mismatch,
        "field_mismatch_count": len(mismatches),
        "field_mismatches": mismatches,
    }


def _require_exact_mapping_keys(
    value: object,
    expected: Sequence[str],
    label: str,
) -> dict[str, object]:
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object")
    actual = set(value)
    required = set(expected)
    if actual != required:
        raise ValueError(
            f"{label} fields differ from the exact schema; "
            f"missing={sorted(required - actual)}, extra={sorted(actual - required)}"
        )
    return value


def load_reviewed_delta_rows(path: Path) -> list[dict[str, str]]:
    """Load a phase-specific, product-unique current-to-generated ledger."""

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != REVIEW_DELTA_FIELDS:
            raise ValueError(
                f"{path}: reviewed-delta fields must exactly equal "
                f"{list(REVIEW_DELTA_FIELDS)}"
            )
        rows = [dict(row) for row in reader]
    products: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        if any(value is None or value == "" for value in row.values()):
            raise ValueError(
                f"{path}:{row_number}: reviewed-delta fields must be nonblank"
            )
        raw_product = row["product"]
        product = portable_release_product(raw_product)
        if product != raw_product:
            raise ValueError(
                f"{path}:{row_number}: product must already be release-relative: "
                f"{raw_product!r}"
            )
        if product in products:
            raise ValueError(
                f"{path}:{row_number}: duplicate reviewed-delta product: {product}"
            )
        products.add(product)
    return rows


def _expected_review_delta_row(
    corrected_row_number: int,
    product: str,
    current_fields: Sequence[str],
    generated_fields: Sequence[str],
    changed_fields: Sequence[str],
) -> dict[str, str]:
    return {
        "corrected_row_number": str(corrected_row_number),
        "source": current_fields[0],
        "telescope": current_fields[1],
        "frequency_ghz": current_fields[2],
        "date_label": current_fields[4],
        "product": product,
        "current_flux_field": current_fields[3],
        "generated_flux_field": generated_fields[3],
        "current_snr": current_fields[6],
        "generated_snr": generated_fields[6],
        "changed_fields": "+".join(changed_fields),
        "review_reason": REVIEW_SCOPE,
    }


def compare_reviewed_current_to_generated(
    generated: Path,
    canonical: Path,
    exclusion_policy: Path,
    reviewed_delta_csv: Path,
    addition_policy: Path | None = None,
) -> dict[str, object]:
    """Require an exact reviewed ledger for every and only reported-field delta."""

    generated_rows = read_sorted_table(generated)
    canonical_rows = read_sorted_table(canonical)
    if len(canonical_rows) != EXPECTED_CANONICAL_ROWS:
        raise ValueError(
            f"Adjudication requires {EXPECTED_CANONICAL_ROWS} canonical rows, "
            f"found {len(canonical_rows)}"
        )
    excluded = load_reviewed_exclusion_products(exclusion_policy)
    if len(excluded) != EXPECTED_REVIEWED_EXCLUSIONS:
        raise ValueError(
            f"Adjudication requires exactly {EXPECTED_REVIEWED_EXCLUSIONS} reviewed "
            f"exclusions, found {len(excluded)}"
        )

    additions = (
        []
        if addition_policy is None
        else load_reviewed_product_additions(
            addition_policy,
            canonical_rows=canonical_rows,
        )
    )
    admitted = addition_products(additions)
    publication_additions = publication_addition_products(additions)
    canonical_by_product = _rows_by_product(canonical_rows, "canonical")
    generated_by_product = _rows_by_product(generated_rows, "generated")
    if not excluded.issubset(canonical_by_product):
        raise ValueError(
            "Exclusion policy contains products absent from the canonical baseline: "
            f"{sorted(excluded - canonical_by_product)}"
        )
    expected_rows = [
        row
        for row in canonical_rows
        if portable_release_product(row[7]) not in excluded
    ]
    if len(expected_rows) != EXPECTED_CORRECTED_ROWS:
        raise ValueError(
            f"Adjudication requires {EXPECTED_CORRECTED_ROWS} canonical products after "
            f"exclusions, found {len(expected_rows)}"
        )
    expected_generated_rows = EXPECTED_CORRECTED_ROWS + len(publication_additions)
    if len(generated_rows) != expected_generated_rows:
        raise ValueError(
            f"Adjudication requires {expected_generated_rows} generated rows, "
            f"found {len(generated_rows)}"
        )

    expected_products = [portable_release_product(row[7]) for row in expected_rows]
    all_generated_products = [
        portable_release_product(row[7]) for row in generated_rows
    ]
    generated_products = [
        product for product in all_generated_products if product not in admitted
    ]
    observed_additions = [
        product for product in all_generated_products if product in admitted
    ]
    missing = sorted(
        (set(expected_products) | publication_additions)
        - set(all_generated_products)
    )
    extra = sorted(
        set(all_generated_products)
        - (set(expected_products) | publication_additions)
    )
    if missing or extra:
        raise ValueError(
            "Generated product identities differ from the exact common-plus-addition "
            f"expectation; missing={missing}, extra={extra}"
        )
    if generated_products != expected_products:
        first = next(
            index
            for index, (current, new) in enumerate(
                zip(expected_products, generated_products), start=1
            )
            if current != new
        )
        raise ValueError(
            "Generated product order differs from the canonical-minus-exclusions "
            f"order at corrected row {first}: expected={expected_products[first - 1]!r}, "
            f"generated={generated_products[first - 1]!r}"
        )
    expected_addition_order = publication_addition_products_in_fitter_order(
        additions
    )
    if observed_additions != expected_addition_order:
        raise ValueError("Generated addition-product order differs from policy order")
    if set(generated_by_product) != set(expected_products) | publication_additions:
        raise ValueError("Generated identity inventory is internally inconsistent")
    generated_common_rows = [
        row
        for row in generated_rows
        if portable_release_product(row[7]) not in admitted
    ]

    allowed_changes = {"flux_field", "snr"}
    expected_delta_rows: list[dict[str, str]] = []
    flux_change_count = 0
    snr_change_count = 0
    classification_change_count = 0
    for corrected_row_number, (current_row, generated_row, product) in enumerate(
        zip(expected_rows, generated_common_rows, expected_products), start=1
    ):
        current_fields = [str(value) for value in current_row[:7]]
        generated_fields = [str(value) for value in generated_row[:7]]
        changed_fields = [
            field
            for field, current, new in zip(
                SORTED_REPORTED_FIELD_NAMES, current_fields, generated_fields
            )
            if current != new
        ]
        disallowed = set(changed_fields) - allowed_changes
        if disallowed:
            raise ValueError(
                f"Unadjudicable non-flux/SNR fields changed for {product}: "
                f"{sorted(disallowed)}"
            )
        if not changed_fields:
            continue
        if "flux_field" in changed_fields:
            flux_change_count += 1
        if "snr" in changed_fields:
            snr_change_count += 1
        current_is_limit = current_fields[3].endswith("*") or current_fields[6] == "NA"
        generated_is_limit = (
            generated_fields[3].endswith("*") or generated_fields[6] == "NA"
        )
        if current_is_limit != generated_is_limit:
            classification_change_count += 1
        expected_delta_rows.append(
            _expected_review_delta_row(
                corrected_row_number,
                product,
                current_fields,
                generated_fields,
                changed_fields,
            )
        )
    reviewed_rows = load_reviewed_delta_rows(reviewed_delta_csv)
    expected_delta_products = [row["product"] for row in expected_delta_rows]
    reviewed_products = [row["product"] for row in reviewed_rows]
    missing_reviews = sorted(set(expected_delta_products) - set(reviewed_products))
    stale_reviews = sorted(set(reviewed_products) - set(expected_delta_products))
    if missing_reviews or stale_reviews:
        raise ValueError(
            "Reviewed-delta inventory is not the exact actual-delta inventory; "
            f"missing={missing_reviews}, extra={stale_reviews}"
        )
    if reviewed_products != expected_delta_products:
        raise ValueError(
            "Reviewed-delta rows are not in canonical corrected-row order"
        )
    for expected, reviewed in zip(expected_delta_rows, reviewed_rows):
        if reviewed != expected:
            differing = [
                field
                for field in REVIEW_DELTA_FIELDS
                if reviewed[field] != expected[field]
            ]
            raise ValueError(
                f"Reviewed delta is not exact for {expected['product']}; "
                f"fields={differing}"
            )

    return {
        "schema_version": 1,
        "status": "exact_reviewed_match",
        "review_scope": REVIEW_SCOPE,
        "canonical_rows": len(canonical_rows),
        "reviewed_exclusion_rows": len(excluded),
        "validated_product_identities": len(expected_products),
        "validated_addition_identities": len(publication_additions),
        "product_order_matches": True,
        "changed_product_count": len(expected_delta_rows),
        "flux_field_change_count": flux_change_count,
        "snr_change_count": snr_change_count,
        "classification_change_count": classification_change_count,
        "reviewed_deltas": expected_delta_rows,
    }


def derive_reviewed_delta_rows(
    generated: Path,
    canonical: Path,
    exclusion_policy: Path,
    addition_policy: Path,
) -> list[dict[str, str]]:
    """Derive the exact common-product flux/SNR delta ledger in baseline order."""

    comparison = compare_current_corrected_table(
        generated,
        canonical,
        exclusion_policy,
        addition_policy,
    )
    if (
        comparison["missing_products"]
        or comparison["extra_products"]
        or not comparison["product_order_matches"]
    ):
        raise ValueError("Cannot render review deltas from a drifting product inventory")
    canonical_rows = read_sorted_table(canonical)
    excluded = load_reviewed_exclusion_products(exclusion_policy)
    additions = load_reviewed_product_additions(
        addition_policy,
        canonical_rows=canonical_rows,
    )
    admitted = addition_products(additions)
    common_rows = [
        row
        for row in canonical_rows
        if portable_release_product(row[7]) not in excluded
    ]
    generated_common = [
        row
        for row in read_sorted_table(generated)
        if portable_release_product(row[7]) not in admitted
    ]
    rendered: list[dict[str, str]] = []
    for row_number, (current, new) in enumerate(
        zip(common_rows, generated_common), start=1
    ):
        current_fields = [str(value) for value in current[:7]]
        generated_fields = [str(value) for value in new[:7]]
        changed = [
            field
            for field, old, fresh in zip(
                SORTED_REPORTED_FIELD_NAMES,
                current_fields,
                generated_fields,
            )
            if old != fresh
        ]
        disallowed = set(changed) - {"flux_field", "snr"}
        if disallowed:
            raise ValueError(
                f"Unadjudicable non-flux/SNR fields changed for {new[7]}: "
                f"{sorted(disallowed)}"
            )
        if changed:
            rendered.append(
                _expected_review_delta_row(
                    row_number,
                    portable_release_product(new[7]),
                    current_fields,
                    generated_fields,
                    changed,
                )
            )
    return rendered


def derive_reviewed_addition_rows(
    generated: Path,
    canonical: Path,
    addition_policy: Path,
) -> list[dict[str, str]]:
    """Derive one exact post-run provenance row for every admitted input."""

    canonical_rows = read_sorted_table(canonical)
    additions = load_reviewed_product_additions(
        addition_policy,
        canonical_rows=canonical_rows,
    )
    generated_rows = read_sorted_table(generated)
    generated_by_product = _rows_by_product(generated_rows, "generated")
    row_numbers = {
        portable_release_product(row[7]): index
        for index, row in enumerate(generated_rows, start=1)
    }
    rendered: list[dict[str, str]] = []
    for policy in additions:
        product = policy["product"]
        generated_row = generated_by_product.get(product)
        expected = policy["expected_run_disposition"]
        if expected == "publication_row_required":
            if generated_row is None:
                raise ValueError(f"Required addition publication row is missing: {product}")
            expected_source = fitter_output_source(policy)
            if (
                str(generated_row[0]) != expected_source
                or str(generated_row[1]) != policy["telescope"]
                or str(generated_row[2]) != policy["expected_frequency_ghz"]
                or str(generated_row[4]) != policy["expected_date_label"]
            ):
                raise ValueError(f"Addition publication identity is stale: {product}")
            observed = "publication_row_present"
            values = [str(value) for value in generated_row[:7]]
            generated_row_number = str(row_numbers[product])
        else:
            if generated_row is not None:
                raise ValueError(
                    f"Input-only addition unexpectedly produced a publication row: {product}"
                )
            observed = "no_publication_row"
            values = ["NA"] * 7
            generated_row_number = "NA"
        rendered.append(
            {
                "admission_id": policy["admission_id"],
                "source": policy["source"],
                "telescope": policy["telescope"],
                "product": product,
                "product_sha256": policy["product_sha256"],
                "expected_run_disposition": expected,
                "observed_run_disposition": observed,
                "generated_row_number": generated_row_number,
                "frequency_ghz": values[2],
                "flux_field": values[3],
                "date_label": values[4],
                "beam_arcsec": values[5],
                "snr": values[6],
                "review_reason": ADDITION_REVIEW_SCOPE,
            }
        )
    return rendered


def load_reviewed_addition_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != REVIEW_ADDITION_FIELDS:
            raise ValueError(f"{path}: reviewed-addition fields differ from exact schema")
        rows = [dict(row) for row in reader]
    if any(value is None or value == "" for row in rows for value in row.values()):
        raise ValueError(f"{path}: reviewed-addition fields must be nonblank")
    return rows


def compare_reviewed_product_additions(
    generated: Path,
    canonical: Path,
    addition_policy: Path,
    reviewed_addition_csv: Path,
) -> dict[str, object]:
    expected = derive_reviewed_addition_rows(generated, canonical, addition_policy)
    reviewed = load_reviewed_addition_rows(reviewed_addition_csv)
    if reviewed != expected:
        raise ValueError("Reviewed-addition ledger is not the exact derived inventory")
    published = sum(
        row["observed_run_disposition"] == "publication_row_present"
        for row in expected
    )
    return {
        "schema_version": 1,
        "status": "exact_reviewed_additions_match",
        "review_scope": ADDITION_REVIEW_SCOPE,
        "reviewed_input_additions": len(expected),
        "publication_additions": published,
        "input_only_additions": len(expected) - published,
        "reviewed_additions": expected,
    }


def load_review_freeze(
    path: Path,
    expected_sha256: str,
    run_id: str,
) -> dict[str, object]:
    """Load an explicitly hash-frozen review bundle with an exact schema."""

    if not SHA256_PATTERN.fullmatch(expected_sha256):
        raise ValueError("--review-freeze-sha256 must be 64 lowercase hex characters")
    actual_sha256 = sha256_file(path)
    if actual_sha256 != expected_sha256:
        raise ValueError(
            f"Review-freeze SHA-256 mismatch: expected {expected_sha256}, "
            f"found {actual_sha256}"
        )
    payload = _require_exact_mapping_keys(
        json.loads(path.read_text(encoding="utf-8")),
        REVIEW_FREEZE_FIELDS,
        "review freeze",
    )
    if payload["schema_version"] != 3:
        raise ValueError("Review-freeze schema_version must equal 3")
    if payload["review_scope"] != REVIEW_SCOPE:
        raise ValueError(f"Review-freeze scope must equal {REVIEW_SCOPE!r}")
    if payload["review_status"] != REVIEW_STATUS:
        raise ValueError(f"Review-freeze status must equal {REVIEW_STATUS!r}")
    if payload["run_id"] != run_id:
        raise ValueError("Review-freeze run_id does not match --run-id")

    counts = _require_exact_mapping_keys(
        payload["expected_counts"],
        REVIEW_FREEZE_COUNT_FIELDS,
        "review freeze expected_counts",
    )
    for field, value in counts.items():
        if type(value) is not int or value < 0:
            raise ValueError(f"Review-freeze count {field} must be a nonnegative integer")
    exact_counts = {
        "canonical_rows": EXPECTED_CANONICAL_ROWS,
        "reviewed_exclusion_rows": EXPECTED_REVIEWED_EXCLUSIONS,
        "ngc4945_photometry_rows": EXPECTED_NGC4945_ROWS,
        "ngc4945_component_rows": 7,
        "boundary_rejection_policy_rows": 3,
        "boundary_rejected_noise_audit_rows": 3,
    }
    for field, expected in exact_counts.items():
        if counts[field] != expected:
            raise ValueError(
                f"Review-freeze count {field} must equal {expected}, "
                f"found {counts[field]}"
            )
    for field in (
        "preapplication_correction_rows",
        "postrun_correction_rows",
        "classification_identity_rows",
    ):
        if counts[field] < 1:
            raise ValueError(f"Review-freeze count {field} must be positive")
    if counts["publication_addition_rows"] > counts["addition_policy_rows"]:
        raise ValueError("Publication additions exceed admitted input additions")
    if counts["reviewed_addition_rows"] != counts["addition_policy_rows"]:
        raise ValueError("Reviewed-addition rows must equal addition-policy rows")
    expected_generated = EXPECTED_CORRECTED_ROWS + counts["publication_addition_rows"]
    if (
        counts["generated_rows"] != expected_generated
        or counts["generic_publication_noise_audit_rows"] != expected_generated
        or counts["generic_noise_audit_rows"]
        != expected_generated + counts["boundary_rejected_noise_audit_rows"]
    ):
        raise ValueError(
            "Generated/publication-noise counts must equal 406 plus publication "
            "additions, and raw noise must add the exact boundary rejections"
        )

    artifacts = _require_exact_mapping_keys(
        payload["artifacts"],
        REVIEW_FREEZE_ARTIFACT_NAMES,
        "review freeze artifacts",
    )
    for name, record in artifacts.items():
        artifact = _require_exact_mapping_keys(
            record, ("path", "sha256"), f"review freeze artifact {name}"
        )
        if not isinstance(artifact["path"], str):
            raise ValueError(f"Review-freeze artifact {name} path must be a string")
        relative = _portable_workspace_path(
            artifact["path"], f"review freeze artifact {name} path"
        )
        if relative.as_posix() != artifact["path"]:
            raise ValueError(
                f"Review-freeze artifact {name} path is not normalized"
            )
        digest = artifact["sha256"]
        if not isinstance(digest, str) or not SHA256_PATTERN.fullmatch(digest):
            raise ValueError(
                f"Review-freeze artifact {name} SHA-256 is not lowercase hex"
            )
    return payload


def _csv_data_row_count(path: Path) -> int:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"CSV artifact is empty: {path}") from exc
        if not header or any(not field for field in header):
            raise ValueError(f"CSV artifact has an invalid header: {path}")
        return sum(1 for row in reader if row)


def _csv_value_count(path: Path, field: str, value: str) -> int:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if field not in (reader.fieldnames or ()):
            raise ValueError(f"CSV artifact lacks {field!r}: {path}")
        rows = list(reader)
    if any(row.get(field) is None for row in rows):
        raise ValueError(f"CSV artifact has an incomplete {field!r} value: {path}")
    return sum(row[field] == value for row in rows)


def _read_exact_csv_rows(
    path: Path,
    expected_fields: Sequence[str],
    *,
    label: str,
) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        actual_fields = tuple(reader.fieldnames or ())
        if actual_fields != tuple(expected_fields):
            raise ValueError(
                f"{label} fields differ from the exact schema: {actual_fields!r}"
            )
        rows = [dict(row) for row in reader]
    if not rows:
        raise ValueError(f"{label} cannot be empty")
    return rows


def _normalized_correction_rows(
    path: Path,
    *,
    allow_legacy: bool,
    label: str,
) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = tuple(reader.fieldnames or ())
        raw_rows = [dict(row) for row in reader]
    if not raw_rows:
        raise ValueError(f"{label} cannot be empty")
    if fields == LEGACY_CORRECTION_FIELDS and allow_legacy:
        rows = [
            {
                "schema_version": CORRECTION_SCHEMA_VERSION,
                **{
                    field: row[field]
                    for field in LEGACY_CORRECTION_FIELDS[:-1]
                },
                "validated_adopted_noise_snr": row[
                    "validated_background_only_snr"
                ],
            }
            for row in raw_rows
        ]
    elif fields == CORRECTION_FIELDS:
        rows = raw_rows
    else:
        raise ValueError(
            f"{label} fields differ from the reviewed legacy/v2 schemas: {fields!r}"
        )
    if any(row["schema_version"] != CORRECTION_SCHEMA_VERSION for row in rows):
        raise ValueError(f"{label} contains a stale schema version")
    products = [portable_release_product(row["product"]) for row in rows]
    if "" in products or len(products) != len(set(products)):
        raise ValueError(f"{label} products must be nonblank and unique")
    for row, product in zip(rows, products):
        row["product"] = product
    return rows


def validate_correction_manifest_transition(
    preapplication_path: Path,
    postrun_path: Path,
    source18_photometry_path: Path,
) -> dict[str, object]:
    """Prove the source-18 stage changed only its reviewed manifest records."""

    before = _normalized_correction_rows(
        preapplication_path,
        allow_legacy=True,
        label="pre-application correction manifest",
    )
    after = _normalized_correction_rows(
        postrun_path,
        allow_legacy=False,
        label="post-run correction manifest",
    )
    expected_preapplication_rows = 66
    expected_new_ngc4945_rows = EXPECTED_NGC4945_ROWS - 27
    expected_postrun_rows = expected_preapplication_rows + expected_new_ngc4945_rows
    if len(before) != expected_preapplication_rows or len(after) != expected_postrun_rows:
        raise ValueError(
            "Correction-manifest transition requires exactly "
            f"{expected_preapplication_rows} pre-application and "
            f"{expected_postrun_rows} post-run rows"
        )
    before_by_product = {row["product"]: row for row in before}
    after_by_product = {row["product"]: row for row in after}
    if not set(before_by_product) < set(after_by_product):
        raise ValueError(
            "Pre-application correction products are not a proper subset of the expanded manifest"
        )

    photometry = _read_exact_csv_rows(
        source18_photometry_path,
        PHOTOMETRY_FIELDS,
        label="source-18 photometry",
    )
    if len(photometry) != EXPECTED_NGC4945_ROWS:
        raise ValueError(
            f"Source-18 photometry must contain exactly {EXPECTED_NGC4945_ROWS} rows"
        )
    photometry_by_product = {
        portable_release_product(row["product"]): row for row in photometry
    }
    if len(photometry_by_product) != len(photometry):
        raise ValueError("Source-18 photometry contains duplicate products")
    previous_ngc4945 = {
        row["product"] for row in before if row["source"] == "NGC4945"
    }
    postrun_ngc4945 = {
        row["product"] for row in after if row["source"] == "NGC4945"
    }
    if not previous_ngc4945 < set(photometry_by_product):
        raise ValueError(
            "Previously released NGC4945 correction products are not a proper "
            "subset of the expanded source-18 photometry"
        )
    if postrun_ngc4945 != set(photometry_by_product):
        raise ValueError(
            "Post-run NGC4945 correction products differ from source-18 photometry"
        )
    new_products = set(after_by_product) - set(before_by_product)
    if new_products != set(photometry_by_product) - previous_ngc4945:
        raise ValueError("Only newly admitted NGC4945 products may extend the correction manifest")

    for product, old in before_by_product.items():
        if old["source"] != "NGC4945" and after_by_product.get(product) != old:
            raise ValueError(f"Non-NGC4945 correction-manifest row changed: {product}")

    stable_fields = (
        "schema_version",
        "source",
        "telescope",
        "frequency_ghz",
        "date_label",
        "product",
        "old_flux_field",
        "old_snr",
    )
    updated_products: list[str] = []
    for new in after:
        if new["source"] != "NGC4945":
            continue
        product = new["product"]
        old = before_by_product.get(product)
        measurement = photometry_by_product[product]
        if old is not None and any(old[field] != new[field] for field in stable_fields):
            raise ValueError(
                f"Stable NGC4945 correction-manifest fields changed: {product}"
            )
        for field in ("old_flux_field", "old_snr"):
            if measurement[field] != new[field] or (
                old is not None and measurement[field] != old[field]
            ):
                raise ValueError(
                    f"NGC4945 source-18 photometry {field} disagrees with the "
                    f"pre/post correction-manifest lineage for {product}"
                )
        expected_values = {
            "new_flux_field": measurement["new_flux_field"],
            "new_snr": measurement["new_snr"],
            "validated_unrounded_flux_mjy": measurement["unrounded_flux_mjy"],
            "validated_adopted_noise_snr": measurement["new_snr"],
        }
        for field, expected in expected_values.items():
            if new[field] != expected:
                raise ValueError(
                    f"NGC4945 correction-manifest {field} disagrees with "
                    f"source-18 photometry for {product}"
                )
        reasons = [] if old is None else [
            reason
            for reason in old["correction_reasons"].split("+")
            if reason not in {"single_background_rms", "single_local_image_scatter"}
        ]
        for required in (
            "fixed_nuclear_component_source18",
            "single_local_scaled_mad",
            "source18_10bmaj_local_noise_aperture",
        ):
            if required not in reasons:
                reasons.append(required)
        if new["correction_reasons"] != "+".join(reasons):
            raise ValueError(
                f"NGC4945 correction reasons are not the exact reviewed transition "
                f"for {product}"
            )
        updated_products.append(product)

    if len(updated_products) != EXPECTED_NGC4945_ROWS:
        raise ValueError(
            "Correction-manifest transition did not review "
            f"{EXPECTED_NGC4945_ROWS} NGC4945 rows"
        )
    return {
        "schema_version": 1,
        "status": "exact_source18_only_transition",
        "preapplication_rows": len(before),
        "postrun_rows": len(after),
        "unchanged_non_ngc4945_rows": len(before) - len(previous_ngc4945),
        "added_ngc4945_rows": len(new_products),
        "reviewed_ngc4945_rows": len(updated_products),
        "reviewed_ngc4945_products": updated_products,
    }


def expected_output_paths(root: Path, prefix: str, flagnum: int) -> list[Path]:
    outputs = [
        root / "images/chosenfreq",
        root / f"chosen_frequency/candidates/{prefix}_getchooser.csv",
        root / "chosen_frequency/chosen_frequency_anchors.csv",
    ]
    for _, pipetype in STAGES:
        outputs.extend(
            (
                root / f"fitsumfiles/witherrors/{prefix}_fitsummary_{pipetype}_witherrors.txt",
                root / f"fitsumfiles/noerrors/{prefix}_fitsummary_{pipetype}_noerrors.txt",
                root / f"fitsumfiles/sorted/{prefix}_fitsummary_{pipetype}_sorted.txt",
                root / f"fitsumfiles/fluxerrors/{prefix}_fluxerr_{pipetype}.txt",
            )
        )
    outputs.append(
        root
        / f"fitsumfiles/noise_audits/{prefix}_final_publication_noise_audit.csv"
    )
    outputs.extend(root.glob(f"images/{flagnum}/*"))
    for telescope in ("ALMA", "VLA"):
        outputs.extend(root.glob(f"{telescope}/*/fitfiles"))
    return list(outputs)


def assert_clean_outputs(root: Path, prefix: str, flagnum: int) -> None:
    conflicts = [path for path in expected_output_paths(root, prefix, flagnum) if path.exists()]
    if conflicts:
        formatted = "\n".join(f"  {path}" for path in conflicts[:20])
        raise FileExistsError(
            "Refusing to overwrite existing continuum workflow outputs:\n" + formatted
        )


def stage_commands(
    *,
    python_bin: Path,
    fitter: Path,
    date: str,
    flagnum: int,
    calnum: int,
    smapaths: Sequence[Path],
    selected_manifest: Path | None,
    use_selected_clean: bool,
) -> list[tuple[str, list[str]]]:
    commands = []
    for mode, _ in STAGES:
        command = [
            str(python_bin),
            str(fitter),
            "--flagnum",
            str(flagnum),
            "--calnum",
            str(calnum),
            "--mode",
            mode,
            "--datestr",
            date,
            "--chooser",
            "0",
            "--choose_science",
            "NGC4258",
            "--printimages",
            "0",
            "--writesummary",
            "1",
            "--nofixmultval",
            "10",
            "--printbadfix",
            "0",
            "--telechooser",
            "0",
            "--telechoice",
            "SMA",
            "--moveniter",
            "0",
            "--niternum",
            "0",
            "--smapath",
            *(str(path) for path in smapaths),
        ]
        if mode == "fixtochosenfreq" and use_selected_clean:
            if selected_manifest is None:
                raise ValueError(
                    "--use-selected-clean requires --selected-clean-manifest"
                )
            command.extend(
                ("--selected_clean_manifest", str(selected_manifest), "--use_selected_clean_manifest")
            )
        commands.append((mode, command))
    root = fitter.parent
    final_table = (
        root
        / f"fitsumfiles/sorted/{date}_{flagnum}_fitsummary_fixedtochosen_sorted.txt"
    )
    correction_script = root / "rerun_ngc4945_nuclear.py"
    correction_common = [
        str(python_bin),
        str(correction_script),
        "--source-table",
        str(final_table),
        "--output",
        str(root / "ngc4945_nuclear_photometry.csv"),
        "--component-output",
        str(root / "ngc4945_component_audit.csv"),
        "--correction-manifest",
        str(root / "validated_continuum_corrections.csv"),
    ]
    commands.append(
        (
            "ngc4945_fixed_nucleus",
            correction_common
            + ["--archive-root", str(root), "--apply"],
        )
    )
    commands.append(("ngc4945_fixed_nucleus_check", correction_common + ["--check"]))
    return commands


def _artifact_snapshot(root: Path, stage: str, flagnum: int) -> dict[str, object]:
    image_stage = {
        "prelimchecker": "prelim",
        "fit_sizenoloc": "noloc",
        "getchoosercoords": "chosenfreq",
        "fixtochosenfreq": "fixedtochosen",
    }.get(stage)
    candidates = {
        "fitsumfiles": root / "fitsumfiles",
        "bad": root / "bad",
        "separations": root / "separations",
        "chosen_frequency": root / "chosen_frequency",
        "chosen_models": root / "images/chosenfreq",
        "guessregions": root / "guessregions",
        "regrids": root / "regrids",
        "ALMA_tree_with_generated_fitfiles": root / "ALMA",
        "VLA_tree_with_generated_fitfiles": root / "VLA",
        "SMA_generated_fitfiles": root / "SMA",
    }
    if image_stage is not None:
        candidates["stage_images"] = root / f"images/{flagnum}/{image_stage}"
    snapshot: dict[str, object] = {"stage": stage, "artifacts": {}}
    artifacts = snapshot["artifacts"]
    assert isinstance(artifacts, dict)
    for name, path in candidates.items():
        if path.is_dir():
            artifacts[name] = {
                "path": path.relative_to(root).as_posix(),
                "tree_sha256": tree_sha256(path),
                "files": sum(1 for item in path.rglob("*") if item.is_file()),
            }
    return snapshot


def _safe_run_id(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", value):
        raise ValueError(f"Unsafe run id: {value!r}")
    return value


def build_environment(analysis_utils_dir: Path) -> dict[str, str]:
    if not (analysis_utils_dir / "analysisUtils.py").is_file():
        raise FileNotFoundError(analysis_utils_dir / "analysisUtils.py")
    environment = dict(os.environ)
    environment.update(
        {
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONHASHSEED": "0",
            "LC_ALL": "C",
            "PYTHONPATH": str(analysis_utils_dir.resolve()),
        }
    )
    environment.pop("LD_LIBRARY_PATH", None)
    return environment


def preflight_python(python_bin: Path, environment: dict[str, str]) -> None:
    # Merely importing CASA creates a timestamped log in the current working
    # directory.  Confine that side effect to a disposable directory so
    # ``plan`` remains read-only with respect to both release and workspace.
    # Import the full third-party dependency surface used by the fitter and
    # source-18 correction, not merely CASA itself; this prevents an hour-long
    # workflow from discovering a missing package at its first science stage.
    with tempfile.TemporaryDirectory(prefix="continuum-casa-preflight-") as temporary:
        subprocess.run(
            [
                str(python_bin),
                "-c",
                (
                    "import analysisUtils, astropy, astroquery.alma, "
                    "astroquery.nvas, casatasks, casatools, matplotlib, numpy, "
                    "reproject, requests; print('continuum fitter imports OK')"
                ),
            ],
            cwd=temporary,
            env=environment,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )


def workflow_input_hashes(root: Path) -> dict[str, str]:
    """Hash every small code/policy input that controls the four-pass handoff."""

    hashes = {}
    for name in WORKFLOW_INPUT_NAMES:
        path = root / name
        if not path.is_file():
            raise FileNotFoundError(path)
        hashes[name] = sha256_file(path)
    return hashes


def _workspace_file(workspace: Path, supplied: Path | None, default: Path) -> Path:
    candidate = default if supplied is None else supplied
    candidate = candidate if candidate.is_absolute() else workspace / candidate
    resolved = _inside(workspace, candidate, "comparison input")
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    return resolved


def run_workflow(args: argparse.Namespace, *, dry_run: bool) -> int:
    workspace = assert_isolated_layout(args.workspace, require_marker=True)
    root = analysis_root(workspace)
    input_hashes = workflow_input_hashes(root)
    prepared_inputs_before = validate_prepared_archive_inputs(workspace)
    validate_preapplication_lineage(root, prepared_inputs_before["source18_bootstrap"])
    smapaths, selected_manifest = validate_inputs(
        workspace, args.smapath, args.selected_clean_manifest
    )
    expected_selected_manifest = (
        workspace / "analysis/sma_calibration/selected_clean_iterations.csv"
    ).resolve()
    if args.use_selected_clean and selected_manifest != expected_selected_manifest:
        raise ValueError(
            "Selected-clean execution must use the prepared canonical manifest: "
            f"{expected_selected_manifest}"
        )
    canonical_reference = _workspace_file(
        workspace,
        args.canonical_reference,
        root / HISTORICAL_CANONICAL_RELATIVE,
    )
    exclusion_policy = _workspace_file(
        workspace,
        args.exclusion_policy,
        root / "continuum_product_exclusion_policy.csv",
    )
    reviewed_exclusions = load_reviewed_exclusion_products(exclusion_policy)
    addition_policy = _workspace_file(
        workspace,
        None,
        root / "continuum_product_addition_policy.csv",
    )
    boundary_policy = _workspace_file(
        workspace,
        None,
        root / "continuum_boundary_rejection_policy.csv",
    )
    canonical_rows = read_sorted_table(canonical_reference)
    reviewed_additions = load_reviewed_product_additions(
        addition_policy,
        analysis_directory=root,
        canonical_rows=canonical_rows,
    )
    expected_publication_additions = publication_addition_products(
        reviewed_additions
    )
    reviewed_boundary_rejections = load_boundary_rejection_policy(
        boundary_policy,
        analysis_root=root,
    )
    prefix = f"{args.date}_{args.flagnum}"
    assert_clean_outputs(root, prefix, args.flagnum)
    commands = stage_commands(
        python_bin=args.python_bin,
        fitter=root / "get_continuum_params.py",
        date=args.date,
        flagnum=args.flagnum,
        calnum=args.calnum,
        smapaths=smapaths,
        selected_manifest=selected_manifest,
        use_selected_clean=args.use_selected_clean,
    )

    analysis_utils_source = args.analysis_utils_dir.resolve()
    analysis_utils_preflight_inventory = analysis_utils_bundle_inventory(
        analysis_utils_source
    )
    environment = build_environment(args.analysis_utils_dir)
    preflight_python(args.python_bin, environment)

    if dry_run:
        print(
            "validated isolated input inventory: "
            f"{prepared_inputs_before['total_discovery_inputs']} discovery inputs "
            f"({prepared_inputs_before['archive_input_mode']} mode)"
        )
        print(
            f"corrected comparison target: {EXPECTED_CANONICAL_ROWS} canonical rows "
            f"minus {len(reviewed_exclusions)} reviewed NGC 4258 exclusions = "
            f"{EXPECTED_CANONICAL_ROWS - len(reviewed_exclusions)} common rows, plus "
            f"{len(expected_publication_additions)} publication additions from "
            f"{len(reviewed_additions)} admitted inputs; "
            f"{len(reviewed_boundary_rejections)} reviewed boundary diagnostics "
            "remain outside the publication table"
        )
        for mode, command in commands:
            print(f"{mode}: {shlex.join(command)}")
        return 0

    run_id = _safe_run_id(args.run_id or prefix)
    evidence = root / "workflow_runs" / run_id
    if evidence.exists() or evidence.is_symlink():
        raise FileExistsError(evidence)
    evidence.mkdir(parents=True)
    cumulative_corrections = root / "validated_continuum_corrections.csv"
    if cumulative_corrections.is_symlink() or not cumulative_corrections.is_file():
        raise FileNotFoundError(cumulative_corrections)
    preapplication_corrections = (
        evidence / "preapplication_validated_continuum_corrections.csv"
    )
    shutil.copy2(cumulative_corrections, preapplication_corrections)
    analysis_utils_execution_dir = evidence / "analysis_scripts"
    analysis_utils_snapshot = evidence / "analysis_utils_snapshot.json"
    analysis_utils_manifest = freeze_analysis_utils_bundle(
        analysis_utils_source,
        analysis_utils_execution_dir,
        analysis_utils_snapshot,
    )
    analysis_utils_inventory = _analysis_utils_inventory_from_manifest(
        analysis_utils_manifest
    )
    if analysis_utils_inventory != analysis_utils_preflight_inventory:
        raise RuntimeError(
            "analysisUtils dependency bundle changed during workflow preflight"
        )
    # Execute every CASA stage against the complete frozen dependency bundle.
    # This closes the gap between the bytes recorded in the plan and every
    # local module that the child interpreter can import through analysisUtils.
    environment = build_environment(analysis_utils_execution_dir)
    preflight_python(args.python_bin, environment)
    plan = {
        "workspace": str(workspace),
        "run_id": run_id,
        "date": args.date,
        "flagnum": args.flagnum,
        "calnum": args.calnum,
        "python_bin": str(args.python_bin),
        "fitter_sha256": sha256_file(root / "get_continuum_params.py"),
        "workflow_input_sha256": input_hashes,
        "selected_clean_manifest_sha256": (
            sha256_file(selected_manifest) if selected_manifest is not None else None
        ),
        "selected_clean_manifest": (
            str(selected_manifest) if selected_manifest is not None else None
        ),
        "analysis_utils_sha256": analysis_utils_inventory[
            "file_inventory_sha256"
        ],
        "analysis_utils_snapshot": str(analysis_utils_snapshot),
        "analysis_utils_execution_dir": str(analysis_utils_execution_dir),
        "prepared_input_inventory": {
            key: value
            for key, value in prepared_inputs_before.items()
            if key != "file_metadata"
        },
        "canonical_reference": str(canonical_reference),
        "canonical_reference_sha256": sha256_file(canonical_reference),
        "exclusion_policy": str(exclusion_policy),
        "exclusion_policy_sha256": sha256_file(exclusion_policy),
        "addition_policy": str(addition_policy),
        "addition_policy_sha256": sha256_file(addition_policy),
        "addition_policy_rows": len(reviewed_additions),
        "publication_addition_rows": len(expected_publication_additions),
        "boundary_rejection_policy": str(boundary_policy),
        "boundary_rejection_policy_sha256": sha256_file(boundary_policy),
        "boundary_rejection_policy_rows": len(reviewed_boundary_rejections),
        "reviewed_exclusion_count": len(reviewed_exclusions),
        "preapplication_cumulative_corrections": str(preapplication_corrections),
        "preapplication_cumulative_corrections_sha256": sha256_file(
            preapplication_corrections
        ),
        "commands": [{"stage": mode, "argv": command} for mode, command in commands],
    }
    (evidence / "plan.json").write_text(json.dumps(plan, indent=2) + "\n")

    for mode, command in commands:
        if workflow_input_hashes(root) != input_hashes:
            raise RuntimeError("Workflow code or policy changed before a CASA stage")
        if analysis_utils_bundle_inventory(
            analysis_utils_execution_dir
        ) != analysis_utils_inventory:
            raise RuntimeError(
                "Frozen analysisUtils dependency bundle changed before a CASA stage"
            )
        log_path = evidence / f"{mode}.log"
        with log_path.open("w") as log:
            completed = subprocess.run(
                command,
                cwd=root,
                env=environment,
                stdin=subprocess.DEVNULL,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
            )
        if completed.returncode:
            raise RuntimeError(
                f"Continuum stage {mode} failed with status {completed.returncode}; "
                f"see {log_path}"
            )
        if analysis_utils_bundle_inventory(
            analysis_utils_execution_dir
        ) != analysis_utils_inventory:
            raise RuntimeError(
                "Frozen analysisUtils dependency bundle changed during a CASA stage"
            )
        if workflow_input_hashes(root) != input_hashes:
            raise RuntimeError("Workflow code or policy changed during a CASA stage")
        snapshot = _artifact_snapshot(root, mode, args.flagnum)
        (evidence / f"{mode}.artifacts.json").write_text(
            json.dumps(snapshot, indent=2, sort_keys=True) + "\n"
        )
        for casa_log in sorted(root.glob("casa-*.log")):
            shutil.copy2(casa_log, evidence / f"{mode}.{casa_log.name}")
        if mode == "getchoosercoords" and not args.skip_reference_check:
            check_command = [
                sys.executable,
                str(root / "chosen_frequency_anchors.py"),
                "check",
                "--analysis-root",
                str(root),
                "--manifest",
                str(root / "chosen_frequency/chosen_frequency_anchors.csv"),
                "--reference",
                str(root / "chosen_frequency_anchor_reference.csv"),
                "--addition-policy",
                str(addition_policy),
                "--component-policy",
                str(root / "chosen_frequency_component_policy.csv"),
            ]
            subprocess.run(check_command, check=True, cwd=root)

    prepared_inputs_after = validate_prepared_archive_inputs(workspace)
    archive_integrity = {
        "status": (
            "unchanged"
            if prepared_inputs_before == prepared_inputs_after
            else "changed"
        ),
        "before_file_metadata_sha256": prepared_inputs_before[
            "file_metadata_sha256"
        ],
        "after_file_metadata_sha256": prepared_inputs_after["file_metadata_sha256"],
        "archive_input_mode": prepared_inputs_after["archive_input_mode"],
        "total_discovery_inputs": prepared_inputs_after["total_discovery_inputs"],
        "addition_policy_rows": prepared_inputs_after["addition_policy_rows"],
        "publication_addition_rows": len(
            prepared_inputs_after["publication_addition_products"]
        ),
        "boundary_rejection_policy_rows": prepared_inputs_after[
            "boundary_rejection_policy_rows"
        ],
    }
    (evidence / "prepared_input_integrity.json").write_text(
        json.dumps(archive_integrity, indent=2, sort_keys=True) + "\n"
    )
    if prepared_inputs_before != prepared_inputs_after:
        raise RuntimeError(
            "Prepared continuum input inventory changed during the run; "
            f"see {evidence / 'prepared_input_integrity.json'}"
        )

    generated_table = (
        root
        / f"fitsumfiles/sorted/{prefix}_fitsummary_fixedtochosen_sorted.txt"
    )
    final_noise_audit = (
        root
        / f"fitsumfiles/noise_audits/{prefix}_final_publication_noise_audit.csv"
    )
    source18_photometry = root / "ngc4945_nuclear_photometry.csv"
    source18_components = root / "ngc4945_component_audit.csv"
    correction_transition = validate_correction_manifest_transition(
        preapplication_corrections,
        cumulative_corrections,
        source18_photometry,
    )
    correction_transition_path = (
        evidence / "correction_manifest_transition_validation.json"
    )
    correction_transition_path.write_text(
        json.dumps(correction_transition, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    inventory_validation = {
        "schema_version": 1,
        "status": "exact_validated_union",
        **validate_final_noise_inventory(
            generated_table,
            final_noise_audit,
            source18_photometry,
            addition_policy,
            boundary_policy,
            root,
        ),
        **validate_ngc4945_component_artifacts(
            source18_photometry,
            source18_components,
        ),
    }
    inventory_validation_path = evidence / "final_noise_inventory_validation.json"
    inventory_validation_path.write_text(
        json.dumps(inventory_validation, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    classification_identity_policy = (
        root / "continuum_classification_identity_policy.csv"
    )
    classification_identity_validation = validate_classification_identity_policy(
        workspace=workspace,
        canonical_table=canonical_reference,
        generated_table=generated_table,
        exclusion_policy=exclusion_policy,
        generic_noise_audit=final_noise_audit,
        identity_policy=classification_identity_policy,
        addition_policy=addition_policy,
        boundary_rejection_policy=boundary_policy,
    )
    classification_identity_validation_path = (
        evidence / "classification_identity_validation.json"
    )
    classification_identity_validation_path.write_text(
        json.dumps(classification_identity_validation, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
    comparison = compare_current_corrected_table(
        generated_table,
        canonical_reference,
        exclusion_policy,
        addition_policy,
    )
    comparison_path = evidence / "current_corrected_table_comparison.json"
    comparison_path.write_text(
        json.dumps(comparison, indent=2, sort_keys=True) + "\n"
    )
    if comparison["status"] != "exact_match":
        raise RuntimeError(
            "Corrected continuum rebuild reached the intentional review_required "
            "stop because it differs from the current canonical "
            f"{comparison['expected_corrected_rows']}-row expectation; "
            f"review {comparison_path}"
        )
    print(
        "completed four generic continuum passes and the mandatory guarded "
        "NGC 4945 source-18 correction; the "
        f"{comparison['expected_corrected_rows']} corrected rows exactly match "
        f"the current canonical expectation; evidence: {evidence}"
    )
    return 0


def _atomic_write_new_json(path: Path, payload: object) -> None:
    """Atomically create one JSON evidence file without overwrite semantics."""

    if path.exists() or path.is_symlink():
        raise FileExistsError(path)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=f".{path.name}.",
            suffix=".tmp",
            dir=path.parent,
            delete=False,
        ) as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
            temporary_path = Path(handle.name)
        try:
            os.link(temporary_path, path)
        except FileExistsError as exc:
            raise FileExistsError(path) from exc
        directory_fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def _atomic_write_new_csv(
    path: Path,
    fields: Sequence[str],
    rows: Sequence[dict[str, str]],
) -> None:
    """Atomically create an exact-schema CSV without overwrite semantics."""

    if path.exists() or path.is_symlink():
        raise FileExistsError(path)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            prefix=f".{path.name}.",
            suffix=".tmp",
            dir=path.parent,
            delete=False,
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
            temporary_path = Path(handle.name)
        try:
            os.link(temporary_path, path)
        except FileExistsError as exc:
            raise FileExistsError(path) from exc
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def _read_json_object(path: Path, label: str) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object: {path}")
    return value


def _workspace_relative_file(
    workspace: Path,
    path: Path,
    label: str,
) -> tuple[Path, str]:
    resolved = _inside(workspace, path, label)
    if path.is_symlink() or not resolved.is_file():
        raise ValueError(f"{label} must be a real file inside the workspace: {path}")
    return resolved, resolved.relative_to(workspace).as_posix()


def render_review_evidence(args: argparse.Namespace) -> int:
    """Render exact ledgers and their hash freeze without CASA or table writes."""

    workspace = assert_isolated_layout(args.workspace, require_marker=True)
    root = analysis_root(workspace)
    run_id = _safe_run_id(args.run_id)
    evidence = root / "workflow_runs" / run_id
    if not evidence.is_dir() or evidence.is_symlink():
        raise FileNotFoundError(evidence)
    plan_path = evidence / "plan.json"
    plan = _read_json_object(plan_path, "workflow plan")
    if plan.get("run_id") != run_id or Path(str(plan.get("workspace", ""))).resolve() != workspace:
        raise ValueError("Workflow plan identity does not match render request")
    if workflow_input_hashes(root) != plan.get("workflow_input_sha256"):
        raise ValueError("Workflow code or policy inputs changed after the CASA run")
    prefix = f"{plan['date']}_{int(plan['flagnum'])}"
    canonical = _inside(workspace, Path(str(plan["canonical_reference"])), "canonical")
    exclusions = _inside(workspace, Path(str(plan["exclusion_policy"])), "exclusions")
    additions = _inside(workspace, Path(str(plan["addition_policy"])), "additions")
    boundaries = _inside(
        workspace,
        Path(str(plan["boundary_rejection_policy"])),
        "boundary rejections",
    )
    generated = root / f"fitsumfiles/sorted/{prefix}_fitsummary_fixedtochosen_sorted.txt"
    noise = root / f"fitsumfiles/noise_audits/{prefix}_final_publication_noise_audit.csv"
    delta_path = evidence / "reviewed_current_to_generated_deltas.csv"
    addition_path = evidence / "reviewed_new_product_additions.csv"
    freeze_path = evidence / "reviewed_local_scaled_mad_freeze.json"
    for output in (delta_path, addition_path, freeze_path):
        if output.exists() or output.is_symlink():
            raise FileExistsError(output)

    initial_path = evidence / "current_corrected_table_comparison.json"
    stored_initial = _read_json_object(initial_path, "initial comparison")
    fresh_initial = compare_current_corrected_table(
        generated, canonical, exclusions, additions
    )
    if stored_initial != fresh_initial or fresh_initial["status"] != "review_required":
        raise ValueError("Initial review-required comparison is absent or stale")

    delta_rows = derive_reviewed_delta_rows(
        generated, canonical, exclusions, additions
    )
    addition_rows = derive_reviewed_addition_rows(generated, canonical, additions)
    _atomic_write_new_csv(delta_path, REVIEW_DELTA_FIELDS, delta_rows)
    _atomic_write_new_csv(addition_path, REVIEW_ADDITION_FIELDS, addition_rows)
    compare_reviewed_current_to_generated(
        generated, canonical, exclusions, delta_path, additions
    )
    compare_reviewed_product_additions(
        generated, canonical, additions, addition_path
    )

    artifact_paths = {
        "workflow_plan": plan_path,
        "prepared_archive_inventory": workspace / ".continuum-fit-inputs.json",
        "prepared_input_integrity": evidence / "prepared_input_integrity.json",
        "initial_comparison": initial_path,
        "final_stage_snapshot": evidence / "ngc4945_fixed_nucleus_check.artifacts.json",
        "canonical_table": canonical,
        "generated_table": generated,
        "exclusion_policy": exclusions,
        "addition_policy": additions,
        "boundary_rejection_policy": boundaries,
        "generic_final_noise_audit": noise,
        "ngc4945_nuclear_photometry": root / "ngc4945_nuclear_photometry.csv",
        "ngc4945_component_audit": root / "ngc4945_component_audit.csv",
        "final_noise_inventory_validation": evidence / "final_noise_inventory_validation.json",
        "correction_manifest_transition_validation": evidence / "correction_manifest_transition_validation.json",
        "classification_identity_policy": root / "continuum_classification_identity_policy.csv",
        "classification_identity_validation": evidence / "classification_identity_validation.json",
        "analysis_utils_snapshot": evidence / "analysis_utils_snapshot.json",
        "preapplication_cumulative_corrections": evidence / "preapplication_validated_continuum_corrections.csv",
        "postrun_cumulative_corrections": root / "validated_continuum_corrections.csv",
        "reviewed_delta_csv": delta_path,
        "reviewed_addition_csv": addition_path,
    }
    if set(artifact_paths) != set(REVIEW_FREEZE_ARTIFACT_NAMES):
        raise AssertionError("Internal freeze-artifact inventory is inconsistent")
    resolved = {
        name: _workspace_relative_file(workspace, path, name)[0]
        for name, path in artifact_paths.items()
    }
    policy_rows = load_reviewed_product_additions(
        additions,
        analysis_directory=root,
        canonical_rows=read_sorted_table(canonical),
    )
    boundary_rows = load_boundary_rejection_policy(
        boundaries,
        analysis_root=root,
    )
    generic_noise_rows = _csv_data_row_count(noise)
    boundary_noise_rows = _csv_value_count(
        noise, "row_disposition", "rejected_on_boundary"
    )
    counts = {
        "canonical_rows": len(read_sorted_table(canonical)),
        "reviewed_exclusion_rows": len(load_reviewed_exclusion_products(exclusions)),
        "generated_rows": len(read_sorted_table(generated)),
        "generic_noise_audit_rows": generic_noise_rows,
        "generic_publication_noise_audit_rows": (
            generic_noise_rows - boundary_noise_rows
        ),
        "boundary_rejection_policy_rows": len(boundary_rows),
        "boundary_rejected_noise_audit_rows": boundary_noise_rows,
        "ngc4945_photometry_rows": _csv_data_row_count(resolved["ngc4945_nuclear_photometry"]),
        "ngc4945_component_rows": _csv_data_row_count(resolved["ngc4945_component_audit"]),
        "preapplication_correction_rows": _csv_data_row_count(resolved["preapplication_cumulative_corrections"]),
        "postrun_correction_rows": _csv_data_row_count(resolved["postrun_cumulative_corrections"]),
        "classification_identity_rows": _csv_data_row_count(resolved["classification_identity_policy"]),
        "reviewed_delta_rows": len(delta_rows),
        "addition_policy_rows": len(policy_rows),
        "publication_addition_rows": len(publication_addition_products(policy_rows)),
        "reviewed_addition_rows": len(addition_rows),
    }
    payload = {
        "schema_version": 3,
        "review_scope": REVIEW_SCOPE,
        "review_status": REVIEW_STATUS,
        "run_id": run_id,
        "expected_counts": counts,
        "artifacts": {
            name: {
                "path": path.relative_to(workspace).as_posix(),
                "sha256": sha256_file(path),
            }
            for name, path in resolved.items()
        },
    }
    _atomic_write_new_json(freeze_path, payload)
    load_review_freeze(freeze_path, sha256_file(freeze_path), run_id)
    print(
        f"rendered {len(delta_rows)} common-product deltas and "
        f"{len(addition_rows)} product-admission rows; freeze SHA-256 "
        f"{sha256_file(freeze_path)}"
    )
    return 0


def adjudicate_workflow(args: argparse.Namespace) -> int:
    """Validate a frozen reviewed delta ledger without invoking CASA or editing tables."""

    workspace = assert_isolated_layout(args.workspace, require_marker=True)
    root = analysis_root(workspace)
    run_id = _safe_run_id(args.run_id)
    evidence = root / "workflow_runs" / run_id
    if not evidence.is_dir() or evidence.is_symlink():
        raise FileNotFoundError(evidence)

    freeze_candidate = (
        args.review_freeze
        if args.review_freeze.is_absolute()
        else workspace / args.review_freeze
    )
    freeze_path = _inside(workspace, freeze_candidate, "review freeze")
    expected_freeze_path = evidence / "reviewed_local_scaled_mad_freeze.json"
    if freeze_path != expected_freeze_path.resolve():
        raise ValueError(
            f"Review freeze must be the run-specific evidence file: {expected_freeze_path}"
        )
    if freeze_candidate.is_symlink() or not freeze_path.is_file():
        raise ValueError(f"Review freeze must be a real file: {freeze_candidate}")
    freeze = load_review_freeze(
        freeze_path,
        args.review_freeze_sha256,
        run_id,
    )

    plan_path = evidence / "plan.json"
    plan = _read_json_object(plan_path, "workflow plan")
    if Path(str(plan.get("workspace", ""))).resolve() != workspace:
        raise ValueError("Workflow plan workspace does not match --workspace")
    if plan.get("run_id") != run_id:
        raise ValueError("Workflow plan run_id does not match --run-id")
    try:
        date = str(plan["date"])
        flagnum = int(plan["flagnum"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Workflow plan lacks a valid date/flagnum") from exc
    prefix = f"{date}_{flagnum}"

    planned_input_hashes = plan.get("workflow_input_sha256")
    if not isinstance(planned_input_hashes, dict):
        raise ValueError("Workflow plan lacks workflow_input_sha256")
    current_input_hashes = workflow_input_hashes(root)
    if current_input_hashes != planned_input_hashes:
        raise ValueError("Workflow code or policy inputs changed after the CASA run")
    executing_runner = Path(__file__).resolve()
    workspace_runner = (root / "run_continuum_four_pass.py").resolve()
    if (
        executing_runner != workspace_runner
        or sha256_file(executing_runner)
        != planned_input_hashes.get("run_continuum_four_pass.py")
    ):
        raise ValueError(
            "Adjudication must be executed by the unchanged workspace runner"
        )

    current_inventory = validate_prepared_archive_inputs(workspace)
    current_inventory_summary = {
        key: value for key, value in current_inventory.items() if key != "file_metadata"
    }
    if current_inventory_summary != plan.get("prepared_input_inventory"):
        raise ValueError("Prepared archive input inventory changed after the CASA run")
    selected_manifest = (
        workspace / "analysis/sma_calibration/selected_clean_iterations.csv"
    ).resolve()
    if (
        plan.get("selected_clean_manifest") != str(selected_manifest)
        or plan.get("selected_clean_manifest_sha256")
        != sha256_file(selected_manifest)
    ):
        raise ValueError("Selected-clean manifest differs from the workflow plan")

    integrity_path = evidence / "prepared_input_integrity.json"
    integrity = _read_json_object(integrity_path, "prepared-input integrity evidence")
    if integrity.get("status") != "unchanged":
        raise ValueError("Prepared-input integrity evidence is not unchanged")
    expected_discovery_inputs = 455 + int(
        current_inventory_summary.get("addition_policy_rows", -1)
    )
    if integrity.get("total_discovery_inputs") != expected_discovery_inputs:
        raise ValueError(
            "Prepared-input integrity evidence does not record the exact "
            "base-plus-addition inventory"
        )
    metadata_sha256 = current_inventory_summary["file_metadata_sha256"]
    if (
        integrity.get("before_file_metadata_sha256") != metadata_sha256
        or integrity.get("after_file_metadata_sha256") != metadata_sha256
    ):
        raise ValueError("Prepared-input metadata hash differs from run integrity evidence")

    canonical_path = _inside(
        workspace,
        Path(str(plan.get("canonical_reference", ""))),
        "canonical reference",
    )
    exclusion_path = _inside(
        workspace,
        Path(str(plan.get("exclusion_policy", ""))),
        "exclusion policy",
    )
    addition_path = _inside(
        workspace,
        Path(str(plan.get("addition_policy", ""))),
        "addition policy",
    )
    boundary_path = _inside(
        workspace,
        Path(str(plan.get("boundary_rejection_policy", ""))),
        "boundary-rejection policy",
    )
    addition_rows = load_reviewed_product_additions(
        addition_path,
        analysis_directory=root,
        canonical_rows=read_sorted_table(canonical_path),
    )
    publication_additions = publication_addition_products(addition_rows)
    boundary_rows = load_boundary_rejection_policy(
        boundary_path,
        analysis_root=root,
    )
    generated_path = root / (
        f"fitsumfiles/sorted/{prefix}_fitsummary_fixedtochosen_sorted.txt"
    )
    noise_audit_path = root / (
        f"fitsumfiles/noise_audits/{prefix}_final_publication_noise_audit.csv"
    )
    reviewed_delta_path = evidence / "reviewed_current_to_generated_deltas.csv"
    reviewed_addition_path = evidence / "reviewed_new_product_additions.csv"
    comparison_path = evidence / "current_corrected_table_comparison.json"
    final_snapshot_path = evidence / "ngc4945_fixed_nucleus_check.artifacts.json"
    expected_artifact_paths = {
        "workflow_plan": plan_path,
        "prepared_archive_inventory": workspace / ".continuum-fit-inputs.json",
        "prepared_input_integrity": integrity_path,
        "initial_comparison": comparison_path,
        "final_stage_snapshot": final_snapshot_path,
        "canonical_table": canonical_path,
        "generated_table": generated_path,
        "exclusion_policy": exclusion_path,
        "addition_policy": addition_path,
        "boundary_rejection_policy": boundary_path,
        "generic_final_noise_audit": noise_audit_path,
        "ngc4945_nuclear_photometry": root / "ngc4945_nuclear_photometry.csv",
        "ngc4945_component_audit": root / "ngc4945_component_audit.csv",
        "final_noise_inventory_validation": (
            evidence / "final_noise_inventory_validation.json"
        ),
        "correction_manifest_transition_validation": (
            evidence / "correction_manifest_transition_validation.json"
        ),
        "classification_identity_policy": (
            root / "continuum_classification_identity_policy.csv"
        ),
        "classification_identity_validation": (
            evidence / "classification_identity_validation.json"
        ),
        "analysis_utils_snapshot": evidence / "analysis_utils_snapshot.json",
        "preapplication_cumulative_corrections": (
            evidence / "preapplication_validated_continuum_corrections.csv"
        ),
        "postrun_cumulative_corrections": (
            root / "validated_continuum_corrections.csv"
        ),
        "reviewed_delta_csv": reviewed_delta_path,
        "reviewed_addition_csv": reviewed_addition_path,
    }
    resolved_artifacts: dict[str, Path] = {}
    frozen_artifacts = freeze["artifacts"]
    assert isinstance(frozen_artifacts, dict)
    for name in REVIEW_FREEZE_ARTIFACT_NAMES:
        expected_path, expected_relative = _workspace_relative_file(
            workspace, expected_artifact_paths[name], f"frozen artifact {name}"
        )
        record = frozen_artifacts[name]
        assert isinstance(record, dict)
        if record["path"] != expected_relative:
            raise ValueError(
                f"Frozen artifact {name} path differs from the run-derived path"
            )
        actual_digest = sha256_file(expected_path)
        if record["sha256"] != actual_digest:
            raise ValueError(
                f"Frozen artifact {name} SHA-256 differs from the current file"
            )
        resolved_artifacts[name] = expected_path

    if plan.get("canonical_reference_sha256") != sha256_file(canonical_path):
        raise ValueError("Canonical table no longer matches the workflow plan")
    if plan.get("exclusion_policy_sha256") != sha256_file(exclusion_path):
        raise ValueError("Exclusion policy no longer matches the workflow plan")
    if (
        plan.get("addition_policy_sha256") != sha256_file(addition_path)
        or plan.get("addition_policy_rows") != len(addition_rows)
        or plan.get("publication_addition_rows") != len(publication_additions)
    ):
        raise ValueError("Addition policy no longer matches the workflow plan")
    if (
        plan.get("boundary_rejection_policy_sha256")
        != sha256_file(boundary_path)
        or plan.get("boundary_rejection_policy_rows") != len(boundary_rows)
    ):
        raise ValueError(
            "Boundary-rejection policy no longer matches the workflow plan"
        )
    if (
        plan.get("preapplication_cumulative_corrections")
        != str(resolved_artifacts["preapplication_cumulative_corrections"])
        or plan.get("preapplication_cumulative_corrections_sha256")
        != sha256_file(resolved_artifacts["preapplication_cumulative_corrections"])
    ):
        raise ValueError(
            "Pre-application cumulative-correction snapshot differs from the plan"
        )
    analysis_utils_manifest = _read_json_object(
        resolved_artifacts["analysis_utils_snapshot"],
        "analysisUtils dependency-bundle manifest",
    )
    analysis_utils_inventory = _analysis_utils_inventory_from_manifest(
        analysis_utils_manifest
    )
    analysis_utils_execution_dir = _inside(
        workspace,
        Path(str(analysis_utils_manifest["execution_directory"])),
        "analysisUtils dependency-bundle execution directory",
    )
    if not analysis_utils_execution_dir.is_dir() or analysis_utils_execution_dir.is_symlink():
        raise ValueError(
            "analysisUtils dependency-bundle execution path must be a real directory"
        )
    if (
        plan.get("analysis_utils_snapshot")
        != str(resolved_artifacts["analysis_utils_snapshot"])
        or plan.get("analysis_utils_execution_dir")
        != str(analysis_utils_execution_dir)
        or analysis_utils_manifest["execution_directory"]
        != str(analysis_utils_execution_dir)
        or plan.get("analysis_utils_sha256")
        != analysis_utils_inventory["file_inventory_sha256"]
        or analysis_utils_bundle_inventory(analysis_utils_execution_dir)
        != analysis_utils_inventory
    ):
        raise ValueError(
            "analysisUtils dependency bundle differs from the workflow plan"
        )

    stored_comparison = _read_json_object(comparison_path, "initial comparison")
    if stored_comparison.get("status") != "review_required":
        raise ValueError("Adjudication requires an initial review_required stop")
    recomputed_comparison = compare_current_corrected_table(
        generated_path,
        canonical_path,
        exclusion_path,
        addition_path,
    )
    if stored_comparison != recomputed_comparison:
        raise ValueError("Initial comparison evidence differs from recomputation")
    if (
        recomputed_comparison.get("expected_corrected_rows")
        != EXPECTED_CORRECTED_ROWS + len(publication_additions)
        or recomputed_comparison.get("generated_rows")
        != EXPECTED_CORRECTED_ROWS + len(publication_additions)
        or recomputed_comparison.get("missing_products")
        or recomputed_comparison.get("extra_products")
        or not recomputed_comparison.get("product_order_matches")
    ):
        raise ValueError(
            "Initial review_required stop is not limited to ordered field deltas"
        )

    stored_snapshot = _read_json_object(final_snapshot_path, "final stage snapshot")
    recomputed_snapshot = _artifact_snapshot(
        root, "ngc4945_fixed_nucleus_check", flagnum
    )
    if stored_snapshot != recomputed_snapshot:
        raise ValueError("Final-stage artifact trees changed after their snapshot")

    recomputed_inventory_validation = {
        "schema_version": 1,
        "status": "exact_validated_union",
        **validate_final_noise_inventory(
            generated_path,
            noise_audit_path,
            resolved_artifacts["ngc4945_nuclear_photometry"],
            addition_path,
            boundary_path,
            root,
        ),
        **validate_ngc4945_component_artifacts(
            resolved_artifacts["ngc4945_nuclear_photometry"],
            resolved_artifacts["ngc4945_component_audit"],
        ),
    }
    stored_inventory_validation = _read_json_object(
        resolved_artifacts["final_noise_inventory_validation"],
        "final-noise inventory validation",
    )
    if stored_inventory_validation != recomputed_inventory_validation:
        raise ValueError(
            "Final-noise inventory validation differs from strict recomputation"
        )

    recomputed_classification_identity = validate_classification_identity_policy(
        workspace=workspace,
        canonical_table=canonical_path,
        generated_table=generated_path,
        exclusion_policy=exclusion_path,
        generic_noise_audit=noise_audit_path,
        identity_policy=resolved_artifacts["classification_identity_policy"],
        addition_policy=addition_path,
        boundary_rejection_policy=boundary_path,
    )
    stored_classification_identity = _read_json_object(
        resolved_artifacts["classification_identity_validation"],
        "classification-identity validation",
    )
    if stored_classification_identity != recomputed_classification_identity:
        raise ValueError(
            "Classification-identity validation differs from strict recomputation"
        )

    recomputed_correction_transition = validate_correction_manifest_transition(
        resolved_artifacts["preapplication_cumulative_corrections"],
        resolved_artifacts["postrun_cumulative_corrections"],
        resolved_artifacts["ngc4945_nuclear_photometry"],
    )
    stored_correction_transition = _read_json_object(
        resolved_artifacts["correction_manifest_transition_validation"],
        "correction-manifest transition validation",
    )
    if stored_correction_transition != recomputed_correction_transition:
        raise ValueError(
            "Correction-manifest transition validation differs from strict "
            "recomputation"
        )

    counts = freeze["expected_counts"]
    assert isinstance(counts, dict)
    generic_noise_rows = _csv_data_row_count(noise_audit_path)
    boundary_noise_rows = _csv_value_count(
        noise_audit_path, "row_disposition", "rejected_on_boundary"
    )
    actual_counts = {
        "canonical_rows": len(read_sorted_table(canonical_path)),
        "reviewed_exclusion_rows": len(
            load_reviewed_exclusion_products(exclusion_path)
        ),
        "generated_rows": len(read_sorted_table(generated_path)),
        "generic_noise_audit_rows": generic_noise_rows,
        "generic_publication_noise_audit_rows": (
            generic_noise_rows - boundary_noise_rows
        ),
        "boundary_rejection_policy_rows": len(boundary_rows),
        "boundary_rejected_noise_audit_rows": boundary_noise_rows,
        "ngc4945_photometry_rows": _csv_data_row_count(
            resolved_artifacts["ngc4945_nuclear_photometry"]
        ),
        "ngc4945_component_rows": _csv_data_row_count(
            resolved_artifacts["ngc4945_component_audit"]
        ),
        "preapplication_correction_rows": _csv_data_row_count(
            resolved_artifacts["preapplication_cumulative_corrections"]
        ),
        "postrun_correction_rows": _csv_data_row_count(
            resolved_artifacts["postrun_cumulative_corrections"]
        ),
        "classification_identity_rows": _csv_data_row_count(
            resolved_artifacts["classification_identity_policy"]
        ),
        "reviewed_delta_rows": _csv_data_row_count(reviewed_delta_path),
        "addition_policy_rows": len(addition_rows),
        "publication_addition_rows": len(publication_additions),
        "reviewed_addition_rows": _csv_data_row_count(reviewed_addition_path),
    }
    if counts != actual_counts:
        raise ValueError(
            f"Frozen artifact row counts differ from current files: "
            f"frozen={counts}, current={actual_counts}"
        )
    if (
        recomputed_classification_identity.get("required_review_rows")
        != actual_counts["classification_identity_rows"]
    ):
        raise ValueError(
            "Classification-identity row count differs from its strict "
            "fresh-result validation"
        )

    delta_comparison = compare_reviewed_current_to_generated(
        generated_path,
        canonical_path,
        exclusion_path,
        reviewed_delta_path,
        addition_path,
    )
    if delta_comparison["changed_product_count"] != counts["reviewed_delta_rows"]:
        raise ValueError("Reviewed-delta count differs from exact computed delta count")
    addition_comparison = compare_reviewed_product_additions(
        generated_path,
        canonical_path,
        addition_path,
        reviewed_addition_path,
    )
    if addition_comparison["reviewed_input_additions"] != counts["reviewed_addition_rows"]:
        raise ValueError("Reviewed-addition count differs from exact policy count")

    output = evidence / "reviewed_delta_adjudication.json"
    result = {
        "schema_version": 2,
        "status": "exact_reviewed_match",
        "run_id": run_id,
        "review_freeze": freeze_path.relative_to(workspace).as_posix(),
        "review_freeze_sha256": args.review_freeze_sha256,
        "adjudicator_sha256": sha256_file(executing_runner),
        "frozen_artifacts": frozen_artifacts,
        "validated_counts": actual_counts,
        "delta_comparison": delta_comparison,
        "addition_comparison": addition_comparison,
        "policy": (
            "read-only adjudication; neither canonical nor generated table was modified"
        ),
    }
    _atomic_write_new_json(output, result)
    print(
        f"adjudicated {delta_comparison['changed_product_count']} reviewed product "
        f"deltas across {EXPECTED_CORRECTED_ROWS} common ordered identities and "
        f"{addition_comparison['reviewed_input_additions']} admitted inputs; "
        f"evidence: {output}"
    )
    return 0


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--release", type=Path, required=True)
    prepare_parser.add_argument("--archive-project-root", type=Path, required=True)
    prepare_parser.add_argument("--workspace", type=Path, required=True)
    prepare_parser.add_argument("--date", required=True)
    prepare_parser.add_argument("--flagnum", type=int, required=True)
    prepare_parser.add_argument(
        "--sma-discovery-relative", type=Path, required=True
    )
    prepare_parser.add_argument(
        "--archive-input-mode",
        choices=("copy", "symlink"),
        default="copy",
        help=(
            "copy gives independent archive inputs (recommended); symlink uses "
            "less disk but cannot provide OS-level write isolation"
        ),
    )
    init_parser = subparsers.add_parser("init")
    init_parser.add_argument("--workspace", type=Path, required=True)

    adjudicate_parser = subparsers.add_parser(
        "adjudicate",
        help=(
            "validate a hash-frozen reviewed current-to-generated delta ledger "
            "without rerunning CASA or modifying either table"
        ),
    )
    adjudicate_parser.add_argument("--workspace", type=Path, required=True)
    adjudicate_parser.add_argument("--run-id", required=True)
    adjudicate_parser.add_argument("--review-freeze", type=Path, required=True)
    adjudicate_parser.add_argument("--review-freeze-sha256", required=True)

    render_parser = subparsers.add_parser(
        "render-review-evidence",
        help=(
            "CASA-free no-clobber rendering of exact common-delta, product-"
            "addition, and hash-freeze evidence"
        ),
    )
    render_parser.add_argument("--workspace", type=Path, required=True)
    render_parser.add_argument("--run-id", required=True)

    for command in ("plan", "run"):
        subparser = subparsers.add_parser(command)
        subparser.add_argument("--workspace", type=Path, required=True)
        subparser.add_argument("--python-bin", type=Path, required=True)
        subparser.add_argument("--analysis-utils-dir", type=Path, required=True)
        subparser.add_argument("--date", required=True)
        subparser.add_argument("--flagnum", type=int, required=True)
        subparser.add_argument("--calnum", type=int, default=1)
        subparser.add_argument("--smapath", action="append", type=Path, required=True)
        subparser.add_argument("--selected-clean-manifest", type=Path)
        subparser.add_argument("--use-selected-clean", action="store_true")
        subparser.add_argument("--skip-reference-check", action="store_true")
        subparser.add_argument("--canonical-reference", type=Path)
        subparser.add_argument("--exclusion-policy", type=Path)
        subparser.add_argument("--run-id")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.command == "prepare":
        sma_path = prepare_workspace(
            release=args.release,
            archive_project_root=args.archive_project_root,
            workspace=args.workspace,
            date=args.date,
            flagnum=args.flagnum,
            sma_discovery_relative=args.sma_discovery_relative,
            archive_input_mode=args.archive_input_mode,
        )
        print(f"prepared and initialized isolated continuum workspace: {args.workspace}")
        print(f"isolated SMA discovery path: {sma_path}")
        return 0
    if args.command == "init":
        marker = initialize_workspace(args.workspace)
        print(f"initialized isolated continuum workspace: {marker}")
        return 0
    if args.command == "adjudicate":
        return adjudicate_workflow(args)
    if args.command == "render-review-evidence":
        return render_review_evidence(args)
    return run_workflow(args, dry_run=args.command == "plan")


if __name__ == "__main__":
    raise SystemExit(main())
