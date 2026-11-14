"""Portable handling of the selected SMA clean-image manifest.

The continuum fitter changes its working directory to this package directory.
Paths stored in ``selected_clean_iterations.csv`` are instead relative to the
repository root.  Keeping resolution here makes the behavior testable without
importing CASA and prevents a selected clean image from silently falling back
to a dirty ``nitercal_0`` product.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Callable


DEFAULT_REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def normalize_source_name(source_name: object, killstrings=None) -> str:
    if source_name is None:
        return ""

    cleaned_source = str(source_name).strip().split(".")[0]
    if killstrings is None:
        killstrings = ["Galaxy"]

    for killstring in killstrings:
        cleaned_source = cleaned_source.replace(killstring, "")

    if cleaned_source and cleaned_source[-1].isalpha():
        cleaned_source = cleaned_source[:-1] + cleaned_source[-1].upper()
    return cleaned_source


def resolve_manifest_path(
    manifest_path: str | Path,
    repository_root: str | Path = DEFAULT_REPOSITORY_ROOT,
) -> Path:
    """Resolve an absolute, current-directory, or repository-relative manifest."""

    raw_path = Path(manifest_path).expanduser()
    if raw_path.is_absolute():
        return raw_path.resolve()

    repository_candidate = Path(repository_root).resolve() / raw_path
    working_candidate = Path.cwd() / raw_path
    for candidate in (repository_candidate, working_candidate):
        if candidate.exists():
            return candidate.resolve()
    # Repository-relative is the documented interpretation when neither exists.
    return repository_candidate.resolve()


def resolve_manifest_entry_path(
    entry_path: str | Path,
    manifest_path: str | Path,
    repository_root: str | Path = DEFAULT_REPOSITORY_ROOT,
) -> Path:
    """Resolve an image path relative to the repository root or manifest.

    Entries beginning with ``analysis/`` use the repository root, matching the
    released manifest.  Other relative entries are first interpreted beside
    the manifest, which is convenient for a relocated user-supplied manifest.
    Existing paths win when both interpretations are possible.
    """

    raw_path = Path(str(entry_path).strip()).expanduser()
    if raw_path.is_absolute():
        return raw_path.resolve()

    manifest = Path(manifest_path).resolve()
    repository_candidate = Path(repository_root).resolve() / raw_path
    manifest_candidate = manifest.parent / raw_path

    normalized_parts = tuple(part for part in raw_path.parts if part != ".")
    if normalized_parts and normalized_parts[0] == "analysis":
        candidates = (repository_candidate, manifest_candidate)
    else:
        candidates = (manifest_candidate, repository_candidate)

    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return candidates[0].resolve()


def load_selected_clean_manifest(
    manifest_path: str | Path,
    date_code: str = "",
    flag_number: int | None = None,
    repository_root: str | Path = DEFAULT_REPOSITORY_ROOT,
) -> tuple[dict[str, str], list[str]]:
    """Load selected images for one track/flag combination."""

    manifest = resolve_manifest_path(manifest_path, repository_root)
    selected_lookup: dict[str, str] = {}
    selected_original_sources: list[str] = []

    if not manifest.exists():
        raise FileNotFoundError(
            f"Selected-clean manifest not found: {manifest}. Refusing to run "
            "without the release's clean-image selection record."
        )

    with manifest.open(newline="") as manifest_file:
        reader = csv.DictReader(manifest_file)
        for row in reader:
            image_path = (row.get("image_path") or "").strip()
            if not image_path:
                continue

            row_date = (row.get("date_code") or "").strip()
            if date_code and row_date and row_date != str(date_code):
                continue

            row_flagnum = (row.get("flagnum") or "").strip()
            if flag_number is not None and row_flagnum:
                try:
                    if int(float(row_flagnum)) != int(flag_number):
                        continue
                except ValueError:
                    continue

            source_name = (row.get("source") or "").strip()
            if not source_name:
                source_name = Path(image_path).name.replace("_test.image", "")

            normalized_source = normalize_source_name(source_name)
            if not normalized_source:
                continue

            resolved_image = resolve_manifest_entry_path(
                image_path,
                manifest,
                repository_root,
            )
            selected_lookup[normalized_source] = str(resolved_image)
            selected_original_sources.append(
                Path(image_path).name.replace("_test.image", "")
            )

    return selected_lookup, list(dict.fromkeys(selected_original_sources))


def selected_image_path(
    selected_lookup: dict[str, str],
    source: str,
    exists: Callable[[str], bool] | None = None,
) -> str | None:
    """Return a selected image or fail instead of silently using nitercal_0."""

    image_path = selected_lookup.get(normalize_source_name(source))
    if image_path is None:
        return None

    if exists is None:
        exists = lambda value: Path(value).exists()
    if not exists(image_path):
        raise FileNotFoundError(
            "The selected-clean manifest requires an image that is missing: "
            f"{image_path}. Retrieve or regenerate that clean image; refusing "
            "to substitute an nitercal_0 upper-limit product."
        )
    return image_path
