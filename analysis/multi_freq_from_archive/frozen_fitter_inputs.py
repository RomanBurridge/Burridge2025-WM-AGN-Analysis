#!/usr/bin/env python3
"""Validate and recover the exact historical continuum-fitter input set.

``frozen_fitter_inputs.csv`` freezes the 544 products presented to the Task 28
discovery pass (222 ALMA FITS files, 307 historical NVAS products, and 15 SMA
dirty-image CASA trees).  It also records the six selected-clean SMA
image/residual trees used by the fitter after discovery.  The latter are
required support products, but are not additional discovery inputs.

The module is deliberately independent of CASA.  Its public API lets the
four-pass runner load the exact whitelist and validate the requested archive
tree before creating a workspace.  The optional recovery command downloads
the 524 directly retrievable files and the five public IQUV parents from which
the historical Stokes-I products were derived.  An explicit optional flag runs
the reviewed deterministic extractor and installs those five outputs only when
their complete frozen file hashes match.  Recovery does not attempt to
reproduce the 21 SMA CASA products.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import re
import shutil
import sys
import tempfile
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Callable, Iterable, Mapping, Sequence
from urllib.parse import unquote, urlsplit


ANALYSIS_DIR = Path(__file__).resolve().parent
DEFAULT_MANIFEST = ANALYSIS_DIR / "frozen_fitter_inputs.csv"
EXPECTED_MANIFEST_SHA256 = (
    "2dd89ea3a30192efdcbcfa95ba3a297cfd59b8d26dd719197931102534c2ed4a"
)

MANIFEST_FIELDS = (
    "schema_version",
    "input_role",
    "discovery_member",
    "telescope",
    "source",
    "relative_path",
    "expected_sha256",
    "sha256_algorithm",
    "recovery_kind",
    "recovery_url",
    "recovery_sha256",
    "transformation",
    "archive_key",
)

EXPECTED_ROLE_COUNTS = {
    "archive_fitter_input": 529,
    "sma_discovery_input": 15,
    "sma_selected_clean_support": 6,
}
EXPECTED_TELESCOPE_COUNTS = {"ALMA": 222, "VLA": 307, "SMA": 21}
EXPECTED_DISCOVERY_COUNT = 544
EXPECTED_SUPPORT_COUNT = 6
EXPECTED_DIRECT_PUBLIC_FILES = 524
EXPECTED_DERIVED_PUBLIC_PARENTS = 5
ALLOWED_RECOVERY_HOSTS = {"almascience.eso.org", "www.vla.nrao.edu"}
SHA256_RE = re.compile(r"[0-9a-f]{64}")
BUFFER_SIZE = 16 * 1024 * 1024
STOKES_I_BUILDER = (
    ANALYSIS_DIR
    / "archive_inventory_audit/task28/build_approved_alma_inputs.py"
)


@dataclass(frozen=True)
class FrozenInputSelection:
    """Hash-validated paths selected from one requested archive tree."""

    manifest_sha256: str
    archive_root: Path
    discovery_paths: tuple[Path, ...]
    selected_clean_paths: tuple[Path, ...]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(BUFFER_SIZE):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_workflow_tree(path: Path) -> str:
    """Return the historical runner's deterministic CASA-tree digest.

    This is intentionally distinct from ``sha256-relpath-tree-v1`` used by the
    publication recovery manifest.  Task 28's runner sorted every relative
    entry, hashed its eight-byte big-endian path length and path, and then
    hashed ``D`` for a directory, ``F`` plus file bytes for a regular file, or
    ``L`` plus an eight-byte target length and the link target for a symlink.
    """

    if not path.is_dir() or path.is_symlink():
        raise FileNotFoundError(f"expected real directory: {path}")
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
                while chunk := handle.read(BUFFER_SIZE):
                    digest.update(chunk)
        elif item.is_dir():
            digest.update(b"D")
        else:
            raise ValueError(f"unsupported entry in CASA tree: {item}")
    return digest.hexdigest()


def manifest_sha256(path: Path = DEFAULT_MANIFEST) -> str:
    return sha256_file(path)


def _normalized_relative_path(value: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if (
        not value
        or path.is_absolute()
        or ".." in path.parts
        or "." in path.parts
        or str(path) != value
        or "\\" in value
    ):
        raise ValueError(f"unsafe or non-normalized relative path: {value!r}")
    return path


def _validate_https_url(value: str, *, allow_sma: bool = False) -> None:
    parsed = urlsplit(value)
    allowed = set(ALLOWED_RECOVERY_HOSTS)
    if allow_sma:
        allowed.add("lweb.cfa.harvard.edu")
    if (
        parsed.scheme != "https"
        or parsed.hostname not in allowed
        or parsed.username is not None
        or parsed.password is not None
        or parsed.fragment
    ):
        raise ValueError(f"unsupported recovery URL: {value!r}")


def validate_manifest_rows(rows: Sequence[Mapping[str, str]]) -> None:
    """Fail closed unless *rows* describe the exact frozen 550-row state."""

    if len(rows) != sum(EXPECTED_ROLE_COUNTS.values()):
        raise ValueError(f"expected 550 frozen rows, found {len(rows)}")

    role_counts = Counter(row.get("input_role", "") for row in rows)
    if dict(role_counts) != EXPECTED_ROLE_COUNTS:
        raise ValueError(
            f"wrong input-role counts: {dict(role_counts)}; "
            f"expected {EXPECTED_ROLE_COUNTS}"
        )
    telescope_counts = Counter(row.get("telescope", "") for row in rows)
    if dict(telescope_counts) != EXPECTED_TELESCOPE_COUNTS:
        raise ValueError(
            f"wrong telescope counts: {dict(telescope_counts)}; "
            f"expected {EXPECTED_TELESCOPE_COUNTS}"
        )

    discovery_count = sum(row.get("discovery_member") == "yes" for row in rows)
    support_count = sum(row.get("discovery_member") == "no" for row in rows)
    if discovery_count != EXPECTED_DISCOVERY_COUNT or support_count != EXPECTED_SUPPORT_COUNT:
        raise ValueError(
            f"wrong discovery/support counts: {discovery_count}/{support_count}"
        )

    seen: set[str] = set()
    direct_count = 0
    derived_count = 0
    for number, row in enumerate(rows, 2):
        missing = [field for field in MANIFEST_FIELDS if field not in row]
        if missing:
            raise ValueError(f"row {number}: missing columns {missing}")
        if row["schema_version"] != "1":
            raise ValueError(f"row {number}: unsupported schema version")
        relative = _normalized_relative_path(row["relative_path"])
        relative_text = relative.as_posix()
        if relative_text in seen:
            raise ValueError(f"row {number}: duplicate path {relative_text}")
        seen.add(relative_text)
        if not SHA256_RE.fullmatch(row["expected_sha256"]):
            raise ValueError(f"row {number}: invalid expected SHA-256")
        if not row["source"] or "/" in row["source"]:
            raise ValueError(f"row {number}: invalid source token")

        role = row["input_role"]
        telescope = row["telescope"]
        if role == "archive_fitter_input":
            if telescope not in {"ALMA", "VLA"} or relative.parts[0] != telescope:
                raise ValueError(f"row {number}: archive path/telescope mismatch")
            if len(relative.parts) < 3 or relative.parts[1] != row["source"]:
                raise ValueError(f"row {number}: archive path/source mismatch")
            if row["discovery_member"] != "yes":
                raise ValueError(f"row {number}: archive input must be a discovery member")
            if row["sha256_algorithm"] != "sha256-file-v1":
                raise ValueError(f"row {number}: archive file uses wrong hash algorithm")
        elif role in {"sma_discovery_input", "sma_selected_clean_support"}:
            if telescope != "SMA" or relative.parts[:2] != ("analysis", "sma_calibration"):
                raise ValueError(f"row {number}: SMA path/telescope mismatch")
            expected_member = "yes" if role == "sma_discovery_input" else "no"
            if row["discovery_member"] != expected_member:
                raise ValueError(f"row {number}: wrong SMA discovery membership")
            if row["sha256_algorithm"] != "sha256-workflow-tree-v1":
                raise ValueError(f"row {number}: SMA tree uses wrong hash algorithm")
            if not relative.name.startswith(row["source"] + "_test."):
                raise ValueError(f"row {number}: SMA path/source mismatch")
        else:
            raise ValueError(f"row {number}: unsupported input role {role!r}")

        recovery_kind = row["recovery_kind"]
        if recovery_kind == "direct_public_file":
            direct_count += 1
            _validate_https_url(row["recovery_url"])
            if row["recovery_sha256"] != row["expected_sha256"]:
                raise ValueError(f"row {number}: direct recovery hash differs from product")
            if row["transformation"] != "none":
                raise ValueError(f"row {number}: direct recovery has a transformation")
            recovery_name = PurePosixPath(
                unquote(urlsplit(row["recovery_url"]).path)
            ).name
            if recovery_name != relative.name:
                raise ValueError(f"row {number}: direct recovery filename mismatch")
        elif recovery_kind == "derive_from_public_parent":
            derived_count += 1
            _validate_https_url(row["recovery_url"])
            if telescope != "ALMA" or role != "archive_fitter_input":
                raise ValueError(f"row {number}: derived parent is not an ALMA input")
            if not SHA256_RE.fullmatch(row["recovery_sha256"]):
                raise ValueError(f"row {number}: invalid parent SHA-256")
            if row["transformation"] != "extract_singleton_stokes_i_from_iquv":
                raise ValueError(f"row {number}: unsupported derived transformation")
            output_name = relative.name
            parent_name = PurePosixPath(unquote(urlsplit(row["recovery_url"]).path)).name
            if output_name == parent_name or ".IQUV." not in parent_name:
                raise ValueError(f"row {number}: derived row does not identify its IQUV parent")
        elif recovery_kind == "sma_archive_reduction_product":
            _validate_https_url(row["recovery_url"], allow_sma=True)
            if telescope != "SMA":
                raise ValueError(f"row {number}: non-SMA reduction product")
            if row["recovery_sha256"]:
                raise ValueError(f"row {number}: SMA archive row has a file recovery hash")
            if row["transformation"] != "calibrate_and_image_with_released_sma_workflow":
                raise ValueError(f"row {number}: unsupported SMA recovery route")
            required_tokens = {
                "program=2022B-H002",
                "date=2023-03-16",
                "track=230316_128",
                "flagnum=2",
                "calnum=1",
            }
            if not required_tokens.issubset(set(row["archive_key"].split(";"))):
                raise ValueError(f"row {number}: incomplete SMA archive key")
        else:
            raise ValueError(f"row {number}: unsupported recovery kind {recovery_kind!r}")

    if direct_count != EXPECTED_DIRECT_PUBLIC_FILES:
        raise ValueError(f"expected 524 direct public files, found {direct_count}")
    if derived_count != EXPECTED_DERIVED_PUBLIC_PARENTS:
        raise ValueError(f"expected five derived-parent routes, found {derived_count}")


def load_frozen_fitter_inputs(
    manifest: Path = DEFAULT_MANIFEST,
) -> list[dict[str, str]]:
    """Load and strictly validate the public-safe frozen input manifest."""

    manifest = Path(manifest)
    if not manifest.is_file() or manifest.is_symlink():
        raise FileNotFoundError(f"manifest must be a real file: {manifest}")
    observed_manifest_sha256 = manifest_sha256(manifest)
    if observed_manifest_sha256 != EXPECTED_MANIFEST_SHA256:
        raise ValueError(
            f"frozen manifest hash mismatch: expected {EXPECTED_MANIFEST_SHA256}, "
            f"observed {observed_manifest_sha256}"
        )
    with manifest.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != MANIFEST_FIELDS:
            raise ValueError(
                f"{manifest}: unexpected columns {reader.fieldnames}; "
                f"expected {list(MANIFEST_FIELDS)}"
            )
        rows = [dict(row) for row in reader]
    validate_manifest_rows(rows)
    return rows


def _product_path(archive_root: Path, row: Mapping[str, str]) -> Path:
    root = archive_root.resolve(strict=False)
    relative = _normalized_relative_path(row["relative_path"])
    if row["input_role"] == "archive_fitter_input":
        candidate = root / "analysis/multi_freq_from_archive" / Path(*relative.parts)
    else:
        candidate = root / Path(*relative.parts)
    resolved = candidate.resolve(strict=False)
    if resolved == root or root not in resolved.parents:
        raise ValueError(f"input path escapes archive root: {row['relative_path']}")
    return candidate


def _verify_product(path: Path, row: Mapping[str, str]) -> str:
    algorithm = row["sha256_algorithm"]
    if algorithm == "sha256-file-v1":
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"expected real file: {path}")
        observed = sha256_file(path)
    elif algorithm == "sha256-workflow-tree-v1":
        observed = sha256_workflow_tree(path)
    else:
        raise ValueError(f"unsupported hash algorithm: {algorithm}")
    if observed != row["expected_sha256"]:
        raise ValueError(
            f"hash mismatch for {row['relative_path']}: "
            f"expected {row['expected_sha256']}, observed {observed}"
        )
    return observed


def validate_and_select_frozen_inputs(
    archive_root: Path,
    manifest: Path = DEFAULT_MANIFEST,
    *,
    workers: int = 4,
) -> FrozenInputSelection:
    """Hash every required product, then return the exact fitter whitelist."""

    if not 1 <= workers <= 16:
        raise ValueError("workers must be in the range 1..16")
    archive_root = Path(archive_root)
    if not archive_root.is_dir() or archive_root.is_symlink():
        raise FileNotFoundError(f"archive root must be a real directory: {archive_root}")
    root = archive_root.resolve()
    rows = load_frozen_fitter_inputs(Path(manifest))
    products = [(row, _product_path(root, row)) for row in rows]

    errors: list[str] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(_verify_product, path, row): row["relative_path"]
            for row, path in products
        }
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as exc:  # preserve all missing/mismatched identities
                errors.append(f"{futures[future]}: {exc}")
    if errors:
        preview = "\n".join(sorted(errors)[:20])
        more = "" if len(errors) <= 20 else f"\n... and {len(errors) - 20} more"
        raise ValueError(f"frozen input validation failed:\n{preview}{more}")

    discovery = tuple(
        path for row, path in products if row["discovery_member"] == "yes"
    )
    support = tuple(
        path for row, path in products if row["discovery_member"] == "no"
    )
    return FrozenInputSelection(
        manifest_sha256=manifest_sha256(Path(manifest)),
        archive_root=root,
        discovery_paths=discovery,
        selected_clean_paths=support,
    )


def _download_verified(
    url: str,
    expected: str,
    destination: Path,
    timeout: int,
    root: Path,
) -> str:
    root = _validate_recovery_destination(destination, root)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if destination.is_file() and sha256_file(destination) == expected:
            return "already_verified"
        raise FileExistsError(f"refusing to overwrite nonmatching file: {destination}")
    legacy_part = destination.with_name(destination.name + ".part")
    if legacy_part.exists() or legacy_part.is_symlink():
        raise FileExistsError(f"refusing existing partial download: {legacy_part}")
    request = urllib.request.Request(url, headers={"User-Agent": "WM-AGN-input-recovery/1"})
    # A unique private staging directory makes cleanup ownership unambiguous:
    # never unlink a shared .part path another process may have created.
    with tempfile.TemporaryDirectory(
        prefix=f".{destination.name}.download-", dir=destination.parent
    ) as temporary:
        part = Path(temporary) / "payload.part"
        with urllib.request.urlopen(request, timeout=timeout) as response, part.open("xb") as out:
            _validate_https_url(response.geturl())
            shutil.copyfileobj(response, out, BUFFER_SIZE)
            out.flush()
            os.fsync(out.fileno())
        observed = sha256_file(part)
        if observed != expected:
            raise ValueError(
                f"download hash mismatch for {url}: expected {expected}, observed {observed}"
            )
        _install_verified_no_clobber(part, destination)
    return "downloaded_verified"


def _install_verified_no_clobber(part: Path, destination: Path) -> None:
    """Atomically install a same-filesystem file without replacing any path.

    A final exists() check followed by replace()/rename() is not no-clobber:
    another process can create a file between those operations. link() fails
    atomically if *destination* exists, including a symlink or concurrent file.
    The source lives in this call's private same-directory staging area.
    """
    if part.is_symlink() or not part.is_file():
        raise ValueError(f"verified staging product is not a real file: {part}")
    os.link(part, destination, follow_symlinks=False)
    part.unlink()


def _validate_recovery_destination(destination: Path, root: Path) -> Path:
    """Return the real recovery root after rejecting path and symlink escapes."""

    root = root.resolve()
    resolved = destination.resolve(strict=False)
    if resolved == root or root not in resolved.parents:
        raise ValueError(f"recovery destination escapes archive root: {destination}")
    if destination.is_symlink():
        raise ValueError(f"refusing symlinked recovery destination: {destination}")
    parent = destination.parent
    while parent != root:
        if parent.is_symlink():
            raise ValueError(f"refusing symlinked recovery parent: {parent}")
        if root not in parent.resolve(strict=False).parents and parent != root:
            raise ValueError(f"recovery parent escapes archive root: {parent}")
        parent = parent.parent
    return root


def _load_stokes_i_extractor() -> Callable[[Path, Path], None]:
    """Load the reviewed deterministic Task 28 extractor only when requested."""

    if not STOKES_I_BUILDER.is_file() or STOKES_I_BUILDER.is_symlink():
        raise FileNotFoundError(f"reviewed Stokes-I builder is missing: {STOKES_I_BUILDER}")
    spec = importlib.util.spec_from_file_location(
        "wm_agn_task28_stokes_i_builder", STOKES_I_BUILDER
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load reviewed Stokes-I builder: {STOKES_I_BUILDER}")
    try:
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    except ImportError as exc:
        raise RuntimeError(
            "--derive-stokes-i requires the release's NumPy and Astropy dependencies"
        ) from exc
    extractor = getattr(module, "extract_stokes_i", None)
    if not callable(extractor):
        raise RuntimeError("reviewed Task 28 builder has no extract_stokes_i function")
    return extractor


def _derived_parent_path(root: Path, row: Mapping[str, str]) -> Path:
    parent_name = PurePosixPath(unquote(urlsplit(row["recovery_url"]).path)).name
    parent = root / ".frozen_input_parents" / row["source"] / parent_name
    resolved = parent.resolve(strict=False)
    if root not in resolved.parents:
        raise ValueError(f"derived-parent cache escapes recovery root: {parent}")
    return parent


def _derive_stokes_i_verified(
    parent: Path,
    parent_expected: str,
    destination: Path,
    output_expected: str,
    root: Path,
    extractor: Callable[[Path, Path], None],
) -> str:
    """Create one exact reviewed Stokes-I product without clobbering data."""

    _validate_recovery_destination(destination, root)
    if not parent.is_file() or parent.is_symlink():
        raise FileNotFoundError(f"expected real IQUV parent: {parent}")
    observed_parent = sha256_file(parent)
    if observed_parent != parent_expected:
        raise ValueError(
            f"IQUV parent hash mismatch for {parent}: expected {parent_expected}, "
            f"observed {observed_parent}"
        )
    if destination.exists() or destination.is_symlink():
        if (
            destination.is_file()
            and not destination.is_symlink()
            and sha256_file(destination) == output_expected
        ):
            return "already_verified_derived"
        raise FileExistsError(
            f"refusing to overwrite nonmatching derived product: {destination}"
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    legacy_part = destination.with_name(destination.name + ".stokes-i.part")
    if legacy_part.exists() or legacy_part.is_symlink():
        raise FileExistsError(f"refusing existing derived partial product: {legacy_part}")
    with tempfile.TemporaryDirectory(
        prefix=f".{destination.name}.stokes-i-", dir=destination.parent
    ) as temporary:
        part = Path(temporary) / "payload.part"
        extractor(parent, part)
        if not part.is_file() or part.is_symlink():
            raise RuntimeError(f"Stokes-I extractor did not create a real file: {part}")
        observed_output = sha256_file(part)
        if observed_output != output_expected:
            raise ValueError(
                f"derived Stokes-I hash mismatch for {destination}: "
                f"expected {output_expected}, observed {observed_output}"
            )
        _install_verified_no_clobber(part, destination)
    return "derived_verified"


def recover_public_sources(
    archive_root: Path,
    manifest: Path = DEFAULT_MANIFEST,
    *,
    workers: int = 4,
    timeout: int = 300,
    derive_stokes_i: bool = False,
) -> dict[str, object]:
    """Recover direct files and the five derived products' public parents.

    Direct products are placed at their exact archive-tree paths.  Public IQUV
    parents are placed under ``.frozen_input_parents``.  With
    ``derive_stokes_i=True``, the reviewed Task 28 extractor is loaded lazily,
    creates each historical Stokes-I product through a temporary file, and must
    reproduce its exact frozen byte hash before installation.  SMA reduction
    products remain an explicit separate archive/calibration step.
    """

    if not 1 <= workers <= 16:
        raise ValueError("workers must be in the range 1..16")
    if timeout < 1:
        raise ValueError("timeout must be positive")
    archive_root = Path(archive_root)
    if archive_root.exists() and (not archive_root.is_dir() or archive_root.is_symlink()):
        raise ValueError(f"recovery root must be a real directory: {archive_root}")
    archive_root.mkdir(parents=True, exist_ok=True)
    root = archive_root.resolve()
    rows = load_frozen_fitter_inputs(Path(manifest))

    jobs: list[tuple[str, str, Path, str]] = []
    for row in rows:
        kind = row["recovery_kind"]
        if kind == "direct_public_file":
            target = _product_path(root, row)
            jobs.append((row["relative_path"], row["recovery_url"], target, row["recovery_sha256"]))
        elif kind == "derive_from_public_parent":
            parent = _derived_parent_path(root, row)
            jobs.append((f"parent:{row['relative_path']}", row["recovery_url"], parent, row["recovery_sha256"]))

    statuses: Counter[str] = Counter()
    errors: list[str] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(_download_verified, url, expected, target, timeout, root): identity
            for identity, url, target, expected in jobs
        }
        for future in as_completed(futures):
            try:
                statuses[future.result()] += 1
            except Exception as exc:
                errors.append(f"{futures[future]}: {exc}")
    if errors:
        preview = "\n".join(sorted(errors)[:20])
        more = "" if len(errors) <= 20 else f"\n... and {len(errors) - 20} more"
        raise RuntimeError(f"public recovery failed:\n{preview}{more}")

    derived_statuses: Counter[str] = Counter()
    if derive_stokes_i:
        extractor = _load_stokes_i_extractor()
        for row in rows:
            if row["recovery_kind"] != "derive_from_public_parent":
                continue
            status = _derive_stokes_i_verified(
                _derived_parent_path(root, row),
                row["recovery_sha256"],
                _product_path(root, row),
                row["expected_sha256"],
                root,
                extractor,
            )
            derived_statuses[status] += 1

    return {
        "manifest_sha256": manifest_sha256(Path(manifest)),
        "archive_root": str(root),
        "direct_historical_products_available": EXPECTED_DIRECT_PUBLIC_FILES,
        "derived_public_parents_available": EXPECTED_DERIVED_PUBLIC_PARENTS,
        "download_status_counts": dict(sorted(statuses.items())),
        "derive_stokes_i_requested": derive_stokes_i,
        "derived_output_status_counts": dict(sorted(derived_statuses.items())),
        "exact_public_fits_available": (
            EXPECTED_DIRECT_PUBLIC_FILES + EXPECTED_DERIVED_PUBLIC_PARENTS
            if derive_stokes_i
            else EXPECTED_DIRECT_PUBLIC_FILES
        ),
        "exact_derived_outputs_still_required": (
            0 if derive_stokes_i else EXPECTED_DERIVED_PUBLIC_PARENTS
        ),
        "sma_reduction_products_still_required": 21,
        "exact_historical_archive_complete": False,
    }


def _summary(rows: Sequence[Mapping[str, str]], manifest: Path) -> dict[str, object]:
    return {
        "manifest": str(manifest),
        "manifest_sha256": manifest_sha256(manifest),
        "rows": len(rows),
        "discovery_inputs": sum(row["discovery_member"] == "yes" for row in rows),
        "selected_clean_support": sum(row["discovery_member"] == "no" for row in rows),
        "roles": dict(sorted(Counter(row["input_role"] for row in rows).items())),
        "telescopes": dict(sorted(Counter(row["telescope"] for row in rows).items())),
        "recovery_kinds": dict(sorted(Counter(row["recovery_kind"] for row in rows).items())),
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("check-manifest", help="validate the frozen CSV only")
    validate = subparsers.add_parser("validate", help="hash the exact requested archive tree")
    validate.add_argument("--archive-root", type=Path, required=True)
    validate.add_argument("--workers", type=int, default=4)
    recover = subparsers.add_parser(
        "recover-public",
        help="download direct public files and the five public IQUV parents",
    )
    recover.add_argument("--archive-root", type=Path, required=True)
    recover.add_argument("--workers", type=int, default=4)
    recover.add_argument("--timeout", type=int, default=300)
    recover.add_argument(
        "--derive-stokes-i",
        action="store_true",
        help="rebuild and exact-hash-verify the five reviewed ALMA Stokes-I products",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    rows = load_frozen_fitter_inputs(args.manifest)
    if args.command == "check-manifest":
        payload = _summary(rows, args.manifest)
    elif args.command == "validate":
        selection = validate_and_select_frozen_inputs(
            args.archive_root, args.manifest, workers=args.workers
        )
        payload = {
            **_summary(rows, args.manifest),
            "archive_root": str(selection.archive_root),
            "validated_discovery_paths": len(selection.discovery_paths),
            "validated_selected_clean_paths": len(selection.selected_clean_paths),
            "status": "exact_historical_inputs_hash_verified",
        }
    elif args.command == "recover-public":
        payload = recover_public_sources(
            args.archive_root,
            args.manifest,
            workers=args.workers,
            timeout=args.timeout,
            derive_stokes_i=args.derive_stokes_i,
        )
    else:  # pragma: no cover - argparse makes this unreachable
        raise AssertionError(args.command)
    json.dump(payload, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
