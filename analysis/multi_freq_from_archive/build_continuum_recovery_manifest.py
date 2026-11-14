#!/usr/bin/env python3
"""Build and verify the portable continuum-product recovery manifest.

The 501-row final continuum table mixes fitted archive products with literature
measurements.  This builder gives every row a stable identity and enough
portable provenance to recover the fitted input (or identify the literature
source) without embedding a workstation path.

Structural generation uses only lightweight release files.  Supplying
``--archive-root`` additionally resolves every fitted product below a local
copy of the historical project tree and recomputes its recorded size and
SHA-256 digest.  CASA image directories use the documented
``sha256-relpath-tree-v1`` algorithm rather than filesystem metadata.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import os
import re
import sys
from collections import Counter
from functools import lru_cache
from pathlib import Path, PurePosixPath
from typing import Iterable, Mapping, Sequence


ANALYSIS_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = ANALYSIS_DIR.parents[1]
FINAL_TABLE = ANALYSIS_DIR / "fitsumfiles/final/fitsummary_final_withextras.txt"
CHECKSUMS = ANALYSIS_DIR / "continuum_product_checksums.csv"
MANIFEST = ANALYSIS_DIR / "machinetables/continuum_recovery_manifest.csv"
SUMMARY = ANALYSIS_DIR / "machinetables/continuum_recovery_manifest_summary.json"
TASK28_ALMA_INPUTS = (
    ANALYSIS_DIR
    / "archive_inventory_audit/task28/approved_alma_input_manifest.csv"
)

EXPECTED_ROWS = 501
ALMA_URL = "https://almascience.org/aq/"
NVAS_URL = "https://www.vla.nrao.edu/astro/nvas/"
SMA_URL = "https://lweb.cfa.harvard.edu/rtdc/SMAdata/archives/"

SMA_TRACKS = {
    "230316_128": {
        "program": "2022B-H002",
        "date": "2023-03-16",
    },
}

# One Cycle-1 image was distributed under its execution-block filename rather
# than the later ``member.uid___...`` product convention.  The FITS HISTORY
# identifies project 2012.1.00912.S and execution block A002/X5b1929/X57d; the
# ALMA archive hierarchy for that project identifies the Member OUS below.
ALMA_LEGACY_PRODUCTS = {
    "uid___A002_X5b1929_X57d.continuum.fits": {
        "project": "2012.1.00912.S",
        "member_ous_uid": "uid://A002/X5a9a13/X528",
        "execution_block_uid": "uid://A002/X5b1929/X57d",
    },
}

LITERATURE = {
    "Kamali, 2017, A&A, 605:A84": {
        "key": "doi:10.1051/0004-6361/201730899",
        "url": "https://doi.org/10.1051/0004-6361/201730899",
    },
    "Doi et al. 2013, ApJ, 765:63, Tables 1-2; reprocessed VLA data": {
        "key": "doi:10.1088/0004-637X/765/1/63",
        "url": "https://doi.org/10.1088/0004-637X/765/1/63",
    },
    "Doi et al. 2013, ApJ, 765:63, Table 2; nuclear-region flux-density upper bound": {
        "key": "doi:10.1088/0004-637X/765/1/63",
        "url": "https://doi.org/10.1088/0004-637X/765/1/63",
    },
    "Doi et al, 2005, MNRAS, 363:2": {
        "key": "doi:10.1111/j.1365-2966.2005.09471.x",
        "url": "https://doi.org/10.1111/j.1365-2966.2005.09471.x",
    },
    "Humphreys et al. 2005, ApJL, 634:L133": {
        "key": "doi:10.1086/498890",
        "url": "https://doi.org/10.1086/498890",
    },
    "Lin et al. 2026, ApJ accepted, arXiv:2604.27159": {
        "key": "arxiv:2604.27159",
        "url": "https://arxiv.org/abs/2604.27159",
    },
}

MANIFEST_FIELDS = (
    "record_id",
    "source",
    "telescope",
    "frequency_ghz",
    "date_label",
    "provenance_kind",
    "observation_id",
    "archive_key",
    "product_filename",
    "product_or_citation",
    "retrieval_url",
    "local_product_kind",
    "local_size_bytes",
    "sha256",
    "sha256_algorithm",
    "recovery_status",
)

CHECKSUM_FIELDS = (
    "record_id",
    "product_path",
    "local_product_kind",
    "local_size_bytes",
    "sha256",
    "sha256_algorithm",
)


def _read_final_table(path: Path = FINAL_TABLE) -> list[list[str]]:
    rows: list[list[str]] = []
    for line_number, raw_line in enumerate(path.read_text().splitlines(), 1):
        if not raw_line.strip():
            continue
        try:
            parsed = ast.literal_eval(raw_line.strip())
        except (SyntaxError, ValueError) as exc:
            raise ValueError(f"{path}:{line_number}: invalid Python-list row") from exc
        if not isinstance(parsed, list) or len(parsed) < 8:
            raise ValueError(f"{path}:{line_number}: expected at least eight fields")
        rows.append([str(value).strip() for value in parsed[:8]])
    if len(rows) != EXPECTED_ROWS:
        raise ValueError(
            f"expected {EXPECTED_ROWS} canonical continuum rows, found {len(rows)}"
        )
    return rows


def _validate_relative_product(value: str) -> str:
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or str(path) != value:
        raise ValueError(f"product path is not a normalized release-relative path: {value}")
    if value.startswith(("/home/", "~/")):
        raise ValueError(f"workstation path is forbidden in release metadata: {value}")
    return value


def stable_record_id(row: Sequence[str]) -> str:
    """Return an order-independent identity for one continuum observation row."""

    identity = [row[index].strip() for index in (0, 1, 2, 4, 7)]
    payload = json.dumps(identity, ensure_ascii=True, separators=(",", ":"))
    return "continuum-" + hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20]


def _alma_uid(filename: str) -> tuple[str, str]:
    match = re.match(r"^(member\.uid___([A-Za-z0-9]+(?:_[A-Za-z0-9]+)+))\.", filename)
    if match is None:
        raise ValueError(f"cannot extract ALMA Member OUS UID from {filename}")
    token = match.group(1)
    official = "uid://" + "/".join(match.group(2).split("_"))
    return token, official


def _execution_uid(filename: str) -> str:
    """Return the official UID encoded by a legacy execution-block filename."""

    match = re.match(r"^(uid___([A-Za-z0-9]+(?:_[A-Za-z0-9]+)+))\.", filename)
    if match is None:
        raise ValueError(f"cannot extract ALMA execution-block UID from {filename}")
    return "uid://" + "/".join(match.group(2).split("_"))


@lru_cache(maxsize=1)
def _load_task28_alma_mappings(
    path: Path = TASK28_ALMA_INPUTS,
) -> dict[str, dict[str, str]]:
    """Load approved products whose filenames encode an execution block, not a MOUS.

    Most ALMA filenames directly contain a ``member.uid___...`` token.  Seven
    approved Task 28 products instead retain older ``uid___A002...`` delivery
    names.  Their official Member OUS identities are therefore taken from the
    hash-validated Task 28 input manifest rather than guessed from the filename.
    """

    required = {
        "member_ous_uid",
        "proposal_ids",
        "output_relative_path",
    }
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing required columns {sorted(missing)}")
        rows = list(reader)

    mappings: dict[str, dict[str, str]] = {}
    for row in rows:
        product = _validate_relative_product(row["output_relative_path"])
        filename = PurePosixPath(product).name
        if filename.startswith("member.uid___"):
            continue
        member_ous_uid = row["member_ous_uid"].strip()
        proposal_ids = row["proposal_ids"].strip()
        if not re.fullmatch(r"uid://[A-Za-z0-9]+(?:/[A-Za-z0-9]+)+", member_ous_uid):
            raise ValueError(f"{path}: invalid Member OUS UID for {product}")
        if not proposal_ids:
            raise ValueError(f"{path}: missing proposal ID for {product}")
        if product in mappings:
            raise ValueError(f"{path}: duplicate output product {product}")
        mappings[product] = {
            "project": proposal_ids,
            "member_ous_uid": member_ous_uid,
            "execution_block_uid": _execution_uid(filename),
        }
    return mappings


def _vla_program(filename: str) -> str:
    parts = filename.split("_")
    if len(parts) < 4 or not parts[1]:
        raise ValueError(f"cannot extract NVAS program token from {filename}")
    return parts[1]


def _classify_row(row: Sequence[str]) -> dict[str, str]:
    source, telescope, frequency, _flux, date_label, _beam, _snr, reference = row
    record_id = stable_record_id(row)
    base = {
        "record_id": record_id,
        "source": source,
        "telescope": telescope,
        "frequency_ghz": frequency,
        "date_label": date_label,
    }

    if telescope == "ALMA" and reference.startswith("ALMA/"):
        product = _validate_relative_product(reference)
        filename = PurePosixPath(product).name
        if filename in ALMA_LEGACY_PRODUCTS:
            legacy = ALMA_LEGACY_PRODUCTS[filename]
            official_uid = legacy["member_ous_uid"]
            archive_key = (
                f"project={legacy['project']};member_ous_uid={official_uid};"
                f"execution_block_uid={legacy['execution_block_uid']}"
            )
        elif (legacy := _load_task28_alma_mappings().get(product)) is not None:
            official_uid = legacy["member_ous_uid"]
            archive_key = (
                f"project={legacy['project']};member_ous_uid={official_uid};"
                f"execution_block_uid={legacy['execution_block_uid']}"
            )
        else:
            member_token, official_uid = _alma_uid(filename)
            archive_key = (
                f"member_ous_uid={official_uid};product_token={member_token}"
            )
        return {
            **base,
            "provenance_kind": "alma_archive_fits",
            "observation_id": official_uid,
            "archive_key": archive_key,
            "product_filename": filename,
            "product_or_citation": product,
            "retrieval_url": ALMA_URL,
            "local_product_kind": "file",
        }

    if telescope == "VLA" and reference.startswith("VLA/"):
        product = _validate_relative_product(reference)
        filename = PurePosixPath(product).name
        program = _vla_program(filename)
        return {
            **base,
            "provenance_kind": "nvas_image_product",
            "observation_id": f"NVAS:{program}:{date_label}",
            "archive_key": f"program={program};filename={filename}",
            "product_filename": filename,
            "product_or_citation": product,
            "retrieval_url": NVAS_URL,
            "local_product_kind": "file",
        }

    if telescope == "SMA" and reference.startswith("analysis/sma_calibration/"):
        product = _validate_relative_product(reference)
        parts = PurePosixPath(product).parts
        if len(parts) < 4:
            raise ValueError(f"SMA product path is incomplete: {product}")
        track = parts[2]
        if track not in SMA_TRACKS:
            raise ValueError(f"SMA track lacks an explicit archive mapping: {track}")
        track_info = SMA_TRACKS[track]
        program = track_info["program"]
        observing_date = track_info["date"]
        return {
            **base,
            "provenance_kind": "sma_casa_image",
            "observation_id": f"SMA:{program}:{observing_date}:{track}",
            "archive_key": (
                f"program={program};date={observing_date};track={track};"
                f"image_path={product}"
            ),
            "product_filename": PurePosixPath(product).name,
            "product_or_citation": product,
            "retrieval_url": SMA_URL,
            "local_product_kind": "directory",
        }

    if reference not in LITERATURE:
        raise ValueError(
            f"unrecognized non-product provenance for {source}/{telescope}: {reference}"
        )
    citation = LITERATURE[reference]
    return {
        **base,
        "provenance_kind": "literature_measurement",
        "observation_id": citation["key"],
        "archive_key": citation["key"],
        "product_filename": "",
        "product_or_citation": reference,
        "retrieval_url": citation["url"],
        "local_product_kind": "none",
    }


def build_structural_rows(final_table: Path = FINAL_TABLE) -> list[dict[str, str]]:
    rows = [_classify_row(row) for row in _read_final_table(final_table)]
    ids = [row["record_id"] for row in rows]
    if len(set(ids)) != len(ids):
        duplicates = sorted(key for key, count in Counter(ids).items() if count > 1)
        raise ValueError(f"stable record IDs are not unique: {duplicates}")
    return rows


def _read_csv_by_id(path: Path, fields: Sequence[str]) -> dict[str, dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise ValueError(
                f"{path}: unexpected columns {reader.fieldnames}; expected {list(fields)}"
            )
        rows = list(reader)
    keyed: dict[str, dict[str, str]] = {}
    for row in rows:
        record_id = row["record_id"]
        if record_id in keyed:
            raise ValueError(f"{path}: duplicate record_id {record_id}")
        keyed[record_id] = row
    return keyed


def _write_csv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _sha256_file(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            size += len(chunk)
            digest.update(chunk)
    return size, digest.hexdigest()


def _directory_entries(root: Path) -> list[tuple[str, Path]]:
    entries: list[tuple[str, Path]] = []
    for directory, dirnames, filenames in os.walk(root, followlinks=False):
        dirnames.sort()
        filenames.sort()
        directory_path = Path(directory)
        for name in dirnames + filenames:
            path = directory_path / name
            relative = path.relative_to(root).as_posix()
            entries.append((relative, path))
    return sorted(entries, key=lambda item: item[0].encode("utf-8"))


def _sha256_directory(path: Path) -> tuple[int, str]:
    """Hash a directory from sorted relative names, types, and file bytes.

    Timestamps, owners, permissions, inode numbers, and the directory's absolute
    location are deliberately excluded.  Directory and symlink entries are
    included so that empty directories and link targets remain reproducible.
    """

    digest = hashlib.sha256()
    digest.update(b"sha256-relpath-tree-v1\0")
    size = 0
    for relative, child in _directory_entries(path):
        encoded_name = relative.encode("utf-8")
        if child.is_symlink():
            target = os.readlink(child).encode("utf-8")
            digest.update(b"L\0" + encoded_name + b"\0" + target + b"\0")
        elif child.is_dir():
            digest.update(b"D\0" + encoded_name + b"\0")
        elif child.is_file():
            file_size, file_digest = _sha256_file(child)
            size += file_size
            digest.update(
                b"F\0"
                + encoded_name
                + b"\0"
                + str(file_size).encode("ascii")
                + b"\0"
                + bytes.fromhex(file_digest)
            )
        else:
            raise ValueError(f"unsupported filesystem entry in CASA tree: {child}")
    return size, digest.hexdigest()


def hash_product(path: Path, expected_kind: str) -> tuple[int, str, str]:
    if expected_kind == "file":
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"expected regular file: {path}")
        size, digest = _sha256_file(path)
        return size, digest, "sha256-file-v1"
    if expected_kind == "directory":
        if not path.is_dir() or path.is_symlink():
            raise FileNotFoundError(f"expected directory: {path}")
        size, digest = _sha256_directory(path)
        return size, digest, "sha256-relpath-tree-v1"
    raise ValueError(f"unsupported local product kind: {expected_kind}")


def _resolve_below_archive(archive_root: Path, product: str) -> Path:
    archive_root = archive_root.resolve()
    relative = PurePosixPath(_validate_relative_product(product))
    if relative.parts[0] in {"ALMA", "VLA"}:
        candidate = archive_root / "analysis/multi_freq_from_archive" / Path(*relative.parts)
    elif relative.parts[:2] == ("analysis", "sma_calibration"):
        candidate = archive_root / Path(*relative.parts)
    else:
        raise ValueError(f"unsupported release product path: {product}")
    resolved = candidate.resolve()
    if archive_root not in resolved.parents:
        raise ValueError(f"resolved product escapes archive root: {product}")
    return resolved


def compute_checksums(
    structural_rows: Sequence[Mapping[str, str]], archive_root: Path
) -> list[dict[str, str]]:
    checksums: list[dict[str, str]] = []
    for number, row in enumerate(structural_rows, 1):
        if row["provenance_kind"] == "literature_measurement":
            continue
        product = row["product_or_citation"]
        local_path = _resolve_below_archive(archive_root, product)
        size, digest, algorithm = hash_product(local_path, row["local_product_kind"])
        checksums.append(
            {
                "record_id": row["record_id"],
                "product_path": product,
                "local_product_kind": row["local_product_kind"],
                "local_size_bytes": str(size),
                "sha256": digest,
                "sha256_algorithm": algorithm,
            }
        )
        if number % 50 == 0:
            print(f"Hashed {number}/{len(structural_rows)} manifest rows...", file=sys.stderr)
    return checksums


def merge_manifest(
    structural_rows: Sequence[Mapping[str, str]],
    checksums: Mapping[str, Mapping[str, str]],
) -> list[dict[str, str]]:
    local_ids = {
        row["record_id"]
        for row in structural_rows
        if row["provenance_kind"] != "literature_measurement"
    }
    if set(checksums) != local_ids:
        missing = sorted(local_ids - set(checksums))
        extra = sorted(set(checksums) - local_ids)
        raise ValueError(
            "checksum coverage does not exactly match fitted products: "
            f"missing={missing[:5]} extra={extra[:5]}"
        )

    output: list[dict[str, str]] = []
    for row in structural_rows:
        merged = dict(row)
        if row["provenance_kind"] == "literature_measurement":
            merged.update(
                {
                    "local_size_bytes": "",
                    "sha256": "",
                    "sha256_algorithm": "not_applicable",
                    "recovery_status": "literature_only_no_local_product",
                }
            )
        else:
            recorded = checksums[row["record_id"]]
            for key in ("product_path", "local_product_kind"):
                expected = (
                    row["product_or_citation"]
                    if key == "product_path"
                    else row["local_product_kind"]
                )
                if recorded[key] != expected:
                    raise ValueError(
                        f"{row['record_id']}: checksum {key}={recorded[key]!r}; "
                        f"expected {expected!r}"
                    )
            if not re.fullmatch(r"[0-9a-f]{64}", recorded["sha256"]):
                raise ValueError(f"{row['record_id']}: invalid SHA-256 digest")
            if int(recorded["local_size_bytes"]) <= 0:
                raise ValueError(f"{row['record_id']}: local size must be positive")
            expected_algorithm = (
                "sha256-file-v1"
                if row["local_product_kind"] == "file"
                else "sha256-relpath-tree-v1"
            )
            if recorded["sha256_algorithm"] != expected_algorithm:
                raise ValueError(
                    f"{row['record_id']}: unexpected checksum algorithm "
                    f"{recorded['sha256_algorithm']}"
                )
            merged.update(
                {
                    "local_size_bytes": recorded["local_size_bytes"],
                    "sha256": recorded["sha256"],
                    "sha256_algorithm": recorded["sha256_algorithm"],
                    "recovery_status": "local_product_hash_verified",
                }
            )
        if set(merged) != set(MANIFEST_FIELDS):
            raise AssertionError(f"internal manifest field mismatch: {set(merged)}")
        output.append({field: merged[field] for field in MANIFEST_FIELDS})
    return output


def _sha256_path(path: Path) -> str:
    return _sha256_file(path)[1]


def build_summary(
    rows: Sequence[Mapping[str, str]],
    *,
    final_table: Path = FINAL_TABLE,
    checksum_path: Path = CHECKSUMS,
) -> dict[str, object]:
    kinds = Counter(row["provenance_kind"] for row in rows)
    statuses = Counter(row["recovery_status"] for row in rows)
    local_rows = [row for row in rows if row["local_product_kind"] != "none"]
    return {
        "schema_version": 1,
        "row_count": len(rows),
        "unique_record_ids": len({row["record_id"] for row in rows}),
        "provenance_kind_counts": dict(sorted(kinds.items())),
        "recovery_status_counts": dict(sorted(statuses.items())),
        "local_product_count": len(local_rows),
        "literature_only_count": len(rows) - len(local_rows),
        "local_product_total_bytes": sum(int(row["local_size_bytes"]) for row in local_rows),
        "all_local_products_have_recorded_hashes": all(row["sha256"] for row in local_rows),
        "contains_absolute_home_paths": any(
            "/home/" in value
            for row in rows
            for value in row.values()
        ),
        "final_table_sha256": _sha256_path(final_table),
        "checksum_lock_sha256": _sha256_path(checksum_path),
        "directory_hash_algorithm": (
            "sha256-relpath-tree-v1: SHA-256 over sorted relative directory, "
            "symlink, and regular-file entries; regular files contribute size "
            "and SHA-256; filesystem metadata and absolute paths are excluded"
        ),
    }


def _csv_bytes(fields: Sequence[str], rows: Sequence[Mapping[str, object]]) -> bytes:
    import io

    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def _summary_bytes(summary: Mapping[str, object]) -> bytes:
    return (json.dumps(summary, indent=2, sort_keys=True) + "\n").encode("utf-8")


def _assert_exact(path: Path, expected: bytes) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)
    actual = path.read_bytes()
    if actual != expected:
        raise ValueError(f"generated file is stale or non-deterministic: {path}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="compare regenerated lightweight artifacts with checked-in files",
    )
    parser.add_argument(
        "--archive-root",
        type=Path,
        help=(
            "local Burridge2025-WM-AGN-Analysis project root; when supplied, "
            "recompute all fitted-product sizes and hashes"
        ),
    )
    parser.add_argument(
        "--refresh-checksums",
        action="store_true",
        help="rewrite the checksum lock from --archive-root before building",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.refresh_checksums and args.archive_root is None:
        raise ValueError("--refresh-checksums requires --archive-root")
    if args.refresh_checksums and args.check:
        raise ValueError("--refresh-checksums and --check are mutually exclusive")

    structural = build_structural_rows()
    if args.archive_root is not None:
        recomputed_rows = compute_checksums(structural, args.archive_root)
        recomputed_by_id = {row["record_id"]: row for row in recomputed_rows}
        if args.refresh_checksums:
            _write_csv(CHECKSUMS, CHECKSUM_FIELDS, recomputed_rows)
        else:
            recorded = _read_csv_by_id(CHECKSUMS, CHECKSUM_FIELDS)
            if recorded != recomputed_by_id:
                mismatches = sorted(
                    record_id
                    for record_id in set(recorded) | set(recomputed_by_id)
                    if recorded.get(record_id) != recomputed_by_id.get(record_id)
                )
                raise ValueError(
                    "archive verification differs from checksum lock for "
                    f"{len(mismatches)} records: {mismatches[:5]}"
                )

    checksums = _read_csv_by_id(CHECKSUMS, CHECKSUM_FIELDS)
    manifest_rows = merge_manifest(structural, checksums)
    summary = build_summary(manifest_rows)
    manifest_bytes = _csv_bytes(MANIFEST_FIELDS, manifest_rows)
    summary_bytes = _summary_bytes(summary)

    if args.check:
        _assert_exact(MANIFEST, manifest_bytes)
        _assert_exact(SUMMARY, summary_bytes)
        print(
            f"OK: {len(manifest_rows)} continuum records; "
            f"{summary['local_product_count']} hash-verified local products and "
            f"{summary['literature_only_count']} literature-only measurements."
        )
        return 0

    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_bytes(manifest_bytes)
    SUMMARY.write_bytes(summary_bytes)
    print(f"Wrote {MANIFEST}")
    print(f"Wrote {SUMMARY}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
