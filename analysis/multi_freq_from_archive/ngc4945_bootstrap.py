"""Authenticate and initialize historical metadata for a fresh source-18 run.

The public final45/84 records contain the complete historical27/66 subset.
Only raw CSV records are selected: no measurements, errors, classifications,
or old-value fields are calculated. The409 reference supplies product KEYS
only; its already-corrected flux/SNR values must never supply OLD fields.
Fresh source-18 fitting remains mandatory after this initialization.
"""
from __future__ import annotations

import ast
import csv
import hashlib
import io
import json
from pathlib import Path

from rerun_ngc4945_nuclear import CORRECTION_FIELDS, PHOTOMETRY_FIELDS

PHOTOMETRY = "ngc4945_nuclear_photometry.csv"
CORRECTIONS = "validated_continuum_corrections.csv"
CANONICAL = "fitsumfiles/historical/fitsummary_fixtochosen_sorted_409_reference.txt"
POLICY = "continuum_product_addition_policy.csv"
SNAPSHOT_DIR = "source18_bootstrap"
EVIDENCE_NAME = "initialization.json"
INITIAL_HISTORY_DIR = "source18_initial_history"
PROVENANCE_NAME = "public_source_provenance.json"
# Authenticated v1.0.10 public data, not hashes of the changing candidate code.
PUBLIC_SOURCE_SHA256 = {
    PHOTOMETRY: "365878f9099c075a8702bf158523782b62c265d50592e3a2e49df773057c3650",
    CORRECTIONS: "93a6095ed704fc94a421c6bda6f81ac5b5f52c16c8aecf33ca1f8864931e3f11",
    CANONICAL: "c1b3f8e5f3829fab4eecc8bdf8d0e537ee77d68ab89f711189f897f2f2bdf4f8",
}
# Independently freeze-bound historical artifact hashes, verification ONLY.
HISTORICAL_SHA256 = {
    PHOTOMETRY: "c54a119f577da7f3c28c4d3fd8c6270e6c50e891211a5f697d93c7968a2de270",
    CORRECTIONS: "89f7c157cc4888c64a26e22a667854e82a8c669364aca6441e7720e8d8353d91",
}
PUBLIC_POLICY_SHA256 = "903cf3813bfd56c8a24593e4898ff1ce7a7c866a755455646aaf6097e663aa5e"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def checked_file(root: Path, relative: str) -> Path:
    path = root / relative
    require(root.is_dir() and not root.is_symlink(), "source18 root is not a real directory")
    require(not Path(relative).is_absolute() and ".." not in Path(relative).parts,
            "source18 path is not relative")
    for item in (path, *path.parents):
        require(not item.is_symlink(), f"source18 symlink is forbidden: {item}")
        if item == root:
            break
    require(path.is_file(), f"source18 input is missing: {path}")
    return path


def csv_records(data: bytes):
    """Parse values while preserving each CSV record's original exact bytes."""
    text = data.decode("utf-8")
    require(not text.startswith("\ufeff"), "source18 CSV BOM is unsupported")
    physical = io.StringIO(text, newline="").readlines()
    reader = csv.reader(io.StringIO(text, newline=""), strict=True)
    records, previous = [], 0
    for values in reader:
        raw = "".join(physical[previous:reader.line_num]).encode("utf-8")
        previous = reader.line_num
        require(bool(values), "source18 blank CSV record")
        records.append((values, raw))
    require(bool(records) and previous == len(physical), "source18 incomplete CSV")
    require(b"".join(raw for _, raw in records) == data, "source18 raw CSV bytes changed")
    fields = tuple(records[0][0])
    require(len(set(fields)) == len(fields), "source18 duplicate CSV header")
    rows = []
    for values, raw in records[1:]:
        require(len(values) == len(fields), "source18 CSV width mismatch")
        rows.append((dict(zip(fields, values, strict=True)), raw))
    return fields, records[0][1], rows


def product(value: str) -> str:
    # Only409 embeds whitespace/newline around an otherwise portable product.
    value = value.strip()
    require(bool(value) and "\\" not in value and not value.startswith("/")
            and ":" not in value and ".." not in Path(value).parts,
            "source18 unsupported product identity")
    return value


def unique_products(rows, label):
    indexed = {}
    for row, _ in rows:
        key = product(row["product"])
        require(row["product"] == key, f"source18 noncanonical product in {label}")
        require(key not in indexed, f"source18 duplicate product in {label}: {key}")
        indexed[key] = row
    return indexed


def initial_history_provenance() -> dict:
    """Deterministic public-only derivation record, independent of fresh outputs."""
    return {
        "schema_version": 1,
        "public_release": "Zenodo v1.0.10 record22121490",
        "public_zip_sha256": "05672f6e083ddbea78375741390cb51fd17dc8186997a0d0086746a91dbd9734",
        "public_member_sha256": {**PUBLIC_SOURCE_SHA256, POLICY: PUBLIC_POLICY_SHA256},
        "historical_sha256": HISTORICAL_SHA256,
        "counts": {"historical_photometry": 27, "historical_corrections": 66,
                   "unchanged_other_source_corrections": 39, "new_NGC4945_products": 18},
        "selection": "409 NGC4945 product keys ONLY; retain raw CSV records in original order",
        "lineage": "OLD fields from public correction ledger, never409 fitted values",
        "serialization": "unchanged headers, values, quoting, newlines and relative record order",
        "qualification": "Historical metadata only; fresh source18 measurement remains mandatory",
    }


def provenance_bytes() -> bytes:
    return (json.dumps(initial_history_provenance(), indent=2, sort_keys=True) + "\n").encode()


def initial_history_payloads(root: Path) -> dict[str, bytes]:
    payloads = {name: checked_file(root, f"{INITIAL_HISTORY_DIR}/{name}").read_bytes()
                for name in HISTORICAL_SHA256}
    require({n: sha256(b) for n, b in payloads.items()} == HISTORICAL_SHA256,
            "source18 immutable initial-history asset changed")
    require(checked_file(root, f"{INITIAL_HISTORY_DIR}/{PROVENANCE_NAME}").read_bytes()
            == provenance_bytes(), "source18 initial-history public provenance changed")
    return payloads


def initialization_evidence(policy_bytes: bytes, release_metadata_sha256: dict) -> dict:
    require(isinstance(release_metadata_sha256, dict) and
            set(release_metadata_sha256) == {*PUBLIC_SOURCE_SHA256, POLICY} and
            all(isinstance(v, str) and len(v) == 64 and set(v) <= set('0123456789abcdef')
                for v in release_metadata_sha256.values()) and
            release_metadata_sha256[CANONICAL] == PUBLIC_SOURCE_SHA256[CANONICAL] and
            release_metadata_sha256[POLICY] == sha256(policy_bytes),
            "source18 release metadata hash record missing/stale")
    return {
        "schema_version": 1,
        "status": "historical_metadata_initialized_not_fitted_or_applied",
        "public_release": "Zenodo v1.0.10 record22121490",
        "initial_history_provenance_sha256": sha256(provenance_bytes()),
        "release_metadata_sha256": release_metadata_sha256,
        "addition_policy_sha256": sha256(policy_bytes),
        "initializer_sha256": sha256(Path(__file__).read_bytes()),
        "historical_sha256": HISTORICAL_SHA256,
        "counts": {"historical_photometry": 27, "historical_corrections": 66,
                   "unchanged_other_source_corrections": 39, "new_NGC4945_products": 18},
        "selection": "409 NGC4945 product keys ONLY; retain raw CSV records in original order",
        "lineage": "OLD fields from public correction ledger, never409 fitted values",
    }


def derive_initial_history_from_public(root: Path) -> dict[str, bytes]:
    """Release-build operation ONLY: derive assets from original public records."""
    payloads = {name: checked_file(root, name).read_bytes() for name in PUBLIC_SOURCE_SHA256}
    for name, expected in PUBLIC_SOURCE_SHA256.items():
        require(sha256(payloads[name]) == expected,
                f"source18 public source SHA256 mismatch: {name}")
    policy_bytes = checked_file(root, POLICY).read_bytes()
    require(sha256(policy_bytes) == PUBLIC_POLICY_SHA256, "source18 public admission SHA256 mismatch")
    pf, ph, pr = csv_records(payloads[PHOTOMETRY])
    cf, ch, cr = csv_records(payloads[CORRECTIONS])
    policy_fields, _, policy_rows = csv_records(policy_bytes)
    require(pf == PHOTOMETRY_FIELDS and cf == CORRECTION_FIELDS, "source18 schema mismatch")
    require({"source", "product", "expected_run_disposition"} <= set(policy_fields),
            "source18 admission schema mismatch")
    require(len(pr) == 45 and len(cr) == 84, "source18 final45/84 counts mismatch")
    phot, corr = unique_products(pr, "photometry"), unique_products(cr, "corrections")
    unique_products(policy_rows, "admissions")
    require(all(row["source"] == "NGC4945" for row, _ in pr), "source18 source mismatch")
    require(all(row["schema_version"] == "2" for row, _ in cr), "source18 correction version")
    require({k for k, row in corr.items() if row["source"] == "NGC4945"} == set(phot),
            "source18 final photometry/correction inventory mismatch")
    canonical = [ast.literal_eval(line) for line in payloads[CANONICAL].decode().splitlines()
                 if line.strip()]
    require(len(canonical) == 409 and all(isinstance(r, list) and len(r) == 8 for r in canonical),
            "source18 historical409 schema/count mismatch")
    require(len({product(r[7]) for r in canonical}) == 409, "source18 duplicate historical product")
    historical = {product(r[7]) for r in canonical if r[0] == "NGC4945"}
    admitted_rows = [r for r, _ in policy_rows if r["source"] == "NGC4945"]
    admitted = {product(r["product"]) for r in admitted_rows}
    require(len(historical) == 27 and len(admitted_rows) == len(admitted) == 18,
            "source18 historical27/admitted18 inventory mismatch")
    require(all(r["expected_run_disposition"] == "publication_row_required" for r in admitted_rows),
            "source18 admission disposition mismatch")
    require(not historical & admitted and historical | admitted == set(phot),
            "source18 historical27/admitted18 is not the exact disjoint final45 partition")
    for key in historical:
        for field in ("old_flux_field", "old_snr", "new_flux_field", "new_snr",
                      "source", "telescope", "frequency_ghz", "date_label"):
            require(phot[key][field] == corr[key][field], f"source18 lineage mismatch: {key} {field}")
        require(bool(corr[key]["old_flux_field"]) and bool(corr[key]["old_snr"]),
                "source18 blank historical old field")
    chosen_phot = [raw for r, raw in pr if product(r["product"]) in historical]
    chosen_corr = [raw for r, raw in cr if r["source"] != "NGC4945"
                   or product(r["product"]) in historical]
    require(len(chosen_phot) == 27 and len(chosen_corr) == 66, "source18 derived counts mismatch")
    derived = {PHOTOMETRY: ph + b"".join(chosen_phot), CORRECTIONS: ch + b"".join(chosen_corr)}
    require({n: sha256(b) for n, b in derived.items()} == HISTORICAL_SHA256,
            "source18 derived historical byte hashes mismatch")
    return derived


def materialize_initial_history(public_root: Path, destination: Path) -> None:
    """Build new immutable candidate assets; never overwrite source or outputs."""
    payloads = derive_initial_history_from_public(public_root)
    require(not destination.exists() and not destination.is_symlink(),
            "source18 initial-history destination already exists")
    destination.mkdir()
    for name, data in {**payloads, PROVENANCE_NAME: provenance_bytes()}.items():
        with (destination / name).open('xb') as stream:
            stream.write(data)


def derive_bootstrap(root: Path) -> tuple[dict[str, bytes], dict[str, object]]:
    """Prepare from immutable history, NOT from current final measurement values."""
    derived = initial_history_payloads(root)
    current = {name: checked_file(root, name).read_bytes()
               for name in (*PUBLIC_SOURCE_SHA256, POLICY)}
    require(sha256(current[CANONICAL]) == PUBLIC_SOURCE_SHA256[CANONICAL],
            "source18 historical409 reference changed")
    pf, _, pr = csv_records(current[PHOTOMETRY])
    cf, _, cr = csv_records(current[CORRECTIONS])
    _, _, admissions = csv_records(current[POLICY])
    require(pf == PHOTOMETRY_FIELDS and cf == CORRECTION_FIELDS,
            "source18 current final metadata schema mismatch")
    require(len(pr) == 45 and len(cr) == 84, "source18 current final45/84 counts mismatch")
    final_phot, final_corr = unique_products(pr, 'current photometry'), unique_products(cr, 'current corrections')
    unique_products(admissions, 'admissions')
    require(all(r['source'] == 'NGC4945' for r, _ in pr) and
            all(r['schema_version'] == '2' for r, _ in cr), 'source18 current final source/schema mismatch')
    require({k for k, r in final_corr.items() if r['source'] == 'NGC4945'} == set(final_phot),
            'source18 current final45 inventory mismatch')
    _, _, hp = csv_records(derived[PHOTOMETRY])
    _, _, hc = csv_records(derived[CORRECTIONS])
    historical_phot, historical_corr = unique_products(hp, 'initial photometry'), unique_products(hc, 'initial corrections')
    canonical = [ast.literal_eval(line) for line in current[CANONICAL].decode().splitlines() if line.strip()]
    historical_keys = {product(r[7]) for r in canonical if r[0] == 'NGC4945'}
    admitted_rows = [r for r, _ in admissions if r['source'] == 'NGC4945']
    admitted = {product(r['product']) for r in admitted_rows}
    require(set(historical_phot) == historical_keys and len(historical_keys) == 27 and
            len(admitted_rows) == len(admitted) == 18 and not historical_keys & admitted and
            historical_keys | admitted == set(final_phot),
            'source18 historical27/admitted18 is not the exact disjoint final45 partition')
    require(all(r['expected_run_disposition'] == 'publication_row_required' for r in admitted_rows),
            'source18 admission disposition mismatch')
    for key, row in final_phot.items():
        for field in ('source','telescope','frequency_ghz','date_label','old_flux_field','old_snr',
                      'new_flux_field','new_snr'):
            require(row[field] == final_corr[key][field], f'source18 current lineage mismatch: {key} {field}')
        if key in historical_keys:
            for field in ('source','telescope','frequency_ghz','date_label','old_flux_field','old_snr'):
                require(row[field] == historical_phot[key][field],
                        f'source18 historical OLD/identity changed: {key} {field}')
    for key, row in historical_corr.items():
        if row['source'] != 'NGC4945':
            require(final_corr.get(key) == row, f'source18 other-source lineage changed: {key}')
    hashes = {name: sha256(data) for name, data in current.items()}
    return derived, initialization_evidence(current[POLICY], hashes)


def initialize_bootstrap(root: Path, derived: dict[str, bytes], evidence: dict) -> dict:
    """Initialize ONLY the independently copied workspace, once, before fitting."""
    require({n: sha256(b) for n, b in derived.items()} == HISTORICAL_SHA256,
            "source18 initialization payload mismatch")
    captured_hashes = evidence.get('release_metadata_sha256')
    require(evidence == initialization_evidence(checked_file(root, POLICY).read_bytes(), captured_hashes),
            "source18 initialization evidence mismatch")
    # Recheck both destination files BEFORE any writes; never reset a prior run.
    initial_history_payloads(root)
    for name, expected in captured_hashes.items():
        require(sha256(checked_file(root, name).read_bytes()) == expected,
                f"source18 initialization requires unchanged release metadata after copy: {name}")
    snapshot = root / SNAPSHOT_DIR
    require(not snapshot.exists() and not snapshot.is_symlink(),
            "source18 bootstrap already exists; prepare a NEW workspace")
    snapshot.mkdir()
    for name, data in derived.items():
        with (snapshot / name).open("xb") as stream:
            stream.write(data)
    with (snapshot / EVIDENCE_NAME).open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    # copytree has made these ordinary independent files. Initialization does
    # not touch the release. A write failure leaves an uninitialized workspace,
    # never a marker/preparation manifest claiming success.
    for name, data in derived.items():
        with (root / name).open("wb") as stream:
            stream.write(data)
    record = {"schema_version": 1, "immutable_files": {
        f"{SNAPSHOT_DIR}/{name}": sha256(checked_file(root, f"{SNAPSHOT_DIR}/{name}").read_bytes())
        for name in (*HISTORICAL_SHA256, EVIDENCE_NAME)}}
    validate_preapplication_lineage(root, record)
    return record


def validate_bootstrap_snapshot(root: Path, record: dict) -> dict:
    """Validate immutable preparation evidence before AND after measurement."""
    require(isinstance(record, dict) and set(record) == {"schema_version", "immutable_files"}
            and record["schema_version"] == 1, "source18 bootstrap record missing/stale")
    expected_paths = {f"{SNAPSHOT_DIR}/{n}" for n in (*HISTORICAL_SHA256, EVIDENCE_NAME)}
    files = record["immutable_files"]
    require(isinstance(files, dict) and set(files) == expected_paths, "source18 snapshot inventory mismatch")
    for relative, expected in files.items():
        require(sha256(checked_file(root, relative).read_bytes()) == expected,
                f"source18 immutable snapshot changed: {relative}")
    for name, expected in HISTORICAL_SHA256.items():
        require(files[f"{SNAPSHOT_DIR}/{name}"] == expected, "source18 historical snapshot hash mismatch")
    evidence = json.loads(checked_file(root, f"{SNAPSHOT_DIR}/{EVIDENCE_NAME}").read_text())
    require(evidence == initialization_evidence(checked_file(root, POLICY).read_bytes(),
                                              evidence.get('release_metadata_sha256')),
            "source18 initialization provenance changed")
    initial_history_payloads(root)
    require(sha256(checked_file(root, CANONICAL).read_bytes()) == PUBLIC_SOURCE_SHA256[CANONICAL],
            "source18 historical409 reference changed")
    return record


def validate_preapplication_lineage(root: Path, record: dict) -> None:
    """Fail before any fitting if working metadata is final, corrupt or reused."""
    validate_bootstrap_snapshot(root, record)
    for name, expected in HISTORICAL_SHA256.items():
        require(sha256(checked_file(root, name).read_bytes()) == expected,
                f"source18 pre-application lineage is not authenticated27/66: {name}")
