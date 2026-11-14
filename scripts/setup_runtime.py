#!/usr/bin/env python3
"""Prepare a NEW isolated CASA runtime after installing the pinned dependencies.

Run with the installed CASA-capable Python. The optional --pull-data invokes
the official casaconfig CLI in this runtime only; automatic updates remain
disabled. Existing and failed runtime directories are never overwritten.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path, PurePosixPath
import shlex
import shutil
import subprocess
import sys
import tarfile
from urllib.request import urlopen

from casa_python import build_environment, bundle_inventory, file_sha256


ANALYSIS_UTILS_URL = "https://zenodo.org/api/records/22085689/files/analysis_scripts.tar/content"
ANALYSIS_UTILS_SHA256 = "63d5a314f7c892e2d0bec48dec6449440a1e2779a4d0ebf602aef94d13d31dcc"
ANALYSIS_UTILS_BYTES = 7188480
PINNED_PACKAGES = {"casaconfig": "1.5.0", "casatools": "6.7.6.14", "casatasks": "6.7.6.14"}


def safe_extract(archive_path, destination):
    """Extract regular analysis_scripts members only, without overwrite/links."""
    destination = Path(destination)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"Extraction directory already exists: {destination}")
    with tarfile.open(archive_path, "r:*") as archive:
        members = archive.getmembers()
        names = set()
        total = 0
        for member in members:
            relative = PurePosixPath(member.name)
            if (relative.is_absolute() or ".." in relative.parts or not relative.parts
                    or relative.parts[0] != "analysis_scripts" or "\\" in member.name
                    or relative.as_posix() != member.name.rstrip("/")
                    or member.name.rstrip("/") in names
                    or not (member.isfile() or member.isdir())):
                raise ValueError(f"Unsafe or duplicate analysisUtils tar member: {member.name}")
            names.add(member.name.rstrip("/"))
            total += member.size
            if member.size > 20_000_000 or total > 50_000_000:
                raise ValueError("analysisUtils tar exceeds the expected size bounds")
        if not {"analysis_scripts/analysisUtils.py", "analysis_scripts/mpfit.py"}.issubset(names):
            raise ValueError("analysisUtils tar lacks required source files")
        destination.mkdir(mode=0o700)
        for member in members:
            output = destination / member.name
            if member.isdir():
                output.mkdir(parents=True, exist_ok=False)
            else:
                output.parent.mkdir(parents=True, exist_ok=True)
                source = archive.extractfile(member)
                if source is None:
                    raise ValueError(f"Cannot read tar member: {member.name}")
                with source, output.open("xb") as handle:
                    shutil.copyfileobj(source, handle)
                output.chmod(0o644)


def download_archive(destination):
    """Download the fixed official archive, keeping partial failures as evidence."""
    with urlopen(ANALYSIS_UTILS_URL, timeout=60) as response, Path(destination).open("xb") as handle:
        count = 0
        while chunk := response.read(1024 * 1024):
            count += len(chunk)
            if count > ANALYSIS_UTILS_BYTES:
                raise ValueError("Downloaded analysisUtils exceeds its pinned byte count")
            handle.write(chunk)
    if count != ANALYSIS_UTILS_BYTES or file_sha256(destination) != ANALYSIS_UTILS_SHA256:
        raise ValueError("Downloaded analysisUtils does not match the pinned archive")


def installed_versions():
    if sys.version_info[:2] != (3, 12):
        raise RuntimeError("The pinned CASA environment requires CPython 3.12")
    versions = {name: importlib.metadata.version(name) for name in PINNED_PACKAGES}
    if versions != PINNED_PACKAGES:
        raise RuntimeError(f"CASA packages differ from the tested pins: {versions}")
    return versions


def config_text():
    # casaconfig evaluates this file in a defaults namespace, so __file__ is
    # not this config's location. The wrapper supplies an explicit root env.
    return '''"""Generated explicit configuration; never uses the user's ~/.casa."""
from os import environ as _environ
from pathlib import Path as _Path
_root = _Path(_environ["WM_AGN_CASA_RUNTIME"])
measurespath = str(_root / "data")
datapath = [measurespath]
data_auto_update = False
measures_auto_update = False
cachedir = str(_root / "cache_casa")
logfile = str(_root / "logs/casa.log")
startupfile = ""
user_site = False
nologfile = True
nologger = True
nogui = True
agg = True
iplog = False
log2term = True
casaconfig_verbose = 2
'''


def wrapper_text(python_executable):
    return '''#!/bin/sh
set -eu
wm_agn_runtime_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
export WM_AGN_CASA_RUNTIME="$wm_agn_runtime_dir"
export WM_AGN_REQUESTED_UTILS_PATH="${PYTHONPATH:-}"
unset PYTHONPATH PYTHONHOME PYTHONSTARTUP LD_LIBRARY_PATH LD_PRELOAD ASDM_LIBRARY_PATH CASASITECONFIG
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
exec ''' + shlex.quote(str(python_executable)) + ''' -B "$wm_agn_runtime_dir/casa_python.py" --runtime "$wm_agn_runtime_dir" -- "$@"
'''


def prepare_runtime(runtime, *, analysis_utils_archive=None, pull_data=False):
    root = Path(runtime).expanduser().absolute()
    if root.exists() or root.is_symlink():
        raise FileExistsError(f"Refusing existing runtime (including partial failures): {root}")
    if not root.parent.is_dir():
        raise ValueError("Create the runtime's parent directory first")
    versions = installed_versions()
    cached = None if analysis_utils_archive is None else Path(analysis_utils_archive).expanduser()
    if cached is not None and (cached.is_symlink() or not cached.is_file()
                               or cached.stat().st_size != ANALYSIS_UTILS_BYTES
                               or file_sha256(cached) != ANALYSIS_UTILS_SHA256):
        raise ValueError("Cached analysisUtils archive differs from the pinned original")
    root.mkdir(mode=0o700)
    status_path = root / "setup_status.json"
    status = {"schema_version": 1, "state": "building",
              "started_utc": datetime.now(timezone.utc).isoformat()}
    status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    try:
        for relative in ("downloads", "data", "cache", "config", "cache_casa", "logs",
                         "matplotlib", "ipython", "jupyter/config", "jupyter/data", "numba"):
            (root / relative).mkdir(parents=True, exist_ok=False)
        archive = root / "downloads/analysis_scripts_2.398.tar"
        if cached is not None:
            shutil.copyfile(cached, archive)
        else:
            download_archive(archive)
        if file_sha256(archive) != ANALYSIS_UTILS_SHA256:
            raise ValueError("Copied analysisUtils archive checksum mismatch")
        safe_extract(archive, root / "analysis_utils")
        inventory = bundle_inventory(root / "analysis_utils/analysis_scripts")
        shutil.copyfile(Path(__file__).with_name("casa_python.py"), root / "casa_python.py")
        (root / "isolated_casa_config.py").write_text(config_text(), encoding="utf-8")
        (root / "casarc").write_text(f"# Explicit runtime resource file.\nmeasures.directory: {root / 'data'}\n", encoding="utf-8")
        wrapper = root / "wm-agn-casa-python"
        wrapper.write_text(wrapper_text(Path(sys.executable).absolute()), encoding="utf-8")
        wrapper.chmod(0o755)
        pull_command = None
        if pull_data:
            pull_command = [sys.executable, "-B", "-m", "casaconfig", "--nositeconfig",
                            "--configfile", str(root / "isolated_casa_config.py"), "--pull-data"]
            with (root / "logs/pull_data.log").open("x", encoding="utf-8") as log:
                subprocess.run(pull_command, cwd=root, env=build_environment(root),
                               stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                               check=True, text=True)
            if not (root / "data/readme.txt").is_file():
                raise RuntimeError("casaconfig returned success but no CASA data readme was installed")
        data_inventory = {p.relative_to(root / "data").as_posix():
                          {"size": p.stat().st_size, "sha256": file_sha256(p)}
                          for p in sorted((root / "data").rglob("*")) if p.is_file()}
        (root / "data_inventory.json").write_text(json.dumps(data_inventory, indent=2) + "\n", encoding="utf-8")
        manifest = {"schema_version": 1, "state": "configured", "python_executable": sys.executable,
                    "python_version": sys.version, "packages": versions,
                    "analysis_utils_version": "2.398", "analysis_utils_url": ANALYSIS_UTILS_URL,
                    "analysis_utils_archive_sha256": ANALYSIS_UTILS_SHA256,
                    "analysis_utils_files": inventory, "pull_data_requested": pull_data,
                    "pull_data_command": pull_command, "data_file_count": len(data_inventory),
                    "data_inventory_sha256": file_sha256(root / "data_inventory.json"),
                    "automatic_updates": False,
                    "runtime_file_sha256": {name: file_sha256(root / name) for name in
                                            ("casa_python.py", "isolated_casa_config.py", "casarc")}}
        (root / "runtime_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        status.update(state="configured", completed_utc=datetime.now(timezone.utc).isoformat())
        status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        return root
    except Exception as exc:
        status.update(state="failed", error=f"{type(exc).__name__}: {exc}")
        status_path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime", type=Path, required=True, help="New runtime directory; must not already exist")
    parser.add_argument("--analysis-utils-archive", type=Path, help="Optional cached, checksum-verified official tar")
    parser.add_argument("--pull-data", action="store_true", help="Populate this runtime via the official casaconfig CLI")
    args = parser.parse_args(argv)
    root = prepare_runtime(args.runtime, analysis_utils_archive=args.analysis_utils_archive,
                           pull_data=args.pull_data)
    print(f"Configured isolated runtime: {root}")
    print(f"CASA Python launcher: {root / 'wm-agn-casa-python'}")
    if not args.pull_data:
        print("CASA reference data were not downloaded; use --pull-data during a new setup for fitting.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
