#!/usr/bin/env python3
"""Run a script, -c command or -m module with an explicit isolated CASA runtime.

Normal use is through the generated ``wm-agn-casa-python`` shell wrapper,
which removes ambient Python/loader settings before starting Python itself.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import runpy
import sys


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bundle_inventory(path):
    root = Path(path)
    if root.is_symlink() or not root.is_dir():
        raise ValueError("analysisUtils must be a real directory")
    files = {}
    for item in sorted(root.rglob("*")):
        relative = item.relative_to(root)
        if "__pycache__" in relative.parts or item.suffix == ".pyc" or item.name.startswith("casa-"):
            continue
        if item.is_symlink():
            raise ValueError("analysisUtils bundle contains a symlink")
        if item.is_dir():
            continue
        if not item.is_file():
            raise ValueError("analysisUtils bundle contains a non-file")
        files[relative.as_posix()] = {"size": item.stat().st_size, "sha256": file_sha256(item)}
    if "analysisUtils.py" not in files or "mpfit.py" not in files:
        raise ValueError("analysisUtils.py and mpfit.py are required")
    return files


def build_environment(runtime, source=None):
    root = Path(runtime).resolve()
    environment = dict(os.environ if source is None else source)
    # Do not change HOME. These variables must not redirect CASA/Python into
    # the user's existing research environment or incompatible shared libs.
    for key in ("PYTHONPATH", "PYTHONHOME", "PYTHONSTARTUP", "LD_LIBRARY_PATH",
                "LD_PRELOAD", "ASDM_LIBRARY_PATH", "CASASITECONFIG"):
        environment.pop(key, None)
    environment.update({
        "WM_AGN_CASA_RUNTIME": str(root), "PYTHONNOUSERSITE": "1",
        "PYTHONDONTWRITEBYTECODE": "1", "CASARCFILES": str(root / "casarc"),
        "XDG_CONFIG_HOME": str(root / "config"), "XDG_CACHE_HOME": str(root / "cache"),
        "MPLCONFIGDIR": str(root / "matplotlib"), "IPYTHONDIR": str(root / "ipython"),
        "JUPYTER_CONFIG_DIR": str(root / "jupyter/config"),
        "JUPYTER_DATA_DIR": str(root / "jupyter/data"),
        "NUMBA_CACHE_DIR": str(root / "numba"),
    })
    return environment


def select_analysis_utils(runtime, manifest, explicit=None):
    expected = manifest["analysis_utils_files"]
    default = runtime / "analysis_utils/analysis_scripts"
    if bundle_inventory(default) != expected:
        raise RuntimeError("Installed analysisUtils bundle changed after runtime setup")
    requested = explicit or os.environ.get("WM_AGN_ANALYSIS_UTILS_DIR")
    if requested:
        selected = Path(requested).expanduser().resolve()
        if bundle_inventory(selected) != expected:
            raise RuntimeError("Explicit analysisUtils bundle differs from the pinned runtime")
        return selected
    # The four-pass orchestrator freezes the complete bundle and supplies its
    # location via PYTHONPATH. The shell captures, then clears that variable.
    # Use a frozen copy only when every non-cache file matches the pinned copy.
    inherited = os.environ.get("WM_AGN_REQUESTED_UTILS_PATH", "")
    if inherited and os.pathsep not in inherited:
        try:
            candidate = Path(inherited).expanduser().resolve()
            if bundle_inventory(candidate) == expected:
                return candidate
        except (OSError, ValueError):
            pass
    return default


def initialize_runtime(runtime, explicit_utils=None):
    root = Path(runtime).expanduser().resolve(strict=True)
    manifest = json.loads((root / "runtime_manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1 or manifest.get("state") != "configured":
        raise RuntimeError("Runtime setup is incomplete; preserve it and use a new directory")
    for name in ("casa_python.py", "isolated_casa_config.py", "casarc"):
        if file_sha256(root / name) != manifest["runtime_file_sha256"][name]:
            raise RuntimeError(f"Runtime configuration changed after setup: {name}")
    selected_utils = select_analysis_utils(root, manifest, explicit_utils)
    environment = build_environment(root, dict(os.environ))
    os.environ.clear()
    os.environ.update(environment)
    # CASA parses configuration from argv on its first import. Load only our
    # explicit config before target arguments or target imports are visible.
    saved_argv = list(sys.argv)
    sys.argv = [str(root / "casa_python.py"), "--nositeconfig", "--configfile",
                str(root / "isolated_casa_config.py")]
    try:
        if "casaconfig.config" in sys.modules:
            raise RuntimeError("CASA configuration was imported before isolated initialization")
        from casaconfig import config
    finally:
        sys.argv = saved_argv
    if config.load_failure():
        raise RuntimeError(f"CASA configuration failed: {config.load_failure()}")
    if (Path(config.measurespath).resolve() != root / "data"
            or [Path(p).resolve() for p in config.datapath] != [root / "data"]
            or config.data_auto_update or config.measures_auto_update):
        raise RuntimeError("CASA runtime escaped its explicit data/update settings")
    sys.path.insert(1 if sys.path else 0, str(selected_utils))
    return root, selected_utils


def dispatch_target(arguments):
    if not arguments:
        raise ValueError("Provide a script, -c CODE, or -m MODULE")
    if arguments[0] in ("-c", "-m"):
        if len(arguments) < 2:
            raise ValueError(f"{arguments[0]} requires an argument")
        sys.path[0] = os.getcwd()
        if arguments[0] == "-c":
            sys.argv = ["-c", *arguments[2:]]
            exec(compile(arguments[1], "<wm-agn-casa-command>", "exec"),
                 {"__name__": "__main__", "__file__": "<wm-agn-casa-command>"})
        else:
            sys.argv = [arguments[1], *arguments[2:]]
            runpy.run_module(arguments[1], run_name="__main__", alter_sys=True)
    else:
        target = Path(arguments[0]).expanduser().resolve(strict=True)
        if not target.is_file():
            raise ValueError("The CASA target must be a script file")
        sys.argv = [str(target), *arguments[1:]]
        sys.path[0] = str(target.parent)
        runpy.run_path(str(target), run_name="__main__")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime", type=Path, required=True)
    parser.add_argument("--analysis-utils-dir", type=Path)
    parser.add_argument("target", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    target = args.target[1:] if args.target[:1] == ["--"] else args.target
    if not target:
        parser.error("a script, -c CODE, or -m MODULE is required after --")
    initialize_runtime(args.runtime, args.analysis_utils_dir)
    dispatch_target(target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
