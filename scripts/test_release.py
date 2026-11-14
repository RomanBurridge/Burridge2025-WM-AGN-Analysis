#!/usr/bin/env python3
"""Run every public source test in its validated Python environment.

Scientific, table and workflow tests run with the explicitly supplied fitting
Python.  The two figure-artifact modules run inside the explicitly supplied
locked Ubuntu renderer, because running those byte checks with a newer
Matplotlib/FreeType stack produces misleading failures.  No failure is waived:
any nonzero invocation, missing module, or source-tree drift makes this command
return nonzero.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

try:  # Direct script execution and package import are both supported.
    from .render_publication import command_prefix
except ImportError:  # pragma: no cover - exercised by direct CLI execution
    from render_publication import command_prefix


RELEASE_ROOT = Path(__file__).resolve().parents[1]
TEST_PACKAGES = (
    (
        "analysis.CountingBHS_and_Spectra.tests",
        Path("analysis/CountingBHS_and_Spectra/tests"),
    ),
    (
        "analysis.multi_freq_from_archive.tests",
        Path("analysis/multi_freq_from_archive/tests"),
    ),
)
RENDERER_MODULES = frozenset(
    {
        "analysis.multi_freq_from_archive.tests.test_build_bhs_flux_figure",
        "analysis.multi_freq_from_archive.tests.test_build_sed_figures",
    }
)
EXPECTED_RENDERER = {
    "python": "3.12.3",
    "matplotlib": "3.6.3",
    "numpy": "1.26.4",
    "scipy": "1.11.4",
    "pillow": "10.2.0",
    "freetype": "2.13.2",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def discover_test_modules(release: Path) -> tuple[list[str], list[str]]:
    """Return complete fitting/renderer inventories from the two test dirs."""

    modules: list[str] = []
    for package, relative in TEST_PACKAGES:
        directory = release / relative
        if not directory.is_dir() or directory.is_symlink():
            raise FileNotFoundError(f"test directory must be real: {directory}")
        files = sorted(directory.glob("test_*.py"), key=lambda path: path.name)
        if not files:
            raise ValueError(f"no test modules found in {directory}")
        for path in files:
            if not path.is_file() or path.is_symlink():
                raise ValueError(f"test module must be a real file: {path}")
            modules.append(f"{package}.{path.stem}")
    if len(modules) != len(set(modules)):
        raise ValueError("duplicate test module identity")
    observed_renderer = set(modules) & RENDERER_MODULES
    if observed_renderer != RENDERER_MODULES:
        raise ValueError(
            "renderer test inventory changed: "
            f"expected {sorted(RENDERER_MODULES)}, observed {sorted(observed_renderer)}"
        )
    fitting = [module for module in modules if module not in RENDERER_MODULES]
    renderer = [module for module in modules if module in RENDERER_MODULES]
    if len(fitting) + len(renderer) != len(modules):
        raise AssertionError("test routing lost a module")
    return fitting, renderer


def _excluded_snapshot_path(relative: Path) -> bool:
    return (
        "__pycache__" in relative.parts
        or relative.suffix in {".pyc", ".pyo"}
        or (relative.name.startswith("casa-") and relative.suffix == ".log")
        or relative.name == "casa.log"
    )


def source_snapshot(release: Path) -> dict[str, str]:
    """Hash every stable release file so tests cannot silently edit sources."""

    snapshot: dict[str, str] = {}
    for path in sorted(
        release.rglob("*"), key=lambda item: item.relative_to(release).as_posix()
    ):
        relative = path.relative_to(release)
        if _excluded_snapshot_path(relative):
            continue
        if path.is_symlink():
            snapshot[relative.as_posix()] = "symlink:" + os.readlink(path)
        elif path.is_file():
            snapshot[relative.as_posix()] = sha256_file(path)
        elif not path.is_dir():
            raise ValueError(f"unsupported release entry: {path}")
    return snapshot


def _run_logged(
    *,
    label: str,
    command: Sequence[str],
    cwd: Path,
    environment: dict[str, str],
    log_dir: Path,
) -> dict[str, object]:
    start = utc_now()
    stdout_path = log_dir / f"{label}.stdout.log"
    stderr_path = log_dir / f"{label}.stderr.log"
    exit_code = 127
    error: str | None = None
    try:
        with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
            completed = subprocess.run(
                list(command),
                cwd=cwd,
                env=environment,
                stdin=subprocess.DEVNULL,
                stdout=stdout,
                stderr=stderr,
                check=False,
            )
        exit_code = completed.returncode
    except OSError as exc:
        error = f"{type(exc).__name__}: {exc}"
        stderr_path.write_text(error + "\n", encoding="utf-8")
    return {
        "label": label,
        "command": list(command),
        "cwd": str(cwd),
        "start_utc": start,
        "end_utc": utc_now(),
        "exit_code": exit_code,
        "launch_error": error,
        "stdout": str(stdout_path),
        "stderr": str(stderr_path),
    }


def _validate_paths(
    release: Path,
    fitting_python: Path,
    renderer_rootfs: Path,
    output: Path,
) -> tuple[Path, Path, Path, Path]:
    release = release.resolve()
    renderer_rootfs = renderer_rootfs.resolve()
    fitting_python = fitting_python.absolute()
    output = output.absolute()
    if not (release / "analysis").is_dir() or not (release / "scripts").is_dir():
        raise ValueError(f"invalid release root: {release}")
    if not fitting_python.is_file() or not os.access(fitting_python, os.X_OK):
        raise ValueError(f"fitting Python is not executable: {fitting_python}")
    if not (renderer_rootfs / "usr/bin/python3").is_file():
        raise ValueError(f"invalid renderer rootfs: {renderer_rootfs}")
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"output already exists: {output}")
    if output.is_relative_to(release) or release.is_relative_to(output):
        raise ValueError("output must be new and external to the release")
    if output.is_relative_to(renderer_rootfs) or renderer_rootfs.is_relative_to(output):
        raise ValueError("output must be external to the renderer rootfs")
    return release, fitting_python, renderer_rootfs, output


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fitting-python", type=Path, required=True)
    parser.add_argument("--renderer-rootfs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--release", type=Path, default=RELEASE_ROOT)
    return parser.parse_args(argv)


def fitting_test_command(
    fitting_python: Path, modules: Sequence[str]
) -> list[str]:
    """Build the argv accepted by both Python and wm-agn-casa-python.

    The CASA wrapper accepts a script, ``-c`` or ``-m`` as its first
    argument; bytecode suppression is supplied through
    ``PYTHONDONTWRITEBYTECODE`` in the child environment instead of ``-B``.
    """

    return [
        str(fitting_python),
        "-m",
        "unittest",
        "-v",
        *modules,
    ]


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    release, fitting_python, renderer_rootfs, output = _validate_paths(
        args.release, args.fitting_python, args.renderer_rootfs, args.output
    )
    fitting_modules, renderer_modules = discover_test_modules(release)
    before = source_snapshot(release)

    for relative in (
        "logs",
        "fitting_runtime/cache",
        "fitting_runtime/config",
        "renderer_runtime/runtime/cache",
        "renderer_runtime/runtime/config",
    ):
        (output / relative).mkdir(parents=True, exist_ok=True)
    log_dir = output / "logs"

    fitting_environment = dict(os.environ)
    fitting_environment.update(
        {
            "PYTHONPATH": "",
            "PYTHONNOUSERSITE": "1",
            "PYTHONDONTWRITEBYTECODE": "1",
            "XDG_CACHE_HOME": str(output / "fitting_runtime/cache"),
            "XDG_CONFIG_HOME": str(output / "fitting_runtime/config"),
        }
    )
    fitting_command = fitting_test_command(fitting_python, fitting_modules)
    results = [
        _run_logged(
            label="fitting_tests",
            command=fitting_command,
            cwd=release,
            environment=fitting_environment,
            log_dir=log_dir,
        )
    ]

    renderer_output = output / "renderer_runtime"
    prefix = command_prefix(renderer_rootfs, release, renderer_output)
    probe = (
        "import json,sys,matplotlib,numpy,scipy,PIL; "
        "from matplotlib import ft2font; "
        "actual={'python':'.'.join(map(str,sys.version_info[:3])),"
        "'matplotlib':matplotlib.__version__,'numpy':numpy.__version__,"
        "'scipy':scipy.__version__,'pillow':PIL.__version__,"
        "'freetype':ft2font.__freetype_version__}; "
        f"expected={EXPECTED_RENDERER!r}; "
        "print(json.dumps(actual,sort_keys=True)); "
        "assert actual==expected, (actual,expected)"
    )
    renderer_environment = {"PATH": os.environ.get("PATH", "")}
    results.append(
        _run_logged(
            label="renderer_probe",
            command=[*prefix, "-c", probe],
            cwd=release,
            environment=renderer_environment,
            log_dir=log_dir,
        )
    )
    results.append(
        _run_logged(
            label="renderer_tests",
            command=[*prefix, "-m", "unittest", "-v", *renderer_modules],
            cwd=release,
            environment=renderer_environment,
            log_dir=log_dir,
        )
    )

    after = source_snapshot(release)
    changed = sorted(
        set(before) ^ set(after)
        | {key for key in before.keys() & after if before[key] != after[key]}
    )
    passed = all(result["exit_code"] == 0 for result in results) and not changed
    report = {
        "schema_version": 1,
        "status": "PASS" if passed else "FAIL",
        "release": str(release),
        "fitting_python": str(fitting_python),
        "renderer_rootfs": str(renderer_rootfs),
        "output": str(output),
        "fitting_test_modules": fitting_modules,
        "renderer_test_modules": renderer_modules,
        "module_counts": {
            "fitting": len(fitting_modules),
            "renderer": len(renderer_modules),
            "total": len(fitting_modules) + len(renderer_modules),
        },
        "failure_policy": (
            "No expected-failure or return-code waiver. A stale reviewed-release "
            "adjudication is a real failure until regenerated by the actual workflow."
        ),
        "results": results,
        "source_snapshot_before_sha256": hashlib.sha256(
            json.dumps(before, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest(),
        "source_snapshot_after_sha256": hashlib.sha256(
            json.dumps(after, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest(),
        "source_changes": changed,
    }
    (output / "TEST_RELEASE_REPORT.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
