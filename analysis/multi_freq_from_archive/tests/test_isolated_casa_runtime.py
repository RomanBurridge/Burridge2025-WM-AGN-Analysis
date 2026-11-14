from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest import mock

SCRIPTS = Path(__file__).resolve().parents[3] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import casa_python as bootstrap
import setup_runtime as setup


def make_archive(path, extra=None):
    with tarfile.open(path, "w") as archive:
        for name in ("analysis_scripts/analysisUtils.py", "analysis_scripts/mpfit.py"):
            data = b"# synthetic dependency fixture\n"
            info = tarfile.TarInfo(name)
            info.size = len(data)
            archive.addfile(info, io.BytesIO(data))
        if extra is not None:
            archive.addfile(extra, io.BytesIO(b"x") if extra.size else None)


class RuntimeSafetyTests(unittest.TestCase):
    def test_safe_regular_tar_extracts_once(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive, output = root / "input.tar", root / "output"
            make_archive(archive)
            setup.safe_extract(archive, output)
            self.assertTrue((output / "analysis_scripts/analysisUtils.py").is_file())
            with self.assertRaises(FileExistsError):
                setup.safe_extract(archive, output)

    def test_unsafe_or_duplicate_archive_rejected_before_destination_creation(self):
        for name, kind in (("/outside", tarfile.REGTYPE), ("analysis_scripts/../outside", tarfile.REGTYPE),
                           ("analysis_scripts/link", tarfile.SYMTYPE), ("analysis_scripts/link", tarfile.LNKTYPE),
                           ("analysis_scripts/analysisUtils.py", tarfile.REGTYPE)):
            with self.subTest(name=name, kind=kind), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                extra = tarfile.TarInfo(name)
                extra.type = kind
                extra.linkname = "/outside"
                extra.size = 0
                make_archive(root / "input.tar", extra)
                with self.assertRaises(ValueError):
                    setup.safe_extract(root / "input.tar", root / "output")
                self.assertFalse((root / "output").exists())

    def test_existing_runtime_refused_before_dependency_or_network_checks(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "user_file").write_text("keep")
            with mock.patch.object(setup, "installed_versions", side_effect=AssertionError("should not run")):
                with self.assertRaises(FileExistsError):
                    setup.prepare_runtime(root)
            self.assertEqual((root / "user_file").read_text(), "keep")

    def test_small_mocked_setup_and_failed_pull_preserve_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = root / "source.tar"
            make_archive(archive)
            patches = [mock.patch.object(setup, "ANALYSIS_UTILS_SHA256", bootstrap.file_sha256(archive)),
                       mock.patch.object(setup, "ANALYSIS_UTILS_BYTES", archive.stat().st_size),
                       mock.patch.object(setup, "installed_versions", return_value=setup.PINNED_PACKAGES)]
            with contextlib.ExitStack() as stack:
                for patch in patches:
                    stack.enter_context(patch)
                configured = setup.prepare_runtime(root / "runtime", analysis_utils_archive=archive)
                manifest = json.loads((configured / "runtime_manifest.json").read_text())
                self.assertEqual(manifest["state"], "configured")
                self.assertFalse(manifest["pull_data_requested"])
                self.assertTrue(os.access(configured / "wm-agn-casa-python", os.X_OK))
                self.assertEqual((configured / "isolated_casa_config.py").read_text(), setup.config_text())
                with mock.patch.object(setup.subprocess, "run", side_effect=subprocess.CalledProcessError(1, ["fake-pull"])):
                    with self.assertRaises(subprocess.CalledProcessError):
                        setup.prepare_runtime(root / "failed", analysis_utils_archive=archive, pull_data=True)
                status = json.loads((root / "failed/setup_status.json").read_text())
                self.assertEqual(status["state"], "failed")
                self.assertFalse((root / "failed/runtime_manifest.json").exists())
                with self.assertRaises(FileExistsError):
                    setup.prepare_runtime(root / "failed", analysis_utils_archive=archive)

    def test_environment_sanitized_without_home_override(self):
        inherited = {"HOME": "/unchanged-user-home", "PYTHONPATH": "/ambient/code",
                     "LD_LIBRARY_PATH": "/ambient/lib", "LD_PRELOAD": "bad.so",
                     "CASASITECONFIG": "/ambient/config", "PATH": "/bin"}
        clean = bootstrap.build_environment("/explicit/runtime", inherited)
        self.assertEqual(clean["HOME"], inherited["HOME"])
        self.assertEqual(clean["PATH"], inherited["PATH"])
        for name in ("PYTHONPATH", "LD_LIBRARY_PATH", "LD_PRELOAD", "CASASITECONFIG"):
            self.assertNotIn(name, clean)
        self.assertEqual(clean["CASARCFILES"], "/explicit/runtime/casarc")
        self.assertEqual(clean["WM_AGN_CASA_RUNTIME"], "/explicit/runtime")
        self.assertNotIn("HOME=", setup.wrapper_text("/python"))
        self.assertIn("unset PYTHONPATH", setup.wrapper_text("/python"))

    def test_exact_frozen_utils_selected_and_mutated_explicit_bundle_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            default = root / "analysis_utils/analysis_scripts"
            frozen = root / "frozen"
            for directory in (default, frozen):
                directory.mkdir(parents=True)
                for name in ("analysisUtils.py", "mpfit.py"):
                    (directory / name).write_text("# pinned\n")
            manifest = {"analysis_utils_files": bootstrap.bundle_inventory(default)}
            with mock.patch.dict(os.environ, {"WM_AGN_REQUESTED_UTILS_PATH": str(frozen)}, clear=True):
                self.assertEqual(bootstrap.select_analysis_utils(root, manifest), frozen)
            (frozen / "mpfit.py").write_text("# changed\n")
            with self.assertRaisesRegex(RuntimeError, "differs"):
                bootstrap.select_analysis_utils(root, manifest, explicit=frozen)

    def test_dispatch_preserves_command_script_and_module_arguments(self):
        saved_argv, saved_path, saved_cwd = list(sys.argv), list(sys.path), Path.cwd()
        try:
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                bootstrap.dispatch_target(["-c", "import json,sys; print(json.dumps(sys.argv))", "argument"])
            self.assertEqual(json.loads(output.getvalue()), ["-c", "argument"])
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                script = root / "probe.py"
                script.write_text("import json,sys; print(json.dumps(sys.argv))\n")
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    bootstrap.dispatch_target([str(script), "--flag", "value"])
                self.assertEqual(json.loads(output.getvalue()), [str(script), "--flag", "value"])
                os.chdir(root)
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    bootstrap.dispatch_target(["-m", "probe", "module-argument"])
                arguments = json.loads(output.getvalue())
                self.assertEqual(Path(arguments[0]).name, "probe.py")
                self.assertEqual(arguments[1:], ["module-argument"])
            with self.assertRaises(ValueError):
                bootstrap.dispatch_target(["-c"])
        finally:
            sys.argv[:] = saved_argv
            sys.path[:] = saved_path
            os.chdir(saved_cwd)


if __name__ == "__main__":
    unittest.main()
