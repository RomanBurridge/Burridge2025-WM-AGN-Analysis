from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from scripts import test_release


RELEASE = Path(__file__).resolve().parents[3]


class ReleaseTestOrchestratorTests(unittest.TestCase):
    def test_inventory_routes_every_module_exactly_once(self) -> None:
        fitting, renderer = test_release.discover_test_modules(RELEASE)
        expected = {
            f"{package}.{path.stem}"
            for package, relative in test_release.TEST_PACKAGES
            for path in (RELEASE / relative).glob("test_*.py")
        }
        self.assertEqual(set(fitting) | set(renderer), expected)
        self.assertFalse(set(fitting) & set(renderer))
        self.assertEqual(set(renderer), test_release.RENDERER_MODULES)
        self.assertIn(
            "analysis.multi_freq_from_archive.tests.test_release_orchestrator",
            fitting,
        )

    def test_cli_requires_both_python_environments_and_output(self) -> None:
        with self.assertRaises(SystemExit):
            test_release.parse_args([])

    def test_fitting_argv_starts_with_wrapper_supported_module_mode(self) -> None:
        python = Path("/validation/wm-agn-casa-python")
        command = test_release.fitting_test_command(
            python,
            ["analysis.example.tests.test_one"],
        )
        self.assertEqual(
            command,
            [
                str(python),
                "-m",
                "unittest",
                "-v",
                "analysis.example.tests.test_one",
            ],
        )
        self.assertNotIn("-B", command)

    def test_output_must_be_new_and_external(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_release = root / "release"
            (fake_release / "analysis").mkdir(parents=True)
            (fake_release / "scripts").mkdir()
            rootfs = root / "rootfs"
            (rootfs / "usr/bin").mkdir(parents=True)
            (rootfs / "usr/bin/python3").write_bytes(b"python")
            python = root / "python"
            python.write_bytes(b"python")
            python.chmod(0o700)
            existing = root / "existing"
            existing.mkdir()
            with self.assertRaises(FileExistsError):
                test_release._validate_paths(
                    fake_release, python, rootfs, existing
                )
            with self.assertRaises(ValueError):
                test_release._validate_paths(
                    fake_release, python, rootfs, fake_release / "test-output"
                )

    def test_source_snapshot_detects_stable_file_change(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.py"
            source.write_text("first\n", encoding="utf-8")
            before = test_release.source_snapshot(root)
            source.write_text("second\n", encoding="utf-8")
            after = test_release.source_snapshot(root)
            self.assertNotEqual(before, after)


if __name__ == "__main__":
    unittest.main()
