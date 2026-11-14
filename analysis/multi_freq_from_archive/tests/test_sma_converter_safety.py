from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[2] / "sma_calibration/sma_data/sma2casa.py"
spec = importlib.util.spec_from_file_location("sma_converter", SCRIPT)
converter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(converter)


class FakeUVData:
    def __init__(self):
        self.reads, self.writes = [], []

    def read_mir(self, path):
        self.reads.append(path)

    def write_ms(self, path, *, clobber):
        assert clobber is False
        self.writes.append(path)
        output = Path(path)
        output.mkdir(exist_ok=False)
        (output / "table.info").write_text("Type = Measurement Set\n")


class ConverterSafetyTests(unittest.TestCase):
    def fixture(self, root):
        folder = root / "input track"
        folder.mkdir()
        (folder / "sch_read").write_bytes(b"protected input")
        return folder, root / "track_master.ms"

    def test_absolute_input_and_output_resolve_without_mangling(self):
        with tempfile.TemporaryDirectory() as temp:
            folder, output = self.fixture(Path(temp))
            fake = FakeUVData()
            result = converter.process_folder(folder, output_ms=output, uvdata_factory=lambda: fake)
            self.assertEqual(fake.reads, [str(folder.resolve())])
            self.assertEqual(fake.writes, [str(output.resolve())])
            self.assertEqual(result, output.resolve())
            self.assertEqual((folder / "sch_read").read_bytes(), b"protected input")

    def test_existing_output_refused_before_any_read(self):
        with tempfile.TemporaryDirectory() as temp:
            folder, output = self.fixture(Path(temp))
            output.mkdir()
            (output / "user_data").write_bytes(b"keep")
            fake = FakeUVData()
            with self.assertRaises(FileExistsError):
                converter.process_folder(folder, output_ms=output, uvdata_factory=lambda: fake)
            self.assertFalse(fake.reads)
            self.assertEqual((output / "user_data").read_bytes(), b"keep")

    def test_preexisting_uvfits_is_untouched(self):
        with tempfile.TemporaryDirectory() as temp:
            folder, output = self.fixture(Path(temp))
            uvfits = output.with_suffix(".uvfits")
            uvfits.write_bytes(b"preexisting UVFITS must survive")
            converter.process_folder(folder, output_ms=output, uvdata_factory=FakeUVData)
            self.assertEqual(uvfits.read_bytes(), b"preexisting UVFITS must survive")

    def test_partial_failure_preserved_and_retry_refused(self):
        class FailingUVData(FakeUVData):
            def write_ms(self, path, *, clobber):
                Path(path).mkdir()
                (Path(path) / "partial").write_bytes(b"diagnostic evidence")
                raise RuntimeError("test conversion failed")
        with tempfile.TemporaryDirectory() as temp:
            folder, output = self.fixture(Path(temp))
            with self.assertRaisesRegex(RuntimeError, "test conversion failed"):
                converter.process_folder(folder, output_ms=output, uvdata_factory=FailingUVData)
            self.assertEqual((output / "partial").read_bytes(), b"diagnostic evidence")
            with self.assertRaises(FileExistsError):
                converter.process_folder(folder, output_ms=output, uvdata_factory=FakeUVData)

    def test_output_inside_input_is_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            folder, _ = self.fixture(Path(temp))
            with self.assertRaisesRegex(ValueError, "inside"):
                converter.process_folder(folder, output_ms=folder / "new.ms", uvdata_factory=FakeUVData)

    def test_broken_output_symlink_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            folder, output = self.fixture(Path(temp))
            output.symlink_to(Path(temp) / "absent.ms")
            with self.assertRaises(FileExistsError):
                converter.process_folder(folder, output_ms=output, uvdata_factory=FakeUVData)
            self.assertTrue(output.is_symlink())

    def test_file_input_and_non_ms_output_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            folder, _ = self.fixture(root)
            with self.assertRaisesRegex(ValueError, "directory"):
                converter.process_folder(folder / "sch_read", output_ms=root / "out.ms", uvdata_factory=FakeUVData)
            with self.assertRaisesRegex(ValueError, "end in .ms"):
                converter.process_folder(folder, output_ms=root / "out.fits", uvdata_factory=FakeUVData)


if __name__ == "__main__":
    unittest.main()
