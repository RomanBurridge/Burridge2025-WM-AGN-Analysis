from __future__ import annotations

import csv
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest import mock

try:
    import numpy as np
    from astropy.io import fits
except ImportError:
    np = fits = None


PACKAGE = Path(__file__).resolve().parents[1]
BUILDER = PACKAGE / "archive_inventory_audit/task28/build_approved_alma_inputs.py"


@unittest.skipIf(fits is None, "Astropy and NumPy are required")
class StokesISerializationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("task28_stokes_builder", BUILDER)
        cls.builder = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.builder)

    def make_parent(self, path: Path, *, beam_extension=False):
        data = np.arange(4 * 1 * 5 * 6, dtype=np.float32).reshape((4, 1, 5, 6))
        header = fits.Header()
        for axis, ctype, crval, delta in (
            (1, "RA---SIN", 10.0, -0.0001),
            (2, "DEC--SIN", -20.0, 0.0001),
            (3, "FREQ", 230e9, 1e9),
            (4, "STOKES", 1.0, 1.0),
        ):
            header[f"CTYPE{axis}"] = ctype
            header[f"CRVAL{axis}"] = crval
            header[f"CDELT{axis}"] = delta
            header[f"CRPIX{axis}"] = 1.0
        header["DATE-OBS"] = "2020-01-02T03:04:05"
        header["TELESCOP"] = "ALMA"
        header["BUNIT"] = "Jy/beam"
        hdus = [fits.PrimaryHDU(data=data, header=header)]
        if beam_extension:
            columns = [
                fits.Column(name="BMAJ", format="D", unit="arcsec", array=[1., 2., 3., 4.]),
                fits.Column(name="BMIN", format="D", unit="arcsec", array=[.5, 1., 1.5, 2.]),
                fits.Column(name="BPA", format="D", unit="deg", array=[10., 20., 30., 40.]),
                fits.Column(name="CHAN", format="J", array=[0, 0, 0, 0]),
                fits.Column(name="POL", format="J", array=[0, 1, 2, 3]),
            ]
            hdus.append(fits.BinTableHDU.from_columns(columns, name="BEAMS"))
        else:
            hdus[0].header["BMAJ"] = 1 / 3600
            hdus[0].header["BMIN"] = .5 / 3600
            hdus[0].header["BPA"] = 10.
        fits.HDUList(hdus).writeto(path)
        return data

    def test_repeated_extraction_uses_no_wall_clock_and_valid_checksums(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            parent, left, right = (root / name for name in ("parent.fits", "a.fits", "b.fits"))
            expected = self.make_parent(parent)
            original = parent.read_bytes()
            with mock.patch.object(fits.PrimaryHDU, "_get_timestamp", side_effect=AssertionError("wall clock used")):
                self.builder.extract_stokes_i(parent, left)
                self.builder.extract_stokes_i(parent, right)
            self.assertEqual(left.read_bytes(), right.read_bytes())
            self.assertEqual(parent.read_bytes(), original)
            with fits.open(left, checksum=True) as hdul:
                np.testing.assert_array_equal(hdul[0].data, expected[:1])
                self.assertEqual(hdul[0].verify_checksum(), 1)
                self.assertEqual(hdul[0].verify_datasum(), 1)
                self.assertEqual(hdul[0].header["CRVAL3"], 230e9)
                self.assertAlmostEqual(hdul[0].header["BMAJ"], 1 / 3600)

    def test_exact_reviewed_parent_identities_cover_five_only(self):
        with (BUILDER.parent / "approved_alma_products_validated.csv").open(newline="") as handle:
            expected = {row["parent_sha256"] for row in csv.DictReader(handle)
                        if row["transformation"] == "extract_singleton_stokes_i_from_iquv"}
        self.assertEqual(set(self.builder.REVIEWED_STOKES_I_CHECKSUM_TIMES), expected)
        self.assertEqual(len(expected), 5)

    def test_changed_pixels_or_wcs_change_whole_file_hash(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            parent = root / "parent.fits"
            self.make_parent(parent)
            self.builder.extract_stokes_i(parent, root / "original.fits")
            expected = hashlib.sha256((root / "original.fits").read_bytes()).hexdigest()
            for label in ("pixel", "wcs", "beam"):
                changed = root / f"{label}_parent.fits"
                with fits.open(parent) as hdul:
                    if label == "pixel":
                        hdul[0].data[0, 0, 0, 0] += 1
                    elif label == "wcs":
                        hdul[0].header["CRVAL1"] += .001
                    else:
                        hdul[0].header["BMAJ"] *= 2
                    hdul.writeto(changed)
                result = root / f"{label}.fits"
                self.builder.extract_stokes_i(changed, result)
                self.assertNotEqual(hashlib.sha256(result.read_bytes()).hexdigest(), expected)

    def test_preserves_beam_for_stokes_i_from_beams_extension(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            parent, output = root / "parent.fits", root / "i.fits"
            self.make_parent(parent, beam_extension=True)
            self.builder.extract_stokes_i(parent, output)
            header = fits.getheader(output)
            self.assertAlmostEqual(header["BMAJ"], 1 / 3600)
            self.assertAlmostEqual(header["BMIN"], .5 / 3600)
            self.assertEqual(header["BPA"], 10.)

    def test_existing_output_is_not_modified(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            parent, output = root / "parent.fits", root / "i.fits"
            self.make_parent(parent)
            output.write_bytes(b"preexisting user data")
            with self.assertRaises(FileExistsError):
                self.builder.extract_stokes_i(parent, output)
            self.assertEqual(output.read_bytes(), b"preexisting user data")

    def test_no_admission_checksum_waiver(self):
        source = (PACKAGE / "continuum_product_additions.py").read_text()
        self.assertIn('if observed != row["product_sha256"]:', source)
        self.assertIn('admitted product SHA-256 mismatch', source)


if __name__ == "__main__":
    unittest.main()
