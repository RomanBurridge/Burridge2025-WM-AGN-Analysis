from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_DIR))

from continuum_units import (  # noqa: E402
    brightness_to_mjy_per_beam,
    brightness_to_mjy_per_beam_factor,
    flux_density_to_mjy,
)


class ContinuumUnitTests(unittest.TestCase):
    def test_active_fitter_uses_strict_unit_helpers(self) -> None:
        source = (PACKAGE_DIR / "get_continuum_params.py").read_text()
        self.assertIn("flux_density_to_mjy", source)
        self.assertIn("image_rms_mjy_per_beam", source)
        self.assertNotIn("['rms'][0]*1000", source)
        self.assertNotIn("data_main *= 1000", source)

    def test_component_flux_units_are_converted_to_mjy(self) -> None:
        self.assertEqual(flux_density_to_mjy(1.25, "Jy"), 1250.0)
        self.assertEqual(flux_density_to_mjy(1.25, "mJy"), 1.25)
        self.assertEqual(flux_density_to_mjy(1250, "uJy"), 1.25)
        self.assertTrue(math.isnan(flux_density_to_mjy(float("nan"), "Jy")))

    def test_image_brightness_units_are_converted_to_mjy_per_beam(self) -> None:
        self.assertEqual(brightness_to_mjy_per_beam(0.002, "Jy/beam"), 2.0)
        self.assertEqual(brightness_to_mjy_per_beam(2.0, "mJy beam-1"), 2.0)
        self.assertEqual(brightness_to_mjy_per_beam(2000, "uJy/beam"), 2.0)
        self.assertEqual(brightness_to_mjy_per_beam_factor("Jy/beam"), 1000.0)

    def test_incompatible_or_missing_units_fail_explicitly(self) -> None:
        for unit in ("", "K", "Jy/pixel"):
            with self.assertRaises(ValueError, msg=unit):
                brightness_to_mjy_per_beam(1.0, unit)
        with self.assertRaises(ValueError):
            flux_density_to_mjy(1.0, "Jy/beam")


if __name__ == "__main__":
    unittest.main()
