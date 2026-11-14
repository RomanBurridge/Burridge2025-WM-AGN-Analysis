from __future__ import annotations

import ast
import copy
import importlib.util
from pathlib import Path
import tempfile
import unittest

CALIBRATION = Path(__file__).resolve().parents[2] / "sma_calibration"
spec = importlib.util.spec_from_file_location("sma_channelization", CALIBRATION / "sma_channelization.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def records():
    result = []
    for spw in range(12):
        sign = 1 if spw % 2 == 0 else -1
        result.append({"NUM_CHAN": 128, "TOTAL_BANDWIDTH": 2288000000.,
                       "CHAN_FREQ": [230e9 + sign * i * 17875000. for i in range(128)],
                       "CHAN_WIDTH": [17875000.] * 128,
                       "EFFECTIVE_BW": [17875000.] * 128,
                       "RESOLUTION": [17875000.] * 128})
    return result


class ChannelizationTests(unittest.TestCase):
    def test_both_released_tracks_and_opposite_sideband_order_pass(self):
        for date in (220131, 230316):
            result = module.validate_spectral_windows(records(), date=date)
            self.assertTrue(result["metadata_valid"])
            self.assertEqual(len(result["spectral_windows"]), 12)
            self.assertEqual(result["spectral_windows"][0]["channel_width_hz"], 17875000.)

    def test_second_averaging_rejected(self):
        altered = records()
        altered[0]["NUM_CHAN"] = 64
        with self.assertRaisesRegex(ValueError, "128 already-binned"):
            module.validate_spectral_windows(altered, date=230316)

    def test_mismatched_spacing_width_bandwidth_or_resolution_rejected(self):
        for field in ("CHAN_WIDTH", "TOTAL_BANDWIDTH", "CHAN_FREQ", "EFFECTIVE_BW", "RESOLUTION"):
            altered = records()
            if field == "TOTAL_BANDWIDTH":
                altered[0][field] *= 2
            else:
                altered[0][field][10] *= 2
            with self.subTest(field=field), self.assertRaises(ValueError):
                module.validate_spectral_windows(altered, date=230316)

    def test_nonfinite_and_short_column_rejected(self):
        for value in (float("nan"), float("inf"), 0.):
            altered = records()
            altered[0]["EFFECTIVE_BW"][0] = value
            with self.assertRaises(ValueError):
                module.validate_spectral_windows(altered, date=230316)
        altered = records()
        altered[0]["CHAN_FREQ"].pop()
        with self.assertRaises(ValueError):
            module.validate_spectral_windows(altered, date=230316)

    def test_unknown_track_or_invalid_averaging_rejected(self):
        with self.assertRaisesRegex(ValueError, "contract"):
            module.validate_spectral_windows(records(), date=20200101)
        for averaging in (True, 0, -1, 3, 128.0):
            with self.assertRaises(ValueError):
                module.validate_spectral_windows(records(), date=230316, averaging=averaging)

    def test_read_only_table_closed_before_validation_failure(self):
        class Table:
            closed = False
            def open(self, path, *, nomodify):
                self.nomodify = nomodify
            def nrows(self):
                return 12
            def getcell(self, field, row):
                return 64 if field == "NUM_CHAN" else records()[row][field]
            def close(self):
                self.closed = True
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "SPECTRAL_WINDOW").mkdir()
            table = Table()
            with self.assertRaises(ValueError):
                module.check_measurement_set(root, date=230316, table_factory=lambda: table)
            self.assertTrue(table.nomodify)
            self.assertTrue(table.closed)

    def test_calibrators_preflight_before_clearcal(self):
        for name in ("calibrate_sma_clean.py", "calibrate_sma_dirty.py"):
            tree = ast.parse((CALIBRATION / name).read_text())
            calls = [(node.func.id, node.lineno) for node in ast.walk(tree)
                     if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)]
            preflight = [line for function, line in calls if function == "check_measurement_set"]
            clearcal = [line for function, line in calls if function == "clearcal"]
            self.assertEqual(len(preflight), 1)
            self.assertLess(preflight[0], min(clearcal))


if __name__ == "__main__":
    unittest.main()
