from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from continuum_output_paths import fit_output_directory


class OutputPathTests(unittest.TestCase):
    def test_archive_folder_spelling_is_retained(self):
        self.assertEqual(fit_output_directory('ALMA', 'Circinus',
                         'ALMA/CircinusGalaxy/example.fits'),
                         'ALMA/CircinusGalaxy/fitfiles')

    def test_sma_lowercase_suffix_has_own_canonical_folder(self):
        self.assertEqual(fit_output_directory('SMA', 'NGC5765B',
                         'analysis/sma_calibration/fields/NGC5765b_test.image'),
                         'SMA/NGC5765B/fitfiles')

    def test_no_matching_other_source_can_be_inherited(self):
        with self.assertRaises(ValueError):
            fit_output_directory('VLA', 'NGC4258', 'VLA/NGC3079/a.fits')

    def test_prefix_is_not_source_identity(self):
        with self.assertRaises(ValueError):
            fit_output_directory('SMA', 'NGC425', 'NGC4258_test.image')

    def test_wrong_telescope_rejected(self):
        with self.assertRaises(ValueError):
            fit_output_directory('ALMA', 'NGC4258', 'VLA/NGC4258/a.fits')

    def test_traversal_rejected(self):
        with self.assertRaises(ValueError):
            fit_output_directory('VLA', 'NGC4258', 'VLA/NGC4258/../a.fits')

    def test_two_sources_never_share_model_directory(self):
        first = fit_output_directory('SMA', 'NGC4258', 'NGC4258_test.image')
        second = fit_output_directory('SMA', 'NGC5765B', 'NGC5765b_test.image')
        self.assertNotEqual(first, second)


if __name__ == '__main__':
    unittest.main()
