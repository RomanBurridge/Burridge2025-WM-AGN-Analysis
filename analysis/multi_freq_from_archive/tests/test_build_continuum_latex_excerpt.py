from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import build_continuum_latex_excerpt as excerpt  # noqa: E402


EXPECTED_ACCEPTED_TABLE = r"""\begin{deluxetable*}{lllllll}{
\tablecaption {Continuum Flux-density Measurements}
\label{tab:mr1}
\tablehead{ 
\colhead{[1]} & \colhead{[2]} & \colhead{[3]}  & \colhead{[4]} & \colhead{[5]} & \colhead{[6]}  & \colhead{[7]} }
}
\startdata
Circinus & ALMA & 39 & 67.8 & 2024/11/24\_1 & 7.52 & 593 \\
Circinus & ALMA & 46 & 48.2 & 2025/01/04\_1 & 3.26 & 116 \\
Circinus & ALMA & 86 & 22.8 & 2016/07/02\_1 & 1.38 & 98 \\
Circinus & ALMA & 93 & 25.7 & 2021/08/01\_1 & 0.219 & 290 \\
Circinus & ALMA & 93 & 42.3 & 2016/01/14\_1 & 2.64 & 143 \\
Circinus & ALMA & 98 & 35.8 & 2024/11/27\_1 & 2.94 & 116 \\
Circinus & ALMA & 107 & 36.0 & 2025/01/26\_1 & 1.35 & 523 \\
Circinus & ALMA & 108 & 40.6 & 2014/12/27\_1 & 2.40 & 158 \\
Circinus & ALMA & 108 & 62.4 & 2014/07/02\_1 & 12.9 & 54 \\
Circinus & ALMA & 108 & 39.8 & 2016/03/18\_1 & 1.54 & 124 \\
\enddata
\tablecomments{
[1]: Source Name
[2]: Telescope
[3]: Frequency [GHz]
[4]: Flux Density [mJy]
[5]: Date [yyyy/mm/dd\_x]
[6]: Beam Size FWHM [arcsec]
[7]: SNR
\\
Ten representative rows are shown; the full 501-row table is provided as \texttt{tab2.txt}. It comprises 485 image-based measurements and 16 literature measurements. Image rows use the beam-shaped fits described in Section~\ref{subsec:casa}; literature rows preserve the reported quantities. An asterisk marks an upper limit, and a date suffix distinguishes additional products from the same day. Full numerical precision remains in the machine-readable table.}
\end{deluxetable*}
"""


class ContinuumLatexExcerptTests(unittest.TestCase):
    def test_half_up_display_rounding_is_explicit(self) -> None:
        self.assertEqual(excerpt.format_decimal_places(42.25, 1), "42.3")
        self.assertEqual(excerpt.format_decimal_places(40.55, 1), "40.6")
        self.assertEqual(excerpt.format_significant_figures(0.21017, 3), "0.210")

    def test_excerpt_is_selected_from_current_canonical_rows(self) -> None:
        self.assertEqual(excerpt.EXPECTED_SOURCE_ROWS, 53)
        rows = excerpt.excerpt_rows()
        self.assertEqual(len(rows), 10)
        self.assertTrue(all(row["Source"] == "Circinus" for row in rows))
        target = [row for row in rows if row["date"] == "2024/11/24_1"]
        self.assertEqual(len(target), 1)
        self.assertEqual(target[0]["snr"], 592.681)

    def test_rendered_table_exactly_matches_accepted_r3_excerpt(self) -> None:
        self.assertEqual(excerpt.render_table(), EXPECTED_ACCEPTED_TABLE)
        self.assertNotIn("Provenance", EXPECTED_ACCEPTED_TABLE)
        data_rows = [
            line
            for line in EXPECTED_ACCEPTED_TABLE.splitlines()
            if line.startswith("Circinus &")
        ]
        self.assertEqual(len(data_rows), 10)
        self.assertTrue(all(line.count("&") == 6 for line in data_rows))

    def test_cli_build_and_check(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "Continuum_from_Pipeline.tex"
            self.assertEqual(excerpt.main(["--output", str(output)]), 0)
            self.assertEqual(excerpt.main(["--output", str(output), "--check"]), 0)


if __name__ == "__main__":
    unittest.main()
