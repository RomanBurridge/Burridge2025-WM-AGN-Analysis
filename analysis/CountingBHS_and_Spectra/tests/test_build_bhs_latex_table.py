from __future__ import annotations

import re
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import build_bhs_latex_table as table  # noqa: E402


NUMBER_RE = re.compile(r"(?<![\d.])[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?![\d.])")
REFERENCE_PREFIX_RE = re.compile(r"^\$\^\{\d+\}\$\\,")
TABLE_NOTE_MARK_RE = re.compile(r"\\tablenotemark\{([a-d])\}$")
INCOMPLETE_SIZE_KEYS = frozenset(
    {"NGC1068", "NGC4388", "Circinus", "NGC2960", "IC2560", "J0437"}
)


def strip_reference_prefix(cell: str) -> str:
    return REFERENCE_PREFIX_RE.sub("", cell)


def strip_table_note_mark(cell: str) -> str:
    return TABLE_NOTE_MARK_RE.sub("", cell)


def clean_numeric_tokens(cell: str) -> list[str]:
    """Return displayed numeric magnitudes, excluding a cell reference mark."""

    return [
        match.lstrip("+-")
        for match in NUMBER_RE.findall(strip_reference_prefix(cell))
    ]


def expected_uncertainty_tokens(
    lower: float | None, upper: float | None
) -> list[str]:
    if lower is None and upper is None:
        return []
    if lower is None or upper is None:
        raise AssertionError("incomplete source interval")
    if table.values_close(lower, upper):
        return [table.format_significant_figures(lower, 2)]
    return [
        table.format_significant_figures(lower, 2),
        table.format_significant_figures(upper, 2),
    ]


class BuildBhsLatexTableTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.params = table.load_params()
        cls.metadata = table.load_metadata()
        table.validate_inputs(cls.params, cls.metadata)
        cls.rendered = table.render_table(cls.params, cls.metadata)
        cls.by_key = {row["source_key"]: row for row in cls.params}
        cls.meta_by_key = {row["source_key"]: row for row in cls.metadata}

        body = cls.rendered.split("\\startdata\n", 1)[1].split("\\enddata", 1)[0]
        all_rows = [line for line in body.splitlines() if line.strip()]
        cls.comparison_rows = all_rows[:2]
        cls.sample_rows: dict[str, list[str]] = {}
        for source_key, line in zip(
            table.DISPLAY_SOURCE_KEYS, all_rows[2:], strict=True
        ):
            cells = line.removesuffix(r" \\").split(" & ")
            if len(cells) != 6:
                raise AssertionError(f"malformed table row: {line}")
            cls.sample_rows[source_key] = cells

    def test_exact_inventory_and_deterministic_render(self) -> None:
        self.assertEqual(len(self.params), 21)
        self.assertEqual(len(self.metadata), 21)
        self.assertEqual(len(self.sample_rows), 21)
        self.assertEqual(
            self.rendered, table.render_table(table.load_params(), table.load_metadata())
        )

    def test_six_column_grid_uses_one_compact_comment_and_lettered_markers(self) -> None:
        grid_open = self.rendered.index("\\begin{longrotatetable}")
        data_open = self.rendered.index("\\begin{deluxetable*}{llllll}")
        data_close = self.rendered.index("\\enddata")
        grid_close = self.rendered.index("\\end{longrotatetable}")
        comments = self.rendered.index("\\tablecomments{")
        references = self.rendered.index("\\tablerefs{")

        self.assertLess(grid_open, data_open)
        self.assertLess(data_open, data_close)
        self.assertLess(data_close, grid_close)
        self.assertLess(grid_close, comments)
        self.assertLess(comments, references)
        self.assertEqual(self.rendered.count("\\begin{longrotatetable}"), 1)
        self.assertEqual(self.rendered.count("\\end{longrotatetable}"), 1)
        self.assertEqual(self.rendered.count("\\begin{deluxetable*}{llllll}"), 1)
        self.assertEqual(self.rendered.count("\\end{deluxetable*}"), 1)
        self.assertIn(r"\setlength{\tabcolsep}{2pt}", self.rendered[grid_open:data_open])
        self.assertEqual(self.rendered.count(r"\tablenotemark"), 12)
        self.assertNotIn(r"\tablenotetext", self.rendered)
        self.assertIn(
            r"\colhead{\shortstack{Systematic/model term\\{}"
            r"[$\mu\mathrm{as}$]}}",
            self.rendered,
        )
        self.assertIn(r"\colhead{$i$ [deg] / geometry}", self.rendered)
        self.assertIn(r"\tabletypesize{\footnotesize}", self.rendered[data_open:data_close])
        comment = self.rendered.split(r"\tablecomments{", 1)[1].split("}\n", 1)[0]
        expected_comment = (
            r"Boldface identifies the six sources discussed individually in "
            r"Section~\ref{sec:individual-sources}.  "
            r"Unmarked uncertainties accompanying $\theta_{\rm BHS}$ are "
            "statistical.  A dagger marks an uncertainty that includes other "
            "effects or was published without separating them.  Lettered entries "
            "in the systematic/model column denote (a) a separate allowance for "
            r"the fitted $M_\bullet/D$, and therefore the BHS size, to be 1\% "
            "higher or lower because of small imperfections in the Keplerian disk "
            "model, (b) a model sensitivity, (c) a position systematic, and (d) a "
            "geometry systematic; these terms remain separate from the uncertainty "
            r"accompanying $\theta_{\rm BHS}$.  A leading $+$ denotes a one-sided "
            "increase in the predicted BHS size.  Asymmetric intervals are displayed "
            "using the average of their upper and lower half-widths; exact bounds "
            "and incomplete components remain machine-readable."
        )
        self.assertEqual(comment, expected_comment)
        self.assertEqual(len(comment.split()), 121)
        self.assertLessEqual(len(comment.split()), 180)
        self.assertIn(
            r"these terms remain separate from the uncertainty accompanying $\theta_{\rm BHS}$",
            comment,
        )

        self.assertEqual(self.rendered.count(r"^{\dagger}"), 4)
        self.assertEqual(
            {
                key
                for key, cells in self.sample_rows.items()
                if r"^{\dagger}" in cells[1]
            },
            {"NGC1194", "NGC2273", "ESO558-G009"},
        )
        self.assertNotIn("not total", self.rendered.lower())
        additional = {
            "NGC4258": (r"$\pm 0.0054$", "a"),
            "NGC5765b": (r"$\pm 0.00038$", "a"),
            "CGCG074-064": (r"$\pm 0.00028$", "a"),
            "NGC6264": (r"$\pm 0.00021$", "a"),
            "UGC3789": (r"$\pm 0.00024$", "a"),
            "NGC6323": (r"$\pm 0.000095$", "a"),
            "NGC3393": (r"$\pm 0.0091$", "b"),
            "Circinus": (r"$+0.0012$", "b"),
            "NGC4945": (r"$+0.0093$", "b"),
            "NGC3079": (r"$+0.0040$", "b"),
            "UGC6093": (r"$\pm 0.00087$", "c"),
            "ESO558-G009": (r"$\pm 0.00097$", "d"),
        }
        for key, cells in self.sample_rows.items():
            if key in additional:
                term, marker = additional[key]
                self.assertEqual(
                    cells[2], rf"{term}\tablenotemark{{{marker}}}", key
                )
            else:
                self.assertEqual(cells[2], "---", key)
            if key in INCOMPLETE_SIZE_KEYS:
                self.assertNotIn(r"\pm", cells[1], key)
            self.assertNotIn(r"\shortstack", cells[2], key)
            self.assertNotIn(r"\ldots", cells[2], key)
            self.assertNotIn("—", cells[2], key)
            numeric_cell = strip_table_note_mark(cells[2])
            self.assertNotRegex(
                numeric_cell.replace(r"\pm", ""), r"[A-Za-z]", key
            )
        self.assertEqual(
            sum(cells[2] != "---" for cells in self.sample_rows.values()), 12
        )
        self.assertIn(
            "incomplete components remain machine-readable",
            comment,
        )
        self.assertIn("exact bounds and incomplete components remain machine-readable", comment)

    def test_every_sample_size_number_comes_from_parameter_csv(self) -> None:
        for param in self.params:
            key = param["source_key"]
            actual = clean_numeric_tokens(
                self.sample_rows[key][1] + " " + self.sample_rows[key][2]
            )
            central = table.optional_float(param["size_uas"])
            if central is None:
                self.assertEqual(actual, [], key)
                self.assertEqual(self.sample_rows[key][1], r"\nodata")
                self.assertEqual(self.sample_rows[key][2], "---")
                continue

            stat_lower, stat_upper = table.display_size_stat_bounds(
                param,
                table.optional_float(param["size_stat_unc_minus_uas"]),
                table.optional_float(param["size_stat_unc_plus_uas"]),
            )
            sys_lower = table.optional_float(param["size_systematic_unc_minus_uas"])
            sys_upper = table.optional_float(param["size_systematic_unc_plus_uas"])
            if sys_lower is not None or sys_upper is not None:
                self.assertIsNotNone(sys_lower, key)
                self.assertIsNotNone(sys_upper, key)
                sys_lower, sys_upper = table.symmetric_display_interval_bounds(
                    sys_lower, sys_upper
                )
            uncertainty_values = [
                value
                for value in (stat_lower, stat_upper, sys_lower, sys_upper)
                if value is not None
            ]
            expected_central = (
                table.format_at_decimal_place(
                    central, table.central_decimal_places(uncertainty_values)
                )
                if uncertainty_values
                else table.format_significant_figures(central, 3)
            )
            expected = [expected_central]
            if (
                param["size_stat_component_status"] == "no_formal_uncertainty"
                or key in INCOMPLETE_SIZE_KEYS
            ):
                if key in {"Circinus", "NGC4945", "NGC3079"}:
                    model_lower = table.optional_float(
                        param["size_model_check_lower_uas"]
                    )
                    model_upper = table.optional_float(
                        param["size_model_check_upper_uas"]
                    )
                    self.assertIsNotNone(model_lower, key)
                    self.assertIsNotNone(model_upper, key)
                    expected.append(
                        table.format_significant_figures(
                            model_upper - model_lower, 2
                        )
                    )
                else:
                    self.assertEqual(self.sample_rows[key][2], "---")
            else:
                expected += expected_uncertainty_tokens(stat_lower, stat_upper)
                expected += expected_uncertainty_tokens(sys_lower, sys_upper)
            self.assertEqual(actual, expected, key)

    def test_every_sample_mass_number_comes_from_parameter_csv(self) -> None:
        for param in self.params:
            key = param["source_key"]
            actual = clean_numeric_tokens(self.sample_rows[key][3])
            mass = table.optional_float(param["mass_msun"])
            if mass is None:
                self.assertEqual(actual, [], key)
                self.assertEqual(self.sample_rows[key][3], r"\nodata")
                continue

            central, transformed_lower, transformed_upper = table.mass_log_display_bounds(
                param
            )
            lower = table.optional_float(param["mass_display_unc_minus_msun"])
            upper = table.optional_float(param["mass_display_unc_plus_msun"])
            display_uncertainty = lower is not None and upper is not None
            if display_uncertainty:
                lower_components = [transformed_lower]
                upper_components = [transformed_upper]
                if (
                    param["mass_display_interval_status"]
                    == "separate_statistical_and_systematic_components"
                ):
                    systematic_lower = table.optional_float(
                        param["log_mass_systematic_unc_minus"]
                    )
                    systematic_upper = table.optional_float(
                        param["log_mass_systematic_unc_plus"]
                    )
                    self.assertIsNotNone(systematic_lower, key)
                    self.assertIsNotNone(systematic_upper, key)
                    lower_components.append(systematic_lower)
                    upper_components.append(systematic_upper)
                elif (
                    param["mass_display_interval_status"]
                    in {
                        "separate_decomposed_and_position_systematic_components",
                        "separate_decomposed_and_geometry_systematic_components",
                    }
                ):
                    systematic_lower = table.optional_float(
                        param["log_mass_systematic_unc_minus"]
                    )
                    systematic_upper = table.optional_float(
                        param["log_mass_systematic_unc_plus"]
                    )
                    self.assertIsNotNone(systematic_lower, key)
                    self.assertIsNotNone(systematic_upper, key)
                    lower_components.append(systematic_lower)
                    upper_components.append(systematic_upper)
                transformed_lower, transformed_upper = table.quadrature_interval_bounds(
                    lower_components, upper_components
                )
                expected = [
                    table.format_at_decimal_place(
                        central,
                        table.central_decimal_places(
                            [transformed_lower, transformed_upper]
                        ),
                    )
                ]
                expected += expected_uncertainty_tokens(
                    transformed_lower, transformed_upper
                )
            else:
                expected = [table.format_significant_figures(central, 3)]
            self.assertEqual(actual, expected, key)

    def test_every_sample_distance_number_comes_from_parameter_csv(self) -> None:
        for param in self.params:
            key = param["source_key"]
            actual = clean_numeric_tokens(self.sample_rows[key][4])
            central = table.optional_float(param["distance_mpc"])
            self.assertIsNotNone(central, key)
            stat_lower = table.optional_float(param["distance_stat_unc_minus_mpc"])
            stat_upper = table.optional_float(param["distance_stat_unc_plus_mpc"])
            sys_lower = table.optional_float(
                param["distance_systematic_unc_minus_mpc"]
            )
            sys_upper = table.optional_float(
                param["distance_systematic_unc_plus_mpc"]
            )
            lower_components = [value for value in (stat_lower, sys_lower) if value is not None]
            upper_components = [value for value in (stat_upper, sys_upper) if value is not None]
            uncertainty_values: list[float] = []
            if lower_components:
                lower, upper = table.quadrature_interval_bounds(
                    lower_components, upper_components
                )
                uncertainty_values = [lower, upper]
            expected_central = (
                table.format_at_decimal_place(
                    central, table.central_decimal_places(uncertainty_values)
                )
                if uncertainty_values
                else table.format_significant_figures(central, 3)
            )
            expected = [expected_central]
            if uncertainty_values:
                expected += expected_uncertainty_tokens(lower, upper)
            self.assertEqual(actual, expected, key)

    def test_publication_rounding_primitives_are_exact(self) -> None:
        self.assertEqual(table.format_significant_figures(0.000614185, 2), "0.00061")
        self.assertEqual(table.format_significant_figures(0.000204728, 2), "0.00020")
        self.assertEqual(table.format_significant_figures(0.8, 2), "0.80")
        self.assertEqual(table.format_significant_figures(3.7, 3), "3.70")
        self.assertEqual(table.format_at_decimal_place(132.1, 0), "132")
        self.assertEqual(
            table.central_decimal_places([0.004103, 0.010258, 0.008206]), 4
        )
        self.assertEqual(
            table.symmetric_display_interval_bounds(1.0, 1.2), (1.1, 1.1)
        )
        self.assertEqual(
            table.symmetric_display_interval_bounds(1.0, 1.2000001),
            (1.1, 1.1),
        )
        self.assertEqual(
            table.symmetric_display_interval_bounds(
                0.010257858260937705, 0.008206286608750165
            ),
            (0.0091, 0.0091),
        )
        self.assertEqual(
            table.symmetric_display_interval_bounds(17.3, 21.2), (19.0, 19.0)
        )
        self.assertEqual(
            table.symmetric_display_interval_bounds(23.4, 34.2), (28.5, 28.5)
        )
        for bounds in ((0.0, 1.0), (-1.0, 1.0), (1.0, float("inf"))):
            with self.subTest(bounds=bounds):
                with self.assertRaises(ValueError):
                    table.symmetric_display_interval_bounds(*bounds)

    def test_representative_cells_lock_publication_rounding(self) -> None:
        self.assertEqual(
            self.sample_rows["NGC4258"][1], r"$0.53860\pm 0.00059$"
        )
        self.assertEqual(
            strip_table_note_mark(self.sample_rows["NGC4258"][2]),
            r"$\pm 0.0054$",
        )
        self.assertEqual(self.sample_rows["NGC3393"][1], r"$0.0636\pm 0.0041$")
        self.assertEqual(
            self.by_key["NGC3393"]["size_systematic_unc_minus_uas"],
            "0.010257858260937705",
        )
        self.assertEqual(
            self.by_key["NGC3393"]["size_systematic_unc_plus_uas"],
            "0.008206286608750165",
        )
        self.assertEqual(
            strip_table_note_mark(self.sample_rows["NGC3393"][2]),
            r"$\pm 0.0091$",
        )
        self.assertEqual(
            self.sample_rows["CGCG074-064"][1], r"$0.02835\pm 0.00019$"
        )
        self.assertEqual(
            strip_table_note_mark(self.sample_rows["CGCG074-064"][2]),
            r"$\pm 0.00028$",
        )
        self.assertEqual(
            strip_reference_prefix(self.sample_rows["NGC6264"][4]),
            r"$132\pm 19$",
        )
        self.assertEqual(self.sample_rows["NGC4945"][1], r"$0.0388$")
        self.assertEqual(
            self.sample_rows["NGC4945"][2], r"$+0.0093$\tablenotemark{b}"
        )
        self.assertEqual(
            strip_reference_prefix(self.sample_rows["NGC4945"][3]), r"$6.15$"
        )
        self.assertEqual(
            strip_reference_prefix(self.sample_rows["NGC4945"][4]), r"$3.70$"
        )
        self.assertEqual(
            self.sample_rows["IC2560"][1], r"$0.00954$"
        )
        self.assertEqual(self.sample_rows["IC2560"][2], "---")
        self.assertEqual(
            self.sample_rows["NGC4388"][1], r"$0.0454$"
        )
        self.assertEqual(self.sample_rows["NGC4388"][2], "---")
        self.assertEqual(
            strip_reference_prefix(self.sample_rows["NGC4388"][3]),
            r"$6.924\pm 0.049$",
        )
        self.assertEqual(
            strip_reference_prefix(self.sample_rows["NGC4388"][4]),
            r"$19.0\pm 2.1$",
        )

    def test_cgcg_mass_components_and_full_row_are_exact_and_fail_closed(self) -> None:
        row = self.by_key["CGCG074-064"]
        cells = self.sample_rows["CGCG074-064"]
        self.assertEqual(cells[1], r"$0.02835\pm 0.00019$")
        self.assertEqual(cells[2], r"$\pm 0.00028$\tablenotemark{a}")
        self.assertIn(r"7.384\pm 0.039", cells[3])
        self.assertIn(r"87.6\pm 7.7", cells[4])
        self.assertEqual(cells[5], r"$^{15}$\,$90.8\pm1.3$ (at center)")
        self.assertEqual(row["mass_systematic_unc_minus_msun"], "500000")
        self.assertEqual(row["mass_systematic_unc_plus_msun"], "500000")
        self.assertEqual(
            row["log_mass_systematic_unc_minus"], "0.0090670199703273974"
        )
        self.assertEqual(
            row["log_mass_systematic_unc_plus"], "0.0088815872792344722"
        )
        self.assertEqual(row["mass_total_unc_minus_msun"], "2061552.8128088303")
        self.assertEqual(row["mass_total_unc_plus_msun"], "2256102.8345356956")
        self.assertEqual(row["mass_display_unc_minus_msun"], "2000000")
        self.assertEqual(row["mass_display_unc_plus_msun"], "2200000")
        self.assertEqual(
            row["mass_display_interval_status"],
            "separate_statistical_and_systematic_components",
        )
        self.assertEqual(
            row["mass_systematic_status"],
            "source_reported_model_selection_systematic",
        )
        self.assertEqual(
            row["nonadopted_fit_key"],
            "villaraos_2022_cgcg074064_gr_ratio_fit",
        )
        self.assertEqual(row["nonadopted_fit_ratio_msun_per_mpc"], "274900")
        self.assertEqual(row["nonadopted_fit_size_uas"], "0.02819885235931776")

        mutations = (
            ("mass_systematic_unc_minus_msun", "500001"),
            ("mass_systematic_unc_plus_msun", "500001"),
            ("log_mass_systematic_unc_minus", "0.0091"),
            ("log_mass_systematic_unc_plus", "0.0089"),
            ("mass_systematic_status", "not_available"),
            ("mass_total_unc_minus_msun", "2061553"),
            ("mass_total_unc_plus_msun", "2256103"),
            ("mass_display_unc_minus_msun", "2000001"),
            ("mass_display_unc_plus_msun", "2200001"),
            ("mass_display_interval_status", "decomposed_components"),
            ("nonadopted_fit_ratio_msun_per_mpc", "274901"),
            ("nonadopted_fit_basis", "none"),
            ("nonadopted_fit_reference_doi", "10.0000/wrong"),
        )
        for field, value in mutations:
            with self.subTest(field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "CGCG074-064" else item
                    for item in self.params
                ]
                with self.assertRaises(ValueError):
                    table.validate_inputs(params, self.metadata)

        metadata = [dict(item) for item in self.metadata]
        for field, value in (
            ("inclination_tex", r"(0):90.8\pm0.6"),
            ("inclination_reference_key", "Pesce2020b"),
            ("mass_reference_key", "Pesce2020b"),
            ("distance_reference_key", "Pesce2020b"),
        ):
            with self.subTest(metadata_field=field):
                candidate_metadata = [dict(item) for item in metadata]
                next(
                    item
                    for item in candidate_metadata
                    if item["source_key"] == "CGCG074-064"
                )[field] = value
                with self.assertRaisesRegex(ValueError, "CGCG074-064"):
                    table.validate_inputs(self.params, candidate_metadata)

    def test_ugc6093_position_mass_systematic_and_full_row_are_exact_and_fail_closed(
        self,
    ) -> None:
        row = self.by_key["UGC6093"]
        cells = self.sample_rows["UGC6093"]
        self.assertEqual(cells[1], r"$0.01745\pm 0.00068$")
        self.assertEqual(cells[2], r"$\pm 0.00087$\tablenotemark{c}")
        self.assertEqual(cells[5], r"$^{17}$\,$93.8\pm4.0$")
        self.assertEqual(
            strip_reference_prefix(self.sample_rows["UGC6093"][3]),
            r"$7.388\pm 0.034$",
        )
        self.assertEqual(
            clean_numeric_tokens(self.sample_rows["UGC6093"][3]),
            ["7.388", "0.034"],
        )

        expected_fields = {
            "mass_msun": "24443424.32138069",
            "mass_fit_unc_minus_msun": "852677.59260630317",
            "mass_fit_unc_plus_msun": "1042161.5020743705",
            "mass_shared_distance_unc_minus_msun": "1147364.9666493207",
            "mass_shared_distance_unc_plus_msun": "1147293.4548368326",
            "mass_systematic_unc_minus_msun": "1222171.2160690345",
            "mass_systematic_unc_plus_msun": "1222171.2160690345",
            "log_mass_systematic_unc_minus": "0.022276394711152211",
            "log_mass_systematic_unc_plus": "0.021189299069938092",
            "mass_display_unc_minus_msun": "1429512.3097151271",
            "mass_display_unc_plus_msun": "1549962.2149966571",
            "mass_total_unc_minus_msun": "1880746.6403039885",
            "mass_total_unc_plus_msun": "1973850.386758076",
            "distance_stat_unc_minus_mpc": "3.386942772210658",
            "distance_stat_unc_plus_mpc": "3.3861051591654241",
            "distance_unc_minus_mpc": "6.746328117856665",
            "distance_unc_plus_mpc": "6.7459076394863375",
            "size_uas": "0.017445797174172238",
            "size_stat_unc_minus_uas": "0.00060857432002926412",
            "size_stat_unc_plus_uas": "0.00074381305781354504",
            "size_systematic_unc_minus_uas": "0.00087228985870861197",
            "size_systematic_unc_plus_uas": "0.00087228985870861197",
            "size_total_unc_minus_uas": "0.0010636034508241178",
            "size_total_unc_plus_uas": "0.0011463627098697107",
            "mass_display_interval_status": (
                "separate_decomposed_and_position_systematic_components"
            ),
            "mass_systematic_status": (
                "source_reported_absolute_position_systematic"
            ),
        }
        for field, expected in expected_fields.items():
            self.assertEqual(row[field], expected, field)

        mutations = (
            ("mass_msun", "24443425"),
            ("mass_systematic_unc_minus_msun", "1222172"),
            ("mass_systematic_unc_plus_msun", "1222172"),
            ("log_mass_systematic_unc_minus", "0.0223"),
            ("log_mass_systematic_unc_plus", "0.0212"),
            ("mass_total_unc_minus_msun", "1880747"),
            ("mass_total_unc_plus_msun", "1973851"),
            ("mass_display_unc_minus_msun", "1429513"),
            ("mass_display_unc_plus_msun", "1549963"),
            ("mass_systematic_status", "not_available"),
            ("mass_display_interval_status", "decomposed_components"),
        )
        for field, value in mutations:
            with self.subTest(field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "UGC6093" else item
                    for item in self.params
                ]
                with self.assertRaisesRegex(ValueError, "UGC6093"):
                    table.validate_inputs(params, self.metadata)

    def test_eso558_geometry_mass_full_row_and_metadata_are_exact_and_fail_closed(
        self,
    ) -> None:
        row = self.by_key["ESO558-G009"]
        cells = self.sample_rows["ESO558-G009"]
        self.assertEqual(cells[1], r"$0.01621\pm 0.00037^{\dagger}$")
        self.assertEqual(cells[2], r"$\pm 0.00097$\tablenotemark{d}")
        self.assertEqual(cells[5], r"$^{18}$\,$98\pm1$ (geometry estimate)")
        self.assertEqual(
            strip_reference_prefix(self.sample_rows["ESO558-G009"][3]),
            r"$7.225\pm 0.036$",
        )
        self.assertEqual(
            clean_numeric_tokens(self.sample_rows["ESO558-G009"][3]),
            ["7.225", "0.036"],
        )

        expected_fields = {
            "mass_msun": "16792963.464140728",
            "mass_fit_unc_minus_msun": "386238.15967523673",
            "mass_fit_unc_plus_msun": "386238.15967523673",
            "mass_shared_distance_unc_minus_msun": "868182.04846912157",
            "mass_shared_distance_unc_plus_msun": "868182.04846912157",
            "mass_systematic_unc_minus_msun": "1008137.4178931694",
            "mass_systematic_unc_plus_msun": "1008137.4178931694",
            "log_mass_systematic_unc_minus": "0.026887542911212189",
            "log_mass_systematic_unc_plus": "0.02531951831781824",
            "mass_display_unc_minus_msun": "950221.0191704632",
            "mass_display_unc_plus_msun": "950221.0191704632",
            "mass_total_unc_minus_msun": "1385373.9706771097",
            "mass_total_unc_plus_msun": "1385373.9706771097",
            "distance_stat_unc_minus_mpc": "3.4026220986146067",
            "distance_stat_unc_plus_mpc": "3.4026220986146067",
            "distance_unc_minus_mpc": "5.4950816714869104",
            "distance_unc_plus_mpc": "5.4950816714869104",
            "size_uas": "0.016206653386239871",
            "size_stat_unc_minus_uas": "0.00037275302788351703",
            "size_stat_unc_plus_uas": "0.00037275302788351703",
            "size_systematic_unc_minus_uas": "0.00097293927497563788",
            "size_systematic_unc_plus_uas": "0.00097293927497563788",
            "size_total_unc_minus_uas": "0.0010419000204369179",
            "size_total_unc_plus_uas": "0.0010419000204369179",
            "mass_display_interval_status": (
                "separate_decomposed_and_geometry_systematic_components"
            ),
            "mass_systematic_status": (
                "source_reported_midline_inclination_systematic"
            ),
            "nonadopted_fit_key": (
                "villaraos_2022_eso558g009_gr_ratio_fit"
            ),
            "nonadopted_fit_ratio_msun_per_mpc": "159400",
            "nonadopted_fit_ratio_unc_minus_msun_per_mpc": "2600",
            "nonadopted_fit_ratio_unc_plus_msun_per_mpc": "2500",
            "nonadopted_fit_size_uas": "0.016351026067934704",
            "nonadopted_fit_size_unc_minus_uas": "0.00026670431478438035",
            "nonadopted_fit_size_unc_plus_uas": "0.00025644645652344266",
            "nonadopted_fit_interval_status": (
                "bayesian_formal_conditional_on_model"
            ),
            "nonadopted_fit_adoption_status": (
                "reported_separately_not_adopted"
            ),
            "nonadopted_fit_basis": (
                "same_data_reanalysis_conditional_on_fixed_gao2017_"
                "inclination_and_selected_scatter_angle"
            ),
            "nonadopted_fit_reference_key": "Villaraos2022",
            "nonadopted_fit_reference_doi": "10.1093/mnras/stac2973",
        }
        for field, expected in expected_fields.items():
            self.assertEqual(row[field], expected, field)

        metadata = self.meta_by_key["ESO558-G009"]
        self.assertEqual(metadata["source_tex"], "ESO~558--G009")
        self.assertEqual(metadata["inclination_tex"], r"98\pm1")
        self.assertEqual(metadata["mass_reference_key"], "Gao2016b")
        self.assertEqual(metadata["distance_reference_key"], "Gao2016b")
        self.assertEqual(metadata["inclination_reference_key"], "Gao2016b")
        self.assertEqual(
            cells[2],
            r"$\pm 0.00097$\tablenotemark{d}",
        )
        self.assertEqual(
            cells[5],
            r"$^{18}$\,$98\pm1$ (geometry estimate)",
        )

        row_mutations = (
            ("mass_msun", "16792964"),
            ("mass_systematic_unc_minus_msun", "1008138"),
            ("mass_systematic_unc_plus_msun", "1008138"),
            ("log_mass_systematic_unc_minus", "0.0269"),
            ("log_mass_systematic_unc_plus", "0.0253"),
            ("mass_total_unc_minus_msun", "1385375"),
            ("mass_total_unc_plus_msun", "1385375"),
            ("mass_display_unc_minus_msun", "950222"),
            ("mass_display_unc_plus_msun", "950222"),
            ("mass_systematic_status", "not_available"),
            ("mass_display_interval_status", "decomposed_components"),
            ("nonadopted_fit_ratio_msun_per_mpc", "159401"),
            ("nonadopted_fit_ratio_unc_minus_msun_per_mpc", "2599"),
            ("nonadopted_fit_ratio_unc_plus_msun_per_mpc", "2501"),
            ("nonadopted_fit_basis", "none"),
            ("nonadopted_fit_reference_doi", "10.0000/wrong"),
            ("mass_msun", "NA"),
            ("size_uas", "NA"),
        )
        for field, value in row_mutations:
            with self.subTest(row_field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "ESO558-G009" else item
                    for item in self.params
                ]
                with self.assertRaisesRegex(ValueError, "ESO558-G009"):
                    table.validate_inputs(params, self.metadata)

        metadata_mutations = (
            ("inclination_tex", r"\mathrm{e.o.}"),
            ("inclination_reference_key", "Villaraos2022"),
            ("mass_reference_key", "Villaraos2022"),
            ("distance_reference_key", "Villaraos2022"),
        )
        for field, value in metadata_mutations:
            with self.subTest(metadata_field=field):
                candidate_metadata = [dict(item) for item in self.metadata]
                next(
                    item
                    for item in candidate_metadata
                    if item["source_key"] == "ESO558-G009"
                )[field] = value
                with self.assertRaisesRegex(ValueError, "ESO558-G009"):
                    table.validate_inputs(self.params, candidate_metadata)

    def test_source_and_inclination_cells_come_from_metadata(self) -> None:
        metadata_by_key = {row["source_key"]: row for row in self.metadata}
        display_metadata = [
            metadata_by_key[source_key] for source_key in table.DISPLAY_SOURCE_KEYS
        ]
        refs = {
            key: index
            for index, key in enumerate(table.reference_order(display_metadata), start=1)
        }
        expected_inclinations = {
            "NGC4258": r"$^{7}$\,$87.05\pm0.09$ (at $r=6.1\,\mathrm{mas}$)",
            "NGC1194": r"$^{8}$\,nearly edge-on",
            "NGC1068": r"$^{9}$\,$\sim90$",
            "NGC3393": r"$^{11}$\,nearly edge-on",
            "NGC4388": r"$^{8}$\,unconstrained",
            "Circinus": r"$^{12}$\,warped disk + outflow",
            "NGC5765b": r"$^{14}$\,$72.4\pm0.5$ (at center)",
            "NGC4945": r"$^{13}$\,nearly edge-on",
            "NGC2273": r"$^{8}$\,nearly edge-on",
            "CGCG074-064": r"$^{15}$\,$90.8\pm1.3$ (at center)",
            "NGC6264": r"$^{14}$\,$91.3\pm2.3$ (at center)",
            "UGC3789": r"$^{14}$\,$84.9\pm0.6$ (at center)",
            "UGC6093": r"$^{17}$\,$93.8\pm4.0$",
            "ESO558-G009": r"$^{18}$\,$98\pm1$ (geometry estimate)",
            "NGC2960": r"$^{8}$\,nearly edge-on",
            "IC2560": r"$^{21}$\,nearly edge-on",
            "NGC3079": r"$^{19}$\,thick, self-gravitating disk",
            "NGC6323": r"$^{14}$\,$91.5\pm0.3$ (at center)",
            "J0437": r"$^{22}$\,unconstrained",
            "Mrk1029": r"$^{18}$\,$79\pm2$ (geometry estimate)",
            "NGC1320": r"$^{18}$\,unconstrained",
        }
        for meta in self.metadata:
            cells = self.sample_rows[meta["source_key"]]
            self.assertEqual(cells[5], expected_inclinations[meta["source_key"]])
            self.assertEqual(cells[5], table.format_inclination_cell(meta, refs))

    def test_only_individually_discussed_source_names_are_bold(self) -> None:
        expected_keys = {
            "NGC4258",
            "NGC1194",
            "NGC1068",
            "NGC3079",
            "NGC4945",
            "Circinus",
        }
        self.assertEqual(table.INDIVIDUALLY_DISCUSSED_SOURCE_KEYS, expected_keys)
        for meta in self.metadata:
            source_cell = self.sample_rows[meta["source_key"]][0]
            expected = meta["source_tex"]
            if meta["source_key"] in expected_keys:
                expected = rf"\textbf{{{expected}}}"
            self.assertEqual(source_cell, expected, meta["source_key"])
        self.assertEqual(self.rendered.count(r"\textbf{"), 6)
        self.assertIn(
            r"Boldface identifies the six sources discussed individually in "
            r"Section~\ref{sec:individual-sources}.",
            self.rendered,
        )

    def test_inclination_display_precision_and_qualifiers_are_audited(self) -> None:
        self.assertEqual(
            self.meta_by_key["NGC5765b"]["inclination_tex"],
            r"(0):72.4\pm0.5",
        )
        self.assertEqual(
            self.meta_by_key["NGC6264"]["inclination_tex"],
            r"(0):91.3\pm2.3",
        )
        self.assertEqual(
            self.meta_by_key["UGC6093"]["inclination_tex"],
            r"93.8_{-3.8}^{+4.1}",
        )
        self.assertEqual(
            self.meta_by_key["ESO558-G009"]["inclination_tex"],
            r"98\pm1",
        )
        self.assertEqual(
            self.meta_by_key["NGC4258"]["inclination_tex"],
            r"(6.1\,\mathrm{mas}):87.05\pm0.09",
        )
        self.assertEqual(
            self.meta_by_key["UGC3789"]["inclination_tex"],
            r"(0):84.9\pm0.6",
        )
        self.assertEqual(
            self.meta_by_key["NGC4388"]["inclination_tex"],
            r"\mathrm{unconstr.}",
        )
        self.assertEqual(
            self.meta_by_key["J0437"]["inclination_tex"],
            r"\mathrm{unconstr.}",
        )

    def test_revised_reference_keys_are_assigned_to_the_right_cells(self) -> None:
        self.assertEqual(
            self.meta_by_key["NGC1068"]["mass_reference_key"], "Gallimore2024"
        )
        self.assertEqual(
            self.meta_by_key["NGC1068"]["distance_reference_key"], "Markham2026"
        )
        self.assertEqual(self.meta_by_key["IC2560"]["mass_reference_key"], "Nakai2024")
        self.assertEqual(
            self.meta_by_key["IC2560"]["distance_reference_key"], "Nakai2024"
        )
        self.assertEqual(
            self.meta_by_key["NGC2960"]["mass_reference_key"], "Kuo2018Disk"
        )
        self.assertEqual(
            self.meta_by_key["J0437"]["mass_reference_key"], "Pesce2021J0437"
        )
        self.assertEqual(
            self.meta_by_key["J0437"]["distance_reference_key"], "Gao2016b"
        )
        self.assertEqual(
            self.meta_by_key["J0437"]["inclination_reference_key"],
            "Pesce2021J0437",
        )
        self.assertEqual(
            self.meta_by_key["NGC5765b"]["inclination_reference_key"],
            "Pesce2020b",
        )
        self.assertEqual(
            self.meta_by_key["NGC6264"]["inclination_reference_key"],
            "Pesce2020b",
        )
        self.assertEqual(
            self.meta_by_key["UGC3789"]["inclination_reference_key"],
            "Pesce2020b",
        )

    def test_ngc5765b_updated_fit_row_and_metadata_fail_closed(self) -> None:
        self.assertEqual(
            " & ".join(self.sample_rows["NGC5765b"]) + r" \\",
            r"NGC~5765b & $0.03792\pm 0.00013$ & "
            r"$\pm 0.00038$\tablenotemark{a} & "
            r"$^{14}$\,$7.618\pm 0.020$ & "
            r"$^{14}$\,$112.2\pm 5.3$ & "
            r"$^{14}$\,$72.4\pm0.5$ (at center) \\",
        )

        mutations = (
            ("inclination_tex", r"(0):94.50\pm0.25"),
            ("inclination_reference_key", "Gao2016a"),
            ("mass_reference_key", "Reid2019"),
            ("distance_reference_key", "Reid2019"),
        )
        for field, value in mutations:
            with self.subTest(field=field):
                metadata = [dict(item) for item in self.metadata]
                next(item for item in metadata if item["source_key"] == "NGC5765b")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "NGC5765b"):
                    table.validate_inputs(self.params, metadata)

    def test_ngc6264_updated_fit_row_and_metadata_fail_closed(self) -> None:
        row = self.by_key["NGC6264"]
        self.assertEqual(row["distance_stat_unc_minus_mpc"], "17.300000000000001")
        self.assertEqual(row["distance_stat_unc_plus_mpc"], "21.199999999999999")
        self.assertEqual(row["distance_unc_minus_mpc"], "17.300000000000001")
        self.assertEqual(row["distance_unc_plus_mpc"], "21.199999999999999")
        self.assertEqual(
            " & ".join(self.sample_rows["NGC6264"]) + r" \\",
            r"NGC~6264 & $0.021456\pm 0.000055$ & "
            r"$\pm 0.00021$\tablenotemark{a} & "
            r"$^{14}$\,$7.441\pm 0.064$ & "
            r"$^{14}$\,$132\pm 19$ & "
            r"$^{14}$\,$91.3\pm2.3$ (at center) \\",
        )

        mutations = (
            ("inclination_tex", r"(0):89.5\pm0.9"),
            ("inclination_reference_key", "Kuo2013"),
            ("mass_reference_key", "Kuo2013"),
            ("distance_reference_key", "Kuo2013"),
        )
        for field, value in mutations:
            with self.subTest(field=field):
                metadata = [dict(item) for item in self.metadata]
                next(item for item in metadata if item["source_key"] == "NGC6264")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "NGC6264"):
                    table.validate_inputs(self.params, metadata)

    def test_ugc3789_updated_fit_row_and_metadata_fail_closed(self) -> None:
        self.assertEqual(
            " & ".join(self.sample_rows["UGC3789"]) + r" \\",
            r"UGC~3789 & $0.023706\pm 0.000081$ & "
            r"$\pm 0.00024$\tablenotemark{a} & "
            r"$^{14}$\,$7.076\pm 0.035$ & "
            r"$^{14}$\,$51.5\pm 4.3$ & "
            r"$^{14}$\,$84.9\pm0.6$ (at center) \\",
        )

        mutations = (
            ("inclination_tex", r"(0.60\,\mathrm{mas}):90.6\pm0.4"),
            ("inclination_reference_key", "Reid2013"),
            ("mass_reference_key", "Reid2013"),
            ("distance_reference_key", "Reid2013"),
        )
        for field, value in mutations:
            with self.subTest(field=field):
                metadata = [dict(item) for item in self.metadata]
                next(item for item in metadata if item["source_key"] == "UGC3789")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "UGC3789"):
                    table.validate_inputs(self.params, metadata)

    def test_six_author_supplied_rows_separate_statistical_and_systematic(self) -> None:
        paired = {
            "NGC4258",
            "NGC5765b",
            "CGCG074-064",
            "NGC6264",
            "UGC3789",
            "NGC6323",
        }
        for key in paired:
            central_cell = self.sample_rows[key][1]
            status_cell = self.sample_rows[key][2]
            self.assertIn(r"\pm", central_cell, key)
            self.assertRegex(
                status_cell,
                r"^\$\\pm [0-9.]+\$\\tablenotemark\{a\}$",
                key,
            )
        self.assertIn(
            "Systematic/model term",
            self.rendered,
        )
        self.assertIn(
            "exact bounds and incomplete components remain machine-readable",
            self.rendered,
        )

    def test_ngc4258_data_selection_range_is_separate_from_the_cell(self) -> None:
        cell = self.sample_rows["NGC4258"][1]
        status = self.sample_rows["NGC4258"][2]
        self.assertEqual(cell, r"$0.53860\pm 0.00059$")
        self.assertEqual(
            status,
            r"$\pm 0.0054$\tablenotemark{a}",
        )
        self.assertNotIn("0.5368", cell)
        self.assertNotIn("0.5489", cell)
        self.assertNotIn(r"\mathrm{(model)}", cell)
        self.assertNotIn(r"\tablenotemark", cell)
        self.assertNotIn(r"0.5368\text{--}0.5489", self.rendered)
        self.assertIn(
            "exact bounds and incomplete components remain machine-readable",
            self.rendered,
        )

        row = dict(self.by_key["NGC4258"])
        row["size_model_check_probability_status"] = "probabilistic"
        params = [row if item["source_key"] == "NGC4258" else item for item in self.params]
        with self.assertRaisesRegex(ValueError, "must be non-probabilistic"):
            table.validate_inputs(params, self.metadata)

        mutations = (
            ("size_model_check_data_relation", "unknown", "data relation"),
            ("size_model_check_range_basis", "unknown", "range basis"),
            ("size_model_check_accepted_configuration_count", "14", "count must be 15"),
            ("size_model_check_excluded_configuration_ids", "2;16", "must be 2;17"),
            ("size_model_check_reference_key", "Reid2019", "citation does not match"),
            ("size_model_check_lower_uas", "0.5367", "reviewed endpoints"),
        )
        for field, value, message in mutations:
            with self.subTest(field=field):
                mutated = dict(self.by_key["NGC4258"])
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "NGC4258" else item
                    for item in self.params
                ]
                with self.assertRaisesRegex(ValueError, message):
                    table.validate_inputs(params, self.metadata)

    def test_ngc3079_turbulent_support_check_is_separate_from_the_cell(self) -> None:
        cell = self.sample_rows["NGC3079"][1]
        self.assertEqual(cell, r"$0.0119$")
        self.assertEqual(
            self.sample_rows["NGC3079"][2], r"$+0.0040$\tablenotemark{b}"
        )
        self.assertEqual(
            " & ".join(self.sample_rows["NGC3079"]) + r" \\",
            r"\textbf{NGC~3079} & $0.0119$ & $+0.0040$\tablenotemark{b} & "
            r"$^{19}$\,$6.30$ & "
            r"$^{19}$\,$17.3\pm 1.5$ & "
            r"$^{19}$\,thick, self-gravitating disk \\",
        )
        self.assertNotIn("0.0159", cell)
        self.assertNotIn(r"\mathrm{(model)}", cell)
        self.assertNotIn(r"0.0119\text{--}0.0159", self.rendered)
        self.assertIn("exact bounds and incomplete components remain machine-readable", self.rendered)

        row = self.by_key["NGC3079"]
        self.assertEqual(row["size_uas"], "0.011858795677384631")
        self.assertEqual(row["size_unc_adopted_uas"], "NA")
        self.assertEqual(row["size_stat_unc_minus_uas"], "NA")
        self.assertEqual(row["size_systematic_unc_minus_uas"], "NA")
        self.assertEqual(row["size_total_unc_minus_uas"], "NA")
        self.assertEqual(row["size_stat_component_status"], "no_formal_uncertainty")
        self.assertEqual(row["size_systematic_status"], "unquantified_model")
        self.assertEqual(
            row["size_total_status"], "incomplete_unquantified_systematic"
        )
        self.assertEqual(
            row["distance_stat_component_status"],
            "reported_legacy_tully_fisher",
        )
        self.assertEqual(row["distance_stat_unc_minus_mpc"], "1.5")
        self.assertEqual(row["distance_stat_unc_plus_mpc"], "1.5")
        self.assertEqual(
            row["mass_shared_distance_unc_minus_msun"],
            "173410.40462427746",
        )
        self.assertEqual(
            row["mass_shared_distance_unc_plus_msun"],
            "173410.40462427746",
        )
        self.assertEqual(row["mass_total_unc_minus_msun"], "NA")
        self.assertEqual(row["mass_total_unc_plus_msun"], "NA")
        self.assertEqual(row["mass_display_unc_minus_msun"], "NA")
        self.assertEqual(row["mass_display_unc_plus_msun"], "NA")
        self.assertEqual(row["mass_display_interval_status"], "not_available")

        metadata = self.meta_by_key["NGC3079"]
        self.assertEqual(metadata["source_tex"], "NGC~3079")
        self.assertEqual(metadata["inclination_tex"], r"\mathrm{t/s}")
        self.assertEqual(metadata["mass_reference_key"], "Kondratko2005")
        self.assertEqual(metadata["distance_reference_key"], "Kondratko2005")
        self.assertEqual(metadata["inclination_reference_key"], "Kondratko2005")

        primary_mutations = (
            ("size_uas", "0.0119"),
            ("mass_msun", "2000001"),
            ("distance_stat_component_status", "statistical"),
            ("mass_fit_unc_minus_msun", "1"),
            ("mass_shared_distance_unc_plus_msun", "173411"),
            ("mass_total_unc_minus_msun", "173410.40462427746"),
            ("mass_display_unc_plus_msun", "173410.40462427746"),
            ("size_stat_unc_minus_uas", "0.0001"),
            ("size_systematic_unc_plus_uas", "0.0001"),
            ("size_total_unc_plus_uas", "0.0001"),
            ("size_stat_component_status", "fixed_distance_formal_fit_interval"),
            ("size_systematic_status", "no_formal_uncertainty"),
            ("size_total_status", "no_formal_uncertainty"),
        )
        for field, value in primary_mutations:
            with self.subTest(primary_field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "NGC3079" else item
                    for item in self.params
                ]
                with self.assertRaisesRegex(ValueError, "NGC3079"):
                    table.validate_inputs(params, self.metadata)

        metadata_mutations = (
            ("source_tex", "NGC 3079"),
            ("inclination_tex", r"84\pm0"),
            ("mass_reference_key", "Yamauchi2004"),
            ("distance_reference_key", "Yamauchi2004"),
            ("inclination_reference_key", "Yamauchi2004"),
        )
        for field, value in metadata_mutations:
            with self.subTest(metadata_field=field):
                mutated_metadata = [dict(item) for item in self.metadata]
                next(
                    item
                    for item in mutated_metadata
                    if item["source_key"] == "NGC3079"
                )[field] = value
                with self.assertRaisesRegex(ValueError, "NGC3079"):
                    table.validate_inputs(self.params, mutated_metadata)

        mutations = (
            ("size_model_check_key", "different", "required turbulent-support check"),
            ("size_model_check_type", "unknown", "check type"),
            ("size_model_check_data_relation", "unknown", "data relation"),
            ("size_model_check_range_basis", "unknown", "range basis"),
            ("size_model_check_accepted_configuration_count", "1", "count must remain NA"),
            ("size_model_check_excluded_configuration_ids", "1", "must remain NA"),
            ("size_model_check_probability_status", "probabilistic", "non-probabilistic"),
            ("size_model_check_adoption_status", "adopted", "must remain separate"),
            ("size_model_check_reference_key", "vanDerBoom2026", "citation does not match"),
            ("size_model_check_reference_doi", "10.0000/wrong", "citation does not match"),
            ("size_model_check_lower_uas", "0.0118", "reviewed endpoints"),
            ("size_model_check_upper_uas", "0.0158", "reviewed endpoints"),
        )
        for field, value, message in mutations:
            with self.subTest(field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "NGC3079" else item
                    for item in self.params
                ]
                with self.assertRaisesRegex(ValueError, message):
                    table.validate_inputs(params, self.metadata)

    def test_ngc4945_turbulent_support_check_is_separate_from_the_cell(self) -> None:
        cell = self.sample_rows["NGC4945"][1]
        self.assertEqual(cell, r"$0.0388$")
        self.assertEqual(
            self.sample_rows["NGC4945"][2], r"$+0.0093$\tablenotemark{b}"
        )
        self.assertNotIn("0.0481", cell)
        self.assertNotIn(r"\mathrm{(model)}", cell)
        self.assertNotIn(r"0.0388\text{--}0.0481", self.rendered)
        self.assertIn("exact bounds and incomplete components remain machine-readable", self.rendered)

        row = self.by_key["NGC4945"]
        self.assertEqual(row["size_uas"], "0.038813517744088623")
        self.assertEqual(row["mass_msun"], "1400000")
        self.assertEqual(row["distance_mpc"], "3.7000000000000002")
        self.assertEqual(row["parameter_reference_doi"], "10.1086/310643")
        self.assertEqual(row["mass_uncertainty_reference_doi"], "NA")
        self.assertEqual(row["size_unc_adopted_uas"], "NA")
        self.assertEqual(row["size_stat_unc_minus_uas"], "NA")
        self.assertEqual(row["size_systematic_unc_minus_uas"], "NA")
        self.assertEqual(row["size_total_unc_minus_uas"], "NA")
        self.assertEqual(row["size_stat_component_status"], "no_formal_uncertainty")
        self.assertEqual(row["size_systematic_status"], "no_formal_uncertainty")
        self.assertEqual(row["size_total_status"], "no_formal_uncertainty")

        mutations = (
            ("size_model_check_key", "different", "required turbulent-support check"),
            ("size_model_check_type", "unknown", "check type"),
            ("size_model_check_data_relation", "unknown", "data relation"),
            ("size_model_check_range_basis", "unknown", "range basis"),
            ("size_model_check_accepted_configuration_count", "1", "count must remain NA"),
            ("size_model_check_excluded_configuration_ids", "1", "must remain NA"),
            ("size_model_check_probability_status", "probabilistic", "non-probabilistic"),
            ("size_model_check_adoption_status", "adopted", "must remain separate"),
            ("size_model_check_reference_key", "vanDerBoom2026", "citation does not match"),
            ("size_model_check_reference_doi", "10.0000/wrong", "citation does not match"),
            ("size_model_check_lower_uas", "0.0388", "reviewed endpoints"),
            ("size_model_check_upper_uas", "0.0481", "reviewed endpoints"),
        )
        for field, value, message in mutations:
            with self.subTest(field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "NGC4945" else item
                    for item in self.params
                ]
                with self.assertRaisesRegex(ValueError, message):
                    table.validate_inputs(params, self.metadata)

    def test_ngc4388_components_and_literature_checks_remain_separate(self) -> None:
        row = self.by_key["NGC4388"]
        cell = self.sample_rows["NGC4388"][1]
        self.assertEqual(cell, r"$0.0454$")
        self.assertEqual(self.sample_rows["NGC4388"][2], "---")
        self.assertNotIn("0.04342", cell)
        self.assertNotIn("0.04540", cell)
        self.assertNotIn(r"\mathrm{(model)}", cell)
        self.assertEqual(
            strip_reference_prefix(self.sample_rows["NGC4388"][4]),
            r"$19.0\pm 2.1$",
        )
        self.assertEqual(
            strip_reference_prefix(self.sample_rows["NGC4388"][5]),
            "unconstrained",
        )

        self.assertEqual(row["size_uas"], "0.045350531258882491")
        self.assertEqual(row["mass_msun"], "8400000")
        self.assertEqual(row["distance_mpc"], "19")
        self.assertEqual(row["distance_stat_unc_minus_mpc"], "NA")
        self.assertEqual(row["distance_stat_component_status"], "not_available")
        self.assertEqual(row["distance_systematic_unc_minus_mpc"], "2.0899999999999999")
        self.assertEqual(
            row["size_stat_component_status"],
            "fit_plus_absolute_position_mixed_component",
        )
        self.assertEqual(row["size_total_status"], "incomplete_unquantified_systematic")
        self.assertEqual(
            row["size_model_check_key"],
            "kuo_2020_ngc4388_turbulent_support_bias",
        )
        self.assertEqual(row["size_model_check_lower_uas"], "0.045350531258882491")
        self.assertEqual(row["size_model_check_upper_uas"], "0.045395881790141368")
        self.assertEqual(row["size_model_check_probability_status"], "non_probabilistic")
        self.assertEqual(
            row["size_model_check_adoption_status"],
            "reported_separately_not_adopted",
        )
        self.assertEqual(
            row["nonadopted_fit_ratio_msun_per_mpc"], "423300"
        )
        self.assertEqual(row["nonadopted_fit_size_uas"], "0.043421514018549316")
        self.assertEqual(
            row["nonadopted_fit_adoption_status"],
            "reported_separately_not_adopted",
        )

        self.assertNotIn(r"0.04535\text{--}0.04540", self.rendered)
        self.assertIn("exact bounds and incomplete components remain machine-readable", self.rendered)

        mutations = (
            ("distance_stat_component_status", "statistical", "requires numerical bounds"),
            ("distance_stat_unc_minus_mpc", "2.09", "incomplete distance interval"),
            ("distance_systematic_unc_minus_mpc", "2.08", "must be 2.09 Mpc"),
            ("size_stat_component_status", "fixed_distance_formal_fit_interval", "fit plus position"),
            ("size_model_check_key", "different", "required turbulent-support check"),
            ("size_model_check_type", "unknown", "check type"),
            ("size_model_check_data_relation", "unknown", "data relation"),
            ("size_model_check_range_basis", "unknown", "range basis"),
            ("size_model_check_accepted_configuration_count", "1", "count must remain NA"),
            ("size_model_check_excluded_configuration_ids", "1", "must remain NA"),
            ("size_model_check_probability_status", "probabilistic", "non-probabilistic"),
            ("size_model_check_adoption_status", "adopted", "must remain separate"),
            ("size_model_check_reference_key", "vanDerBoom2026", "citation does not match"),
            ("size_model_check_reference_doi", "10.0000/wrong", "citation does not match"),
            ("size_model_check_lower_uas", "0.0453", "reviewed endpoints"),
            ("size_model_check_upper_uas", "0.0454", "reviewed endpoints"),
            ("nonadopted_fit_ratio_msun_per_mpc", "423301", "reviewed non-adopted fit"),
            ("nonadopted_fit_basis", "none", "non-adoption basis is invalid"),
            ("nonadopted_fit_reference_doi", "10.0000/wrong", "citation does not match"),
            ("treatment", "no_robust_mass", "requires shared_distance"),
            ("mass_msun", "NA", "requires retained primary cells"),
            ("size_uas", "NA", "requires retained primary cells"),
        )
        for field, value, message in mutations:
            with self.subTest(field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "NGC4388" else item
                    for item in self.params
                ]
                with self.assertRaisesRegex(ValueError, message):
                    table.validate_inputs(params, self.metadata)

        metadata = [dict(item) for item in self.metadata]
        next(item for item in metadata if item["source_key"] == "NGC4388")[
            "inclination_tex"
        ] = r"\mathrm{e.o.}"
        with self.assertRaisesRegex(ValueError, "inclination must remain unconstrained"):
            table.validate_inputs(self.params, metadata)

    def test_numeric_model_components_and_caveats_are_explicit(self) -> None:
        self.assertEqual(
            strip_table_note_mark(self.sample_rows["NGC3393"][2]),
            r"$\pm 0.0091$",
        )
        self.assertEqual(
            strip_table_note_mark(self.sample_rows["UGC6093"][2]),
            r"$\pm 0.00087$",
        )
        self.assertIn(
            "Systematic/model term",
            self.rendered,
        )
        self.assertEqual(self.rendered.count(r"^{\dagger}"), 4)
        self.assertNotIn("pending", self.rendered.lower())

    def test_primary_uncertainty_markers_are_exact_and_fail_closed(self) -> None:
        self.assertEqual(
            {
                key: row["size_uncertainty_marker"]
                for key, row in self.meta_by_key.items()
                if row["size_uncertainty_marker"]
            },
            {
                "NGC1194": "dagger",
                "NGC2273": "dagger",
                "ESO558-G009": "dagger",
            },
        )

        mutations = (
            ("NGC1194", "", "dagger set must match"),
            ("NGC4258", "dagger", "dagger set must match"),
            ("NGC4258", "asterisk", "unsupported size uncertainty markers"),
        )
        for key, marker, message in mutations:
            with self.subTest(key=key, marker=marker):
                metadata = [dict(item) for item in self.metadata]
                next(item for item in metadata if item["source_key"] == key)[
                    "size_uncertainty_marker"
                ] = marker
                with self.assertRaisesRegex(ValueError, message):
                    table.validate_inputs(self.params, metadata)

        with self.assertRaisesRegex(
            ValueError, "marker requires a displayed primary uncertainty"
        ):
            table.format_size_cell(self.by_key["NGC4945"], "dagger")

    def test_ic2560_full_row_metadata_and_note_are_exact_and_fail_closed(self) -> None:
        self.assertEqual(
            " & ".join(self.sample_rows["IC2560"]) + r" \\",
            r"IC~2560 & $0.00954$ & --- & "
            r"$^{20}$\,$6.617\pm 0.063$ & "
            r"$^{20}$\,$44.5\pm 6.1$ & "
            r"$^{21}$\,nearly edge-on \\",
        )
        metadata = self.meta_by_key["IC2560"]
        self.assertEqual(metadata["mass_reference_key"], "Nakai2024")
        self.assertEqual(metadata["distance_reference_key"], "Nakai2024")
        self.assertEqual(metadata["inclination_reference_key"], "Yamauchi2012")
        self.assertEqual(metadata["inclination_tex"], r"\mathrm{e.o.}")
        self.assertNotIn("IC2560", {row["source_key"] for row in table.load_params() if row["nonadopted_fit_key"] not in {"", "NA"}})

        for field, value in (
            ("mass_reference_key", "Yamauchi2012"),
            ("distance_reference_key", "Yamauchi2012"),
            ("inclination_reference_key", "Nakai2024"),
            ("inclination_tex", r"90\pm0"),
        ):
            with self.subTest(metadata_field=field):
                mutated_metadata = [dict(item) for item in self.metadata]
                next(
                    item for item in mutated_metadata if item["source_key"] == "IC2560"
                )[field] = value
                with self.assertRaisesRegex(ValueError, "IC2560"):
                    table.validate_inputs(self.params, mutated_metadata)

    def test_circinus_uses_reviewed_radius_mass_distance_and_model_semantics(self) -> None:
        row = self.by_key["Circinus"]
        self.assertEqual(
            self.sample_rows["Circinus"][1],
            r"$0.04152$",
        )
        self.assertEqual(
            self.sample_rows["Circinus"][2], r"$+0.0012$\tablenotemark{b}"
        )
        self.assertEqual(
            strip_reference_prefix(self.sample_rows["Circinus"][3]),
            r"$6.230\pm 0.078$",
        )
        self.assertEqual(
            strip_reference_prefix(self.sample_rows["Circinus"][4]),
            r"$4.20\pm 0.80$",
        )
        self.assertEqual(row["mass_display_interval_status"], "source_reported_total")
        self.assertEqual(row["mass_display_unc_minus_msun"], "300000")
        self.assertEqual(row["mass_display_unc_plus_msun"], "300000")
        self.assertEqual(row["mass_total_unc_minus_msun"], "325394.28943115612")
        self.assertEqual(
            row["distance_stat_component_status"], "reported_legacy_multimethod"
        )
        self.assertEqual(row["size_stat_component_status"], "angular_radius_only_component")
        self.assertNotIn(r"0.0415\text{--}0.0428", self.rendered)
        self.assertIn(
            "incomplete components remain machine-readable",
            self.rendered,
        )
        self.assertIn("exact bounds and incomplete components remain machine-readable", self.rendered)

        mutations = (
            ("mass_display_interval_status", "decomposed_components", "decomposed mass display"),
            ("mass_display_unc_minus_msun", "325394.28943115612", "reviewed bounds"),
            ("distance_stat_component_status", "statistical", "must be labeled reported"),
            ("size_model_check_key", "different", "required turbulent-support check"),
            ("size_model_check_data_relation", "unknown", "data relation"),
            ("size_model_check_upper_uas", "0.0427", "reviewed endpoints"),
            ("size_model_check_probability_status", "probabilistic", "non-probabilistic"),
            ("size_model_check_adoption_status", "adopted", "must remain separate"),
        )
        for field, value, message in mutations:
            with self.subTest(field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "Circinus" else item
                    for item in self.params
                ]
                with self.assertRaisesRegex(ValueError, message):
                    table.validate_inputs(params, self.metadata)

        with self.assertRaisesRegex(ValueError, "unsupported size statistical"):
            table.size_stat_label("angular_radius_formal_term")

    def test_eso558_mixed_component_is_named_and_defined(self) -> None:
        self.assertEqual(
            strip_table_note_mark(self.sample_rows["ESO558-G009"][2]),
            r"$\pm 0.00097$",
        )
        self.assertEqual(
            self.by_key["ESO558-G009"]["size_stat_component_status"],
            "fit_plus_absolute_position_mixed_component",
        )
        self.assertIn(
            "Systematic/model term",
            self.rendered,
        )

    def test_ngc2960_formal_astrometric_disk_and_comparison_contract_is_exact(
        self,
    ) -> None:
        row = self.by_key["NGC2960"]
        self.assertEqual(
            " & ".join(self.sample_rows["NGC2960"]) + r" \\",
            r"NGC~2960 & $0.017575$ & --- & "
            r"$^{16}$\,$7.087\pm 0.027$ & "
            r"$^{8}$\,$71.3\pm 4.5$ & "
            r"$^{8}$\,nearly edge-on \\",
        )

        self.assertEqual(row["size_uas"], "0.017574751618808785")
        self.assertEqual(row["size_stat_unc_minus_uas"], "8.5245359509177624e-05")
        self.assertEqual(row["size_stat_unc_plus_uas"], "9.9452919427373872e-05")
        self.assertEqual(row["size_systematic_unc_minus_uas"], "NA")
        self.assertEqual(row["size_systematic_unc_plus_uas"], "NA")
        self.assertEqual(row["size_total_unc_minus_uas"], "NA")
        self.assertEqual(row["size_total_unc_plus_uas"], "NA")
        self.assertEqual(row["size_systematic_status"], "unquantified_astrometric")
        self.assertEqual(
            row["size_total_status"], "incomplete_unquantified_systematic"
        )
        self.assertEqual(
            row["size_model_check_key"],
            "kuo_2018_ngc2960_physically_allowed_disk_mass_models",
        )
        self.assertEqual(row["size_model_check_type"], "disk_mass_model_choice")
        self.assertEqual(
            row["size_model_check_data_relation"],
            "same_fit_data_alternative_accretion_disk_models",
        )
        self.assertEqual(
            row["size_model_check_range_basis"],
            "published_physically_allowed_model_envelope",
        )
        self.assertEqual(row["size_model_check_accepted_configuration_count"], "3")
        self.assertEqual(row["size_model_check_excluded_configuration_ids"], "2;4")
        self.assertEqual(
            row["size_model_check_lower_uas"], "0.017503713819217805"
        )
        self.assertEqual(
            row["size_model_check_upper_uas"], "0.017574751618808785"
        )
        self.assertEqual(
            row["size_model_check_probability_status"], "non_probabilistic"
        )
        self.assertEqual(
            row["size_model_check_adoption_status"],
            "reported_separately_not_adopted",
        )
        self.assertEqual(row["size_model_check_reference_key"], "Kuo2018Disk")
        self.assertEqual(
            row["size_model_check_reference_doi"], "10.3847/1538-4357/aabff1"
        )

        self.assertEqual(
            row["nonadopted_fit_key"], "villaraos_2022_ngc2960_gr_ratio_fit"
        )
        self.assertEqual(row["nonadopted_fit_ratio_msun_per_mpc"], "173800")
        self.assertEqual(
            row["nonadopted_fit_ratio_unc_minus_msun_per_mpc"], "1600"
        )
        self.assertEqual(
            row["nonadopted_fit_ratio_unc_plus_msun_per_mpc"], "1500"
        )
        self.assertEqual(row["nonadopted_fit_size_uas"], "0.017828157657509735")
        self.assertEqual(
            row["nonadopted_fit_size_unc_minus_uas"], "0.00016412573217500329"
        )
        self.assertEqual(
            row["nonadopted_fit_size_unc_plus_uas"], "0.00015386787391406559"
        )
        self.assertEqual(
            row["nonadopted_fit_basis"],
            "older_published_data_reanalysis_conditional_on_fixed_kuo2011_"
            "inclination_vertical_center_and_selected_scatter_angle",
        )
        self.assertEqual(row["nonadopted_fit_reference_key"], "Villaraos2022")
        self.assertEqual(
            row["nonadopted_fit_reference_doi"], "10.1093/mnras/stac2973"
        )
        self.assertNotEqual(row["size_uas"], row["nonadopted_fit_size_uas"])

        self.assertNotIn(r"0.017504\text{--}0.017575", self.rendered)
        self.assertIn("exact bounds and incomplete components remain machine-readable", self.rendered)

    def test_ngc2960_table_model_comparison_and_metadata_mutations_fail_closed(
        self,
    ) -> None:
        row = self.by_key["NGC2960"]
        row_mutations = (
            ("size_model_check_key", "different"),
            ("size_model_check_type", "turbulent_support_mass_bias"),
            ("size_model_check_data_relation", "unknown"),
            ("size_model_check_range_basis", "unknown"),
            ("size_model_check_accepted_configuration_count", "2"),
            ("size_model_check_excluded_configuration_ids", "2;5"),
            ("size_model_check_lower_uas", "0.01750"),
            ("size_model_check_upper_uas", "0.01758"),
            ("size_model_check_probability_status", "probabilistic"),
            ("size_model_check_adoption_status", "adopted"),
            ("size_model_check_reference_key", "Kuo2020"),
            ("size_model_check_reference_doi", "10.0000/wrong"),
            ("size_systematic_status", "unquantified_model"),
            ("size_total_status", "formal_component_only"),
            ("nonadopted_fit_ratio_msun_per_mpc", "173801"),
            ("nonadopted_fit_ratio_unc_minus_msun_per_mpc", "1601"),
            ("nonadopted_fit_ratio_unc_plus_msun_per_mpc", "1501"),
            ("nonadopted_fit_basis", "none"),
            ("nonadopted_fit_reference_doi", "10.0000/wrong"),
            ("treatment", "no_robust_mass"),
            ("mass_msun", "NA"),
            ("size_uas", "NA"),
        )
        for field, value in row_mutations:
            with self.subTest(row_field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "NGC2960" else item
                    for item in self.params
                ]
                with self.assertRaises(ValueError):
                    table.validate_inputs(params, self.metadata)

        metadata_mutations = (
            ("mass_reference_key", "Kuo2010"),
            ("distance_reference_key", "Kuo2018Disk"),
            ("inclination_reference_key", "Kuo2018Disk"),
            ("inclination_tex", r"89"),
        )
        for field, value in metadata_mutations:
            with self.subTest(metadata_field=field):
                metadata = [dict(item) for item in self.metadata]
                next(item for item in metadata if item["source_key"] == "NGC2960")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "NGC2960"):
                    table.validate_inputs(self.params, metadata)

    def test_ngc6323_paired_primary_model_check_and_comparison_contract_is_exact(
        self,
    ) -> None:
        row = self.by_key["NGC6323"]
        meta = self.meta_by_key["NGC6323"]
        self.assertEqual(
            " & ".join(self.sample_rows["NGC6323"]) + r" \\",
            r"NGC~6323 & $0.009549\pm 0.000059$ & "
            r"$\pm 0.000095$\tablenotemark{a} & "
            r"$^{14}$\,$7.01\pm 0.12$ & "
            r"$^{14}$\,$109\pm 29$ & "
            r"$^{14}$\,$91.5\pm0.3$ (at center) \\",
        )
        self.assertEqual(
            meta,
            {
                "source_order": "18",
                "source_key": "NGC6323",
                "source_tex": "NGC~6323",
                "inclination_tex": r"(0):91.5\pm0.3",
                "mass_reference_key": "Pesce2020b",
                "distance_reference_key": "Pesce2020b",
                "inclination_reference_key": "Pesce2020b",
                "size_uncertainty_marker": "",
            },
        )

        expected_fields = {
            "size_uas": "0.0095494889224502499",
            "size_posterior_q16_uas": "0.0094933704762850176",
            "size_posterior_q50_uas": "0.0095494889224502499",
            "size_posterior_q84_uas": "0.0096108035210381859",
            "size_stat_unc_minus_uas": "5.611844616523233e-05",
            "size_stat_unc_plus_uas": "6.1314598587935976e-05",
            "size_systematic_unc_minus_uas": "9.5494889224502508e-05",
            "size_systematic_unc_plus_uas": "9.5494889224502508e-05",
            "size_total_unc_minus_uas": "0.00011076350422409036",
            "size_total_unc_plus_uas": "0.00011348459749234575",
            "mass_msun": "10200000",
            "distance_mpc": "109.40000000000001",
            "distance_stat_unc_minus_mpc": "23.399999999999999",
            "distance_stat_unc_plus_mpc": "34.200000000000003",
            "distance_unc_minus_mpc": "23.399999999999999",
            "distance_unc_plus_mpc": "34.200000000000003",
            "treatment": "paired_posterior",
            "size_stat_component_status": "paired_chain_quantile_interval",
            "size_systematic_status": "numeric_one_percent_non_keplerian_allowance",
            "size_total_status": "quadrature_statistical_plus_systematic",
            "size_model_check_key": (
                "kuo_2018_ngc6323_physically_allowed_disk_mass_models"
            ),
            "size_model_check_type": "disk_mass_model_choice",
            "size_model_check_data_relation": (
                "predecessor_fit_data_alternative_accretion_disk_models"
            ),
            "size_model_check_range_basis": (
                "published_physically_allowed_model_envelope"
            ),
            "size_model_check_accepted_configuration_count": "3",
            "size_model_check_excluded_configuration_ids": "2;4",
            "size_model_check_lower_uas": "0.0094599180213333686",
            "size_model_check_upper_uas": "0.0095368279239458342",
            "size_model_check_probability_status": "non_probabilistic",
            "size_model_check_adoption_status": "reported_separately_not_adopted",
            "size_model_check_reference_key": "Kuo2018Disk",
            "size_model_check_reference_doi": "10.3847/1538-4357/aabff1",
            "nonadopted_fit_key": "villaraos_2022_ngc6323_gr_ratio_fit",
            "nonadopted_fit_ratio_msun_per_mpc": "91600",
            "nonadopted_fit_ratio_unc_minus_msun_per_mpc": "500",
            "nonadopted_fit_ratio_unc_plus_msun_per_mpc": "500",
            "nonadopted_fit_size_uas": "0.00939619816701894",
            "nonadopted_fit_size_unc_minus_uas": "5.1289291304688531e-05",
            "nonadopted_fit_size_unc_plus_uas": "5.1289291304688531e-05",
            "nonadopted_fit_basis": (
                "older_partial_data_reanalysis_conditional_on_fixed_kuo2015_"
                "inclination_and_selected_scatter_angle"
            ),
            "nonadopted_fit_reference_key": "Villaraos2022",
            "nonadopted_fit_reference_doi": "10.1093/mnras/stac2973",
        }
        for field, expected in expected_fields.items():
            self.assertEqual(row[field], expected, field)
        self.assertAlmostEqual(
            100.0
            * (
                float(row["nonadopted_fit_size_uas"]) / float(row["size_uas"])
                - 1.0
            ),
            -1.6052247054911262,
            places=14,
        )

        self.assertNotIn(r"0.0094599\text{--}0.0095368", self.rendered)
        self.assertIn("exact bounds and incomplete components remain machine-readable", self.rendered)
        self.assertNotIn("Kuo2015", table.reference_order(self.metadata))

    def test_ngc6323_table_contract_mutations_fail_closed(self) -> None:
        row = self.by_key["NGC6323"]
        row_mutations = (
            ("size_uas", "0.00955"),
            ("size_stat_unc_minus_uas", "0.00005"),
            ("size_systematic_unc_plus_uas", "0.00010"),
            ("mass_msun", "10200001"),
            ("distance_mpc", "109.5"),
            ("size_model_check_key", "different"),
            ("size_model_check_type", "turbulent_support_mass_bias"),
            ("size_model_check_data_relation", "unknown"),
            ("size_model_check_range_basis", "unknown"),
            ("size_model_check_accepted_configuration_count", "2"),
            ("size_model_check_excluded_configuration_ids", "2;5"),
            ("size_model_check_lower_uas", "0.00945"),
            ("size_model_check_upper_uas", "0.00954"),
            ("size_model_check_probability_status", "probabilistic"),
            ("size_model_check_adoption_status", "adopted"),
            ("size_model_check_reference_key", "Kuo2020"),
            ("size_model_check_reference_doi", "10.0000/wrong"),
            ("nonadopted_fit_key", "different"),
            ("nonadopted_fit_ratio_msun_per_mpc", "91601"),
            ("nonadopted_fit_ratio_unc_minus_msun_per_mpc", "501"),
            ("nonadopted_fit_size_uas", "0.00940"),
            ("nonadopted_fit_basis", "none"),
            ("nonadopted_fit_reference_doi", "10.0000/wrong"),
            ("treatment", "shared_distance"),
        )
        for field, value in row_mutations:
            with self.subTest(row_field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "NGC6323" else item
                    for item in self.params
                ]
                with self.assertRaises(ValueError):
                    table.validate_inputs(params, self.metadata)

        metadata_mutations = (
            ("source_tex", "NGC 6323"),
            ("mass_reference_key", "Kuo2015"),
            ("distance_reference_key", "Kuo2015"),
            ("inclination_reference_key", "Kuo2015"),
            ("inclination_tex", r"(0):88.5_{-0.6}^{+0.7}"),
        )
        for field, value in metadata_mutations:
            with self.subTest(metadata_field=field):
                metadata = [dict(item) for item in self.metadata]
                next(item for item in metadata if item["source_key"] == "NGC6323")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "NGC6323"):
                    table.validate_inputs(self.params, metadata)

    def test_j0437_primary_comparison_full_row_and_metadata_are_exact(self) -> None:
        row = self.by_key["J0437"]
        meta = self.meta_by_key["J0437"]
        self.assertEqual(
            " & ".join(self.sample_rows["J0437"]) + r" \\",
            r"J0437+2456 & $0.00419$ & --- & "
            r"$^{22}$\,$6.421\pm 0.042$ & "
            r"$^{18}$\,$64.5\pm 4.3$ & "
            r"$^{22}$\,unconstrained \\",
        )
        self.assertEqual(
            meta,
            {
                "source_order": "19",
                "source_key": "J0437",
                "source_tex": "J0437+2456",
                "inclination_tex": r"\mathrm{unconstr.}",
                "mass_reference_key": "Pesce2021J0437",
                "distance_reference_key": "Gao2016b",
                "inclination_reference_key": "Pesce2021J0437",
                "size_uncertainty_marker": "",
            },
        )

        expected_fields = {
            "size_uas": "0.0041910678037545487",
            "size_stat_unc_minus_uas": "0.0002930816645982202",
            "size_stat_unc_plus_uas": "0.0002930816645982202",
            "size_unc_adopted_uas": "0.0002930816645982202",
            "mass_msun": "2635479.2190218437",
            "mass_fit_unc_minus_msun": "184299.24608544362",
            "mass_fit_unc_plus_msun": "184299.24608544362",
            "mass_shared_distance_unc_minus_msun": "176180.85157776746",
            "mass_shared_distance_unc_plus_msun": "176180.85157776746",
            "mass_total_unc_minus_msun": "254962.55523180307",
            "mass_total_unc_plus_msun": "254962.55523180307",
            "mass_display_unc_minus_msun": "254962.55523180307",
            "mass_display_unc_plus_msun": "254962.55523180307",
            "distance_mpc": "64.504736129905268",
            "distance_stat_unc_minus_mpc": "3.4259780517380682",
            "distance_stat_unc_plus_mpc": "3.4259780517380682",
            "distance_systematic_unc_minus_mpc": "2.6185955127160461",
            "distance_systematic_unc_plus_mpc": "2.6185955127160461",
            "distance_unc_minus_mpc": "4.3121187449103928",
            "distance_unc_plus_mpc": "4.3121187449103928",
            "treatment": "shared_distance",
            "mass_display_interval_status": "decomposed_components",
            "distance_stat_component_status": "statistical",
            "size_stat_component_status": "fixed_distance_formal_fit_interval",
            "size_systematic_status": "unquantified_model",
            "size_total_status": "incomplete_unquantified_systematic",
            "legacy_uncertainty_basis": (
                "symmetrized_formal_statistical_term_total_incomplete"
            ),
            "size_model_check_key": (
                "deepshikha_2024_j0437_mestel_disk_mass_model"
            ),
            "size_model_check_type": "disk_mass_model_choice",
            "size_model_check_data_relation": (
                "reanalysis_of_published_gao2017_vlbi_data_data_selection_"
                "unreported"
            ),
            "size_model_check_range_basis": (
                "reported_symmetric_mass_error_span"
            ),
            "size_model_check_accepted_configuration_count": "1",
            "size_model_check_excluded_configuration_ids": "NA",
            "size_model_check_lower_uas": "0.0024341189955645356",
            "size_model_check_upper_uas": "0.0031252885868976758",
            "size_model_check_probability_status": "non_probabilistic",
            "size_model_check_adoption_status": (
                "reported_separately_not_adopted"
            ),
            "size_model_check_reference_key": "Deepshikha2024",
            "size_model_check_reference_doi": "10.1017/S1743921323002491",
            "nonadopted_fit_key": "villaraos_2022_j0437_gr_ratio_fit",
            "nonadopted_fit_ratio_msun_per_mpc": "42500",
            "nonadopted_fit_ratio_unc_minus_msun_per_mpc": "1200",
            "nonadopted_fit_ratio_unc_plus_msun_per_mpc": "1000",
            "nonadopted_fit_size_uas": "0.0043595897608985253",
            "nonadopted_fit_size_unc_minus_uas": "0.00012309429913125248",
            "nonadopted_fit_size_unc_plus_uas": "0.00010257858260937706",
            "nonadopted_fit_interval_status": (
                "bayesian_formal_conditional_on_model"
            ),
            "nonadopted_fit_adoption_status": "reported_separately_not_adopted",
            "nonadopted_fit_basis": (
                "same_data_reanalysis_conditional_on_fixed_gao2017_inclination_"
                "vertical_center_and_zero_scatter_angle"
            ),
            "nonadopted_fit_reference_key": "Villaraos2022",
            "nonadopted_fit_reference_doi": "10.1093/mnras/stac2973",
        }
        for field, expected in expected_fields.items():
            self.assertEqual(row[field], expected, field)
        for field in (
            "size_systematic_unc_minus_uas",
            "size_systematic_unc_plus_uas",
            "size_total_unc_minus_uas",
            "size_total_unc_plus_uas",
        ):
            self.assertEqual(row[field], "NA", field)
        self.assertAlmostEqual(
            100.0
            * (
                float(row["nonadopted_fit_size_uas"]) / float(row["size_uas"])
                - 1.0
            ),
            4.020979020979021,
            places=14,
        )
        self.assertNotIn("0.00436", self.sample_rows["J0437"][1])
        self.assertNotIn("0.002434", self.sample_rows["J0437"][1])
        self.assertNotIn("0.003125", self.sample_rows["J0437"][1])
        self.assertNotIn(r"0.002434\text{--}0.003125", self.rendered)
        self.assertIn("exact bounds and incomplete components remain machine-readable", self.rendered)

    def test_j0437_table_contract_mutations_fail_closed(self) -> None:
        row = self.by_key["J0437"]
        row_mutations = (
            ("size_uas", "0.00420"),
            ("size_stat_unc_minus_uas", "0.00030"),
            ("mass_msun", "2635480"),
            ("mass_shared_distance_unc_minus_msun", "176181"),
            ("mass_display_unc_minus_msun", "254963"),
            ("distance_stat_unc_minus_mpc", "3.4261"),
            ("distance_unc_minus_mpc", "4.3122"),
            ("size_systematic_status", "numeric_systematic"),
            ("size_total_status", "formal_component_only"),
            ("size_total_unc_minus_uas", "0.001"),
            ("size_model_check_key", "different"),
            ("size_model_check_type", "data_selection_component_tracking"),
            ("size_model_check_data_relation", "reanalysis_of_other_data"),
            ("size_model_check_range_basis", "accepted_configuration_medians"),
            ("size_model_check_accepted_configuration_count", "2"),
            ("size_model_check_excluded_configuration_ids", "1"),
            ("size_model_check_lower_uas", "0.0025"),
            ("size_model_check_upper_uas", "0.0032"),
            ("size_model_check_probability_status", "probabilistic"),
            ("size_model_check_adoption_status", "adopted"),
            ("size_model_check_reference_key", "Gao2016b"),
            ("size_model_check_reference_doi", "10.0000/wrong"),
            ("nonadopted_fit_key", "different"),
            ("nonadopted_fit_ratio_msun_per_mpc", "42501"),
            ("nonadopted_fit_ratio_unc_minus_msun_per_mpc", "1201"),
            ("nonadopted_fit_size_uas", "0.00436"),
            ("nonadopted_fit_basis", "none"),
            ("nonadopted_fit_reference_doi", "10.0000/wrong"),
            ("treatment", "no_robust_mass"),
        )
        for field, value in row_mutations:
            with self.subTest(row_field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "J0437" else item
                    for item in self.params
                ]
                with self.assertRaisesRegex(ValueError, "J0437"):
                    table.validate_inputs(params, self.metadata)

        metadata_mutations = (
            ("source_tex", "J0437"),
            ("mass_reference_key", "Gao2016b"),
            ("distance_reference_key", "Pesce2021J0437"),
            ("inclination_reference_key", "Gao2016b"),
            ("inclination_tex", r"\mathrm{e.o.}"),
        )
        for field, value in metadata_mutations:
            with self.subTest(metadata_field=field):
                metadata = [dict(item) for item in self.metadata]
                next(item for item in metadata if item["source_key"] == "J0437")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "J0437"):
                    table.validate_inputs(self.params, metadata)

    def test_mrk1029_no_primary_full_row_diagnostics_and_metadata_are_exact(
        self,
    ) -> None:
        row = self.by_key["Mrk1029"]
        meta = self.meta_by_key["Mrk1029"]
        self.assertEqual(
            " & ".join(self.sample_rows["Mrk1029"]) + r" \\",
            r"Mrk~1029 & \nodata & --- & \nodata & "
            r"$^{18}$\,$119.3\pm 5.9$ & "
            r"$^{18}$\,$79\pm2$ (geometry estimate) \\",
        )
        self.assertEqual(
            meta,
            {
                "source_order": "20",
                "source_key": "Mrk1029",
                "source_tex": "Mrk~1029",
                "inclination_tex": r"79\pm2\;\mathrm{(est.)}",
                "mass_reference_key": "Gao2016b",
                "distance_reference_key": "Gao2016b",
                "inclination_reference_key": "Gao2016b",
                "size_uncertainty_marker": "",
            },
        )

        expected_fields = {
            "distance_mpc": "119.32882273342354",
            "distance_stat_unc_minus_mpc": "3.4105504657071966",
            "distance_stat_unc_plus_mpc": "3.4105504657071966",
            "distance_systematic_unc_minus_mpc": "4.8442011935083977",
            "distance_systematic_unc_plus_mpc": "4.8442011935083977",
            "distance_unc_minus_mpc": "5.9243682939469391",
            "distance_unc_plus_mpc": "5.9243682939469391",
            "treatment": "no_robust_mass",
            "adopted_uncertainty_model": "not_available",
            "mass_display_interval_status": "not_available",
            "distance_stat_component_status": "statistical",
            "size_stat_component_status": "not_available",
            "size_systematic_status": "not_applicable",
            "size_total_status": "not_available",
            "legacy_uncertainty_basis": "not_available",
            "ratio_provenance": "no_robust_primary_mass",
            "size_model_check_key": (
                "deepshikha_2024_mrk1029_mestel_disk_mass_model"
            ),
            "size_model_check_type": "disk_mass_model_choice",
            "size_model_check_data_relation": (
                "reanalysis_of_published_gao2017_vlbi_data_data_selection_"
                "unreported"
            ),
            "size_model_check_range_basis": "reported_symmetric_mass_error_span",
            "size_model_check_accepted_configuration_count": "1",
            "size_model_check_excluded_configuration_ids": "NA",
            "size_model_check_lower_uas": "0.00033171554552000282",
            "size_model_check_upper_uas": "0.00059061548348683429",
            "size_model_check_probability_status": "non_probabilistic",
            "size_model_check_adoption_status": "reported_separately_not_adopted",
            "size_model_check_reference_key": "Deepshikha2024",
            "size_model_check_reference_doi": "10.1017/S1743921323002491",
            "nonadopted_fit_key": "gonzalez_juarez_2024_mrk1029_gr_ratio_fit",
            "nonadopted_fit_ratio_msun_per_mpc": "14400",
            "nonadopted_fit_ratio_unc_minus_msun_per_mpc": "1100",
            "nonadopted_fit_ratio_unc_plus_msun_per_mpc": "1100",
            "nonadopted_fit_size_uas": "0.0014771315895750299",
            "nonadopted_fit_size_unc_minus_uas": "0.00011283644087031476",
            "nonadopted_fit_size_unc_plus_uas": "0.00011283644087031476",
            "nonadopted_fit_interval_status": (
                "bayesian_formal_conditional_on_model"
            ),
            "nonadopted_fit_adoption_status": "reported_separately_not_adopted",
            "nonadopted_fit_basis": (
                "restricted_same_data_fit_conditional_on_fixed_gao2017_"
                "inclination_vertical_center_and_selected_scatter_angle_"
                "does_not_resolve_primary_complex_geometry_and_dynamics"
            ),
            "nonadopted_fit_reference_key": "GonzalezJuarez2024",
            "nonadopted_fit_reference_doi": "10.1051/0004-6361/202450098",
        }
        for field, expected in expected_fields.items():
            self.assertEqual(row[field], expected, field)
        for field in (
            "size_uas",
            "size_unc_adopted_uas",
            "log_mass",
            "mass_msun",
            "mass_fit_unc_minus_msun",
            "mass_fit_unc_plus_msun",
            "mass_total_unc_minus_msun",
            "mass_total_unc_plus_msun",
            "mass_display_unc_minus_msun",
            "mass_display_unc_plus_msun",
            "size_stat_unc_minus_uas",
            "size_stat_unc_plus_uas",
            "size_systematic_unc_minus_uas",
            "size_systematic_unc_plus_uas",
            "size_total_unc_minus_uas",
            "size_total_unc_plus_uas",
        ):
            self.assertEqual(row[field], "NA", field)

        self.assertEqual(self.sample_rows["Mrk1029"][1], r"\nodata")
        self.assertEqual(self.sample_rows["Mrk1029"][2], "---")
        self.assertEqual(self.sample_rows["Mrk1029"][3], r"\nodata")
        for token in ("0.0003317", "0.0005906", "0.001477"):
            self.assertNotIn(token, self.sample_rows["Mrk1029"][1])
            self.assertNotIn(token, self.sample_rows["Mrk1029"][3])
        self.assertEqual(self.rendered.count(r"\nodata"), 4)
        self.assertIn("exact bounds and incomplete components remain machine-readable", self.rendered)

    def test_mrk1029_table_and_metadata_mutations_fail_closed(self) -> None:
        row = self.by_key["Mrk1029"]
        row_mutations = (
            ("distance_mpc", "119.4"),
            ("distance_stat_unc_minus_mpc", "3.4106"),
            ("distance_stat_unc_plus_mpc", "3.4106"),
            ("distance_systematic_unc_minus_mpc", "4.8443"),
            ("distance_unc_minus_mpc", "5.9244"),
            ("size_uas", "0.0004"),
            ("mass_msun", "570000"),
            ("treatment", "shared_distance"),
            ("size_stat_component_status", "fixed_distance_formal_fit_interval"),
            ("size_systematic_status", "unquantified_model"),
            ("size_model_check_key", "different"),
            ("size_model_check_type", "data_selection_component_tracking"),
            ("size_model_check_data_relation", "reanalysis_of_other_data"),
            ("size_model_check_range_basis", "accepted_configuration_medians"),
            ("size_model_check_accepted_configuration_count", "2"),
            ("size_model_check_excluded_configuration_ids", "1"),
            ("size_model_check_lower_uas", "0.00034"),
            ("size_model_check_upper_uas", "0.00060"),
            ("size_model_check_probability_status", "probabilistic"),
            ("size_model_check_adoption_status", "adopted"),
            ("size_model_check_reference_key", "Gao2016b"),
            ("size_model_check_reference_doi", "10.0000/wrong"),
            ("nonadopted_fit_key", "different"),
            ("nonadopted_fit_ratio_msun_per_mpc", "14401"),
            ("nonadopted_fit_ratio_unc_minus_msun_per_mpc", "1101"),
            ("nonadopted_fit_size_uas", "0.00148"),
            ("nonadopted_fit_basis", "none"),
            ("nonadopted_fit_reference_key", "Villaraos2022"),
            ("nonadopted_fit_reference_doi", "10.0000/wrong"),
        )
        for field, value in row_mutations:
            with self.subTest(row_field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "Mrk1029" else item
                    for item in self.params
                ]
                with self.assertRaisesRegex(ValueError, "Mrk1029"):
                    table.validate_inputs(params, self.metadata)

        metadata_mutations = (
            ("source_tex", "Mrk 1029"),
            ("mass_reference_key", "Deepshikha2024"),
            ("distance_reference_key", "Deepshikha2024"),
            ("inclination_reference_key", "Deepshikha2024"),
            ("inclination_tex", r"\mathrm{e.o.}"),
        )
        for field, value in metadata_mutations:
            with self.subTest(metadata_field=field):
                metadata = [dict(item) for item in self.metadata]
                next(item for item in metadata if item["source_key"] == "Mrk1029")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "Mrk1029"):
                    table.validate_inputs(self.params, metadata)

    def test_ngc1320_no_primary_full_row_diagnostics_and_metadata_are_exact(
        self,
    ) -> None:
        row = self.by_key["NGC1320"]
        meta = self.meta_by_key["NGC1320"]
        self.assertEqual(
            " & ".join(self.sample_rows["NGC1320"]) + r" \\",
            r"NGC~1320 & \nodata & --- & \nodata & "
            r"$^{18}$\,$33.8\pm 3.7$ & "
            r"$^{18}$\,unconstrained \\",
        )
        self.assertEqual(
            meta,
            {
                "source_order": "21",
                "source_key": "NGC1320",
                "source_tex": "NGC~1320",
                "inclination_tex": r"\mathrm{unconstr.}",
                "mass_reference_key": "Gao2016b",
                "distance_reference_key": "Gao2016b",
                "inclination_reference_key": "Gao2016b",
                "size_uncertainty_marker": "",
            },
        )

        expected_fields = {
            "distance_mpc": "33.783491204330176",
            "distance_stat_unc_minus_mpc": "3.3898711337347889",
            "distance_stat_unc_plus_mpc": "3.3898711337347889",
            "distance_systematic_unc_minus_mpc": "1.3714543114071789",
            "distance_systematic_unc_plus_mpc": "1.3714543114071789",
            "distance_unc_minus_mpc": "3.6567900174341048",
            "distance_unc_plus_mpc": "3.6567900174341048",
            "treatment": "no_robust_mass",
            "adopted_uncertainty_model": "not_available",
            "mass_display_interval_status": "not_available",
            "distance_stat_component_status": "statistical",
            "size_stat_component_status": "not_available",
            "size_systematic_status": "not_applicable",
            "size_total_status": "not_available",
            "legacy_uncertainty_basis": "not_available",
            "ratio_provenance": "no_robust_primary_mass",
            "size_model_check_key": (
                "deepshikha_2024_ngc1320_mestel_disk_mass_model"
            ),
            "size_model_check_type": "disk_mass_model_choice",
            "size_model_check_data_relation": (
                "reanalysis_of_published_gao2017_vlbi_data_data_selection_"
                "unreported"
            ),
            "size_model_check_range_basis": "reported_symmetric_mass_error_span",
            "size_model_check_accepted_configuration_count": "1",
            "size_model_check_excluded_configuration_ids": "NA",
            "size_model_check_lower_uas": "0.011032540265018939",
            "size_model_check_upper_uas": "0.014131268281343874",
            "size_model_check_probability_status": "non_probabilistic",
            "size_model_check_adoption_status": "reported_separately_not_adopted",
            "size_model_check_reference_key": "Deepshikha2024",
            "size_model_check_reference_doi": "10.1017/S1743921323002491",
            "nonadopted_fit_key": "gonzalez_juarez_2024_ngc1320_gr_ratio_fit",
            "nonadopted_fit_ratio_msun_per_mpc": "149700",
            "nonadopted_fit_ratio_unc_minus_msun_per_mpc": "6900",
            "nonadopted_fit_ratio_unc_plus_msun_per_mpc": "7100",
            "nonadopted_fit_size_uas": "0.015356013816623746",
            "nonadopted_fit_size_unc_minus_uas": "0.0007077922200047018",
            "nonadopted_fit_size_unc_plus_uas": "0.00072830793652657719",
            "nonadopted_fit_interval_status": (
                "bayesian_formal_conditional_on_model"
            ),
            "nonadopted_fit_adoption_status": "reported_separately_not_adopted",
            "nonadopted_fit_basis": (
                "restricted_same_data_fit_conditional_on_fixed_edge_on_"
                "inclination_and_selected_scatter_angle_does_not_resolve_"
                "primary_complex_geometry_and_dynamics"
            ),
            "nonadopted_fit_reference_key": "GonzalezJuarez2024",
            "nonadopted_fit_reference_doi": "10.1051/0004-6361/202450098",
        }
        for field, expected in expected_fields.items():
            self.assertEqual(row[field], expected, field)
        for field in (
            "size_uas",
            "size_unc_adopted_uas",
            "log_mass",
            "mass_msun",
            "mass_fit_unc_minus_msun",
            "mass_fit_unc_plus_msun",
            "mass_shared_distance_unc_minus_msun",
            "mass_shared_distance_unc_plus_msun",
            "mass_total_unc_minus_msun",
            "mass_total_unc_plus_msun",
            "mass_display_unc_minus_msun",
            "mass_display_unc_plus_msun",
            "size_stat_unc_minus_uas",
            "size_stat_unc_plus_uas",
            "size_systematic_unc_minus_uas",
            "size_systematic_unc_plus_uas",
            "size_total_unc_minus_uas",
            "size_total_unc_plus_uas",
        ):
            self.assertEqual(row[field], "NA", field)

        self.assertEqual(self.sample_rows["NGC1320"][1], r"\nodata")
        self.assertEqual(self.sample_rows["NGC1320"][2], "---")
        self.assertEqual(self.sample_rows["NGC1320"][3], r"\nodata")
        for token in ("0.01103", "0.01413", "0.01536"):
            self.assertNotIn(token, self.sample_rows["NGC1320"][1])
            self.assertNotIn(token, self.sample_rows["NGC1320"][3])

        self.assertNotIn(r"0.01103\text{--}0.01413", self.rendered)
        self.assertEqual(self.rendered.count(r"\nodata"), 4)
        self.assertIn("exact bounds and incomplete components remain machine-readable", self.rendered)

    def test_ngc1320_table_and_metadata_mutations_fail_closed(self) -> None:
        row = self.by_key["NGC1320"]
        row_mutations = (
            ("distance_mpc", "33.8"),
            ("distance_stat_unc_minus_mpc", "3.39"),
            ("distance_stat_unc_plus_mpc", "3.39"),
            ("distance_systematic_unc_minus_mpc", "1.38"),
            ("distance_systematic_unc_plus_mpc", "1.38"),
            ("distance_unc_minus_mpc", "3.66"),
            ("distance_unc_plus_mpc", "3.66"),
            ("size_uas", "0.012"),
            ("mass_msun", "4710000"),
            ("treatment", "shared_distance"),
            ("mass_display_interval_status", "decomposed_components"),
            ("distance_stat_component_status", "reported"),
            ("size_stat_component_status", "fixed_distance_formal_fit_interval"),
            ("size_systematic_status", "unquantified_model"),
            ("size_total_status", "incomplete_unquantified_systematic"),
            ("size_model_check_key", "different"),
            ("size_model_check_type", "turbulent_support_mass_bias"),
            ("size_model_check_data_relation", "reanalysis_of_other_data"),
            ("size_model_check_range_basis", "accepted_configuration_medians"),
            ("size_model_check_accepted_configuration_count", "2"),
            ("size_model_check_excluded_configuration_ids", "1"),
            ("size_model_check_lower_uas", "0.0111"),
            ("size_model_check_upper_uas", "0.0142"),
            ("size_model_check_probability_status", "probabilistic"),
            ("size_model_check_adoption_status", "adopted"),
            ("size_model_check_reference_key", "Kuo2020"),
            ("size_model_check_reference_doi", "10.0000/wrong"),
            ("nonadopted_fit_key", "different"),
            ("nonadopted_fit_ratio_msun_per_mpc", "149701"),
            ("nonadopted_fit_ratio_unc_minus_msun_per_mpc", "6901"),
            ("nonadopted_fit_ratio_unc_plus_msun_per_mpc", "7101"),
            ("nonadopted_fit_size_uas", "0.0154"),
            ("nonadopted_fit_size_unc_minus_uas", "0.0008"),
            ("nonadopted_fit_size_unc_plus_uas", "0.0008"),
            ("nonadopted_fit_interval_status", "unconditional"),
            ("nonadopted_fit_adoption_status", "adopted"),
            ("nonadopted_fit_basis", "none"),
            ("nonadopted_fit_reference_key", "Villaraos2022"),
            ("nonadopted_fit_reference_doi", "10.0000/wrong"),
        )
        for field, value in row_mutations:
            with self.subTest(row_field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "NGC1320" else item
                    for item in self.params
                ]
                with self.assertRaisesRegex(ValueError, "NGC1320"):
                    table.validate_inputs(params, self.metadata)

        metadata_mutations = (
            ("source_tex", "NGC 1320"),
            ("mass_reference_key", "Deepshikha2024"),
            ("distance_reference_key", "Deepshikha2024"),
            ("inclination_reference_key", "Deepshikha2024"),
            ("inclination_tex", r"90\;\mathrm{(fixed)}"),
            ("inclination_tex", r"\mathrm{e.o.}"),
        )
        for field, value in metadata_mutations:
            with self.subTest(metadata_field=field):
                metadata = [dict(item) for item in self.metadata]
                next(item for item in metadata if item["source_key"] == "NGC1320")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "NGC1320"):
                    table.validate_inputs(self.params, metadata)

    def test_kuo_central_values_and_gao_reported_provenance_are_explicit(self) -> None:
        self.assertEqual(
            self.sample_rows["NGC1194"][1],
            r"$0.1253\pm 0.0077^{\dagger}$",
        )
        self.assertEqual(self.sample_rows["NGC1194"][2], "---")
        self.assertEqual(
            self.sample_rows["NGC2273"][1],
            r"$0.0299\pm 0.0020^{\dagger}$",
        )
        self.assertEqual(self.sample_rows["NGC2273"][2], "---")
        for key in ("NGC1194", "NGC2273"):
            self.assertEqual(
                self.by_key[key]["size_stat_component_status"],
                "reported_undecomposed_uncertainty",
            )
        self.assertNotIn("resid.", self.rendered)

    def test_unknown_size_component_status_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "unsupported size statistical"):
            table.size_stat_label("unknown_status")

    def test_no_errors_are_invented_for_approximate_mass_rows(self) -> None:
        expected_model_terms = {
            "NGC4945": r"$+0.0093$\tablenotemark{b}",
            "NGC3079": r"$+0.0040$\tablenotemark{b}",
        }
        for key, model_term in expected_model_terms.items():
            self.assertEqual(len(clean_numeric_tokens(self.sample_rows[key][1])), 1, key)
            self.assertEqual(len(clean_numeric_tokens(self.sample_rows[key][3])), 1, key)
            self.assertEqual(self.sample_rows[key][2], model_term)
            self.assertNotIn(r"\tablenotemark", self.sample_rows[key][1], key)

    def test_two_sources_have_neither_mass_nor_shadow_size(self) -> None:
        expected = {
            "Mrk1029": (14_400.0, 0.0014771315895750299),
            "NGC1320": (149_700.0, 0.015356013816623746),
        }
        for key, (ratio, implied_size) in expected.items():
            self.assertEqual(self.sample_rows[key][1], r"\nodata")
            self.assertEqual(self.sample_rows[key][2], "---")
            self.assertEqual(self.sample_rows[key][3], r"\nodata")
            self.assertEqual(
                table.optional_float(
                    self.by_key[key]["nonadopted_fit_ratio_msun_per_mpc"]
                ),
                ratio,
            )
            self.assertEqual(
                table.optional_float(self.by_key[key]["nonadopted_fit_size_uas"]),
                implied_size,
            )
            self.assertEqual(
                self.by_key[key]["nonadopted_fit_interval_status"],
                "bayesian_formal_conditional_on_model",
            )
            self.assertEqual(
                self.by_key[key]["nonadopted_fit_adoption_status"],
                "reported_separately_not_adopted",
            )
            self.assertEqual(
                self.by_key[key]["nonadopted_fit_reference_doi"],
                "10.1051/0004-6361/202450098",
            )
        self.assertEqual(self.rendered.count(r"\nodata"), 4)

        mutations = (
            (
                "nonadopted_fit_adoption_status",
                "adopted",
                "Mrk1029",
            ),
            (
                "nonadopted_fit_interval_status",
                "unconditional",
                "Mrk1029",
            ),
            (
                "nonadopted_fit_ratio_msun_per_mpc",
                "14500",
                "Mrk1029",
            ),
            (
                "nonadopted_fit_reference_doi",
                "10.0000/wrong",
                "Mrk1029",
            ),
        )
        for field, value, message in mutations:
            with self.subTest(field=field):
                mutated = dict(self.by_key["Mrk1029"])
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "Mrk1029" else item
                    for item in self.params
                ]
                with self.assertRaisesRegex(ValueError, message):
                    table.validate_inputs(params, self.metadata)

    def test_ngc1194_nonadopted_fit_preserves_the_primary_table_row(self) -> None:
        row = self.by_key["NGC1194"]
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_ratio_msun_per_mpc"]),
            1_310_000.0,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_ratio_unc_minus_msun_per_mpc"]),
            20_900.0,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_ratio_unc_plus_msun_per_mpc"]),
            21_100.0,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_size_uas"]),
            0.13437794321828395,
        )
        self.assertEqual(
            table.optional_float(row["size_uas"]), 0.12533097499265997
        )
        self.assertEqual(
            self.sample_rows["NGC1194"][1],
            r"$0.1253\pm 0.0077^{\dagger}$",
        )
        self.assertEqual(self.sample_rows["NGC1194"][2], "---")
        self.assertNotEqual(self.sample_rows["NGC1194"][1], r"\nodata")

        mutations = (
            (
                "nonadopted_fit_ratio_msun_per_mpc",
                "1310001",
                "reviewed non-adopted fit",
            ),
            ("nonadopted_fit_basis", "none", "non-adoption basis is invalid"),
            (
                "nonadopted_fit_reference_doi",
                "10.0000/wrong",
                "citation does not match",
            ),
            ("treatment", "no_robust_mass", "requires shared_distance"),
            ("mass_msun", "NA", "requires retained primary cells"),
            ("size_uas", "NA", "requires retained primary cells"),
        )
        for field, value, message in mutations:
            with self.subTest(field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "NGC1194" else item
                    for item in self.params
                ]
                with self.assertRaisesRegex(ValueError, message):
                    table.validate_inputs(params, self.metadata)

    def test_ngc2273_nonadopted_fit_preserves_the_primary_table_row(self) -> None:
        row = self.by_key["NGC2273"]
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_ratio_msun_per_mpc"]),
            328_300.0,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_ratio_unc_minus_msun_per_mpc"]),
            5_100.0,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_ratio_unc_plus_msun_per_mpc"]),
            4_500.0,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_size_uas"]),
            0.033676548670658495,
        )
        self.assertEqual(
            table.optional_float(row["size_uas"]), 0.029935384029973856
        )
        self.assertEqual(
            self.sample_rows["NGC2273"][1],
            r"$0.0299\pm 0.0020^{\dagger}$",
        )
        self.assertEqual(self.sample_rows["NGC2273"][2], "---")

        mutations = (
            (
                "nonadopted_fit_ratio_msun_per_mpc",
                "328301",
                "reviewed non-adopted fit",
            ),
            ("nonadopted_fit_basis", "none", "non-adoption basis is invalid"),
            (
                "nonadopted_fit_reference_doi",
                "10.0000/wrong",
                "citation does not match",
            ),
            ("treatment", "no_robust_mass", "requires shared_distance"),
            ("mass_msun", "NA", "requires retained primary cells"),
            ("size_uas", "NA", "requires retained primary cells"),
        )
        for field, value, message in mutations:
            with self.subTest(field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "NGC2273" else item
                    for item in self.params
                ]
                with self.assertRaisesRegex(ValueError, message):
                    table.validate_inputs(params, self.metadata)

    def test_ngc5765b_nonadopted_fit_preserves_the_primary_table_row(self) -> None:
        row = self.by_key["NGC5765b"]
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_ratio_msun_per_mpc"]),
            372_700.0,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_ratio_unc_minus_msun_per_mpc"]),
            1_300.0,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_ratio_unc_plus_msun_per_mpc"]),
            1_300.0,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_size_uas"]),
            0.038231037738514835,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_size_unc_minus_uas"]),
            0.00013335215739219018,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_size_unc_plus_uas"]),
            0.00013335215739219018,
        )
        self.assertEqual(
            table.optional_float(row["size_uas"]), 0.037924638072366623
        )
        self.assertNotEqual(
            table.optional_float(row["size_uas"]),
            table.optional_float(row["nonadopted_fit_size_uas"]),
        )

        mutations = (
            (
                "nonadopted_fit_ratio_msun_per_mpc",
                "372701",
                "reviewed non-adopted fit",
            ),
            ("nonadopted_fit_basis", "none", "non-adoption basis is invalid"),
            (
                "nonadopted_fit_reference_doi",
                "10.0000/wrong",
                "citation does not match",
            ),
            ("treatment", "shared_distance", "requires paired_posterior"),
            ("mass_msun", "NA", "requires retained primary cells"),
            ("size_uas", "NA", "requires retained primary cells"),
        )
        for field, value, message in mutations:
            with self.subTest(field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "NGC5765b" else item
                    for item in self.params
                ]
                with self.assertRaisesRegex(ValueError, message):
                    table.validate_inputs(params, self.metadata)

    def test_ngc6264_nonadopted_fit_preserves_the_primary_table_row(self) -> None:
        row = self.by_key["NGC6264"]
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_ratio_msun_per_mpc"]),
            212_600.0,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_ratio_unc_minus_msun_per_mpc"]),
            700.0,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_ratio_unc_plus_msun_per_mpc"]),
            700.0,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_size_uas"]),
            0.021808206662753567,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_size_unc_minus_uas"]),
            0.00007180500782656395,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_size_unc_plus_uas"]),
            0.00007180500782656395,
        )
        self.assertEqual(table.optional_float(row["size_uas"]), 0.021455952583840227)
        self.assertNotEqual(
            table.optional_float(row["size_uas"]),
            table.optional_float(row["nonadopted_fit_size_uas"]),
        )

        mutations = (
            (
                "nonadopted_fit_ratio_msun_per_mpc",
                "212601",
                "reviewed non-adopted fit",
            ),
            ("nonadopted_fit_basis", "none", "non-adoption basis is invalid"),
            (
                "nonadopted_fit_reference_doi",
                "10.0000/wrong",
                "citation does not match",
            ),
            ("treatment", "shared_distance", "requires paired_posterior"),
            ("mass_msun", "NA", "requires retained primary cells"),
            ("size_uas", "NA", "requires retained primary cells"),
        )
        for field, value, message in mutations:
            with self.subTest(field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "NGC6264" else item
                    for item in self.params
                ]
                with self.assertRaisesRegex(ValueError, message):
                    table.validate_inputs(params, self.metadata)

    def test_ugc3789_nonadopted_fit_preserves_the_primary_table_row(self) -> None:
        row = self.by_key["UGC3789"]
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_ratio_msun_per_mpc"]),
            227_700.0,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_ratio_unc_minus_msun_per_mpc"]),
            900.0,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_ratio_unc_plus_msun_per_mpc"]),
            900.0,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_size_uas"]),
            0.02335714326015516,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_size_unc_minus_uas"]),
            0.00009232072434843937,
        )
        self.assertEqual(
            table.optional_float(row["nonadopted_fit_size_unc_plus_uas"]),
            0.00009232072434843937,
        )
        self.assertEqual(table.optional_float(row["size_uas"]), 0.023705886582872197)
        self.assertNotEqual(
            table.optional_float(row["size_uas"]),
            table.optional_float(row["nonadopted_fit_size_uas"]),
        )

        mutations = (
            (
                "nonadopted_fit_ratio_msun_per_mpc",
                "227701",
                "reviewed non-adopted fit",
            ),
            ("nonadopted_fit_basis", "none", "non-adoption basis is invalid"),
            (
                "nonadopted_fit_reference_doi",
                "10.0000/wrong",
                "citation does not match",
            ),
            ("treatment", "shared_distance", "requires paired_posterior"),
            ("mass_msun", "NA", "requires retained primary cells"),
            ("size_uas", "NA", "requires retained primary cells"),
        )
        for field, value, message in mutations:
            with self.subTest(field=field):
                mutated = dict(row)
                mutated[field] = value
                params = [
                    mutated if item["source_key"] == "UGC3789" else item
                    for item in self.params
                ]
                with self.assertRaisesRegex(ValueError, message):
                    table.validate_inputs(params, self.metadata)

    def test_comparison_rows_and_inclination_qualifiers_are_retained(self) -> None:
        self.assertEqual(
            self.comparison_rows[0].removesuffix(r" \\"),
            r"Sgr~A* & $^{1}$\,$51.8 \pm 2.3$ & --- & "
            r"$^{2}$\,$6.6332 \pm 0.0012$ & "
            r"$^{1}$\,$0.00815 \pm 0.00015$ & $^{3}$\,$\leq 30$",
        )
        self.assertEqual(
            self.comparison_rows[1].removesuffix(r" \\"),
            r"M87* & $^{4}$\,$42 \pm 3^{\dagger}$ & --- & "
            r"$^{5}$\,$9.74 \pm 0.03$ & "
            r"$^{4}$\,$16.8 \pm 0.8$ & $^{6}$\,$\sim 17$",
        )
        self.assertIn(r"\leq 30", self.comparison_rows[0])
        self.assertIn(r"\sim 17", self.comparison_rows[1])

    def test_cited_inclination_precision_and_asymmetry_are_retained(self) -> None:
        self.assertIn(r"\sim90", self.sample_rows["NGC1068"][5])
        self.assertEqual(
            self.sample_rows["CGCG074-064"][5],
            r"$^{15}$\,$90.8\pm1.3$ (at center)",
        )
        self.assertIn(r"72.4\pm0.5", self.sample_rows["NGC5765b"][5])
        self.assertEqual(
            self.sample_rows["UGC6093"][5], r"$^{17}$\,$93.8\pm4.0$"
        )
        self.assertIn(r"91.5\pm0.3$ (at center)", self.sample_rows["NGC6323"][5])

    def test_reference_numbers_are_sequential(self) -> None:
        metadata_by_key = {row["source_key"]: row for row in self.metadata}
        display_metadata = [
            metadata_by_key[source_key] for source_key in table.DISPLAY_SOURCE_KEYS
        ]
        keys = table.reference_order(display_metadata)
        expected_keys = [
            "EventHorizonTelescopeCollaboration2022d",
            "GRAVITYCollaboration2022",
            "EventHorizonTelescopeCollaboration2022f",
            "EventHorizonTelescopeCollaboration2019f",
            "Simon2024",
            "EventHorizonTelescopeCollaboration2021b",
            "Reid2019",
            "Kuo2010",
            "Gallimore2024",
            "Markham2026",
            "Kondratko2008",
            "Greenhill2003",
            "Greenhill1997",
            "Pesce2020b",
            "Pesce2020a",
            "Kuo2018Disk",
            "Zhao2018",
            "Gao2016b",
            "Kondratko2005",
            "Nakai2024",
            "Yamauchi2012",
            "Pesce2021J0437",
        ]
        self.assertEqual(keys, expected_keys)

        body = self.rendered.split("\\startdata\n", 1)[1].split("\\enddata", 1)[0]
        cell_reference_tokens = re.findall(r"\$\^\{([^}]*)\}\$\\,", body)
        self.assertTrue(cell_reference_tokens)
        self.assertTrue(all(token.isdigit() for token in cell_reference_tokens))
        cell_reference_numbers = [int(token) for token in cell_reference_tokens]
        self.assertEqual(
            list(dict.fromkeys(cell_reference_numbers)),
            list(range(1, len(expected_keys) + 1)),
        )

        reference_text = self.rendered.split(r"\tablerefs{", 1)[1]
        mappings = re.findall(r"\((\d+)\)~\\citet\{([^}]+)\}", reference_text)
        self.assertEqual(
            mappings,
            [(str(number), key) for number, key in enumerate(expected_keys, start=1)],
        )
        for number, key in enumerate(expected_keys, start=1):
            self.assertEqual(
                self.rendered.count(rf"({number})~\citet{{{key}}}"),
                1,
            )

        displayed_sources = [
            row.split(" & ", 1)[0] for row in self.comparison_rows
        ] + [cells[0] for cells in self.sample_rows.values()]
        self.assertEqual(
            displayed_sources,
            [
                "Sgr~A*",
                "M87*",
                r"\textbf{NGC~4258}",
                r"\textbf{NGC~1194}",
                r"\textbf{NGC~1068}",
                "NGC~3393",
                "NGC~4388",
                r"\textbf{Circinus}",
                r"\textbf{NGC~4945}",
                "NGC~5765b",
                "NGC~2273",
                "CGCG~074--064",
                "UGC~3789",
                "NGC~6264",
                "NGC~2960",
                "UGC~6093",
                "ESO~558--G009",
                r"\textbf{NGC~3079}",
                "NGC~6323",
                "IC~2560",
                "J0437+2456",
                "Mrk~1029",
                "NGC~1320",
            ],
        )

    def test_obsolete_uncertainty_language_is_absent(self) -> None:
        lowered = self.rendered.lower()
        for forbidden in (r"\rho", r"\ddagger", "independent", "correlat"):
            self.assertNotIn(forbidden, lowered)
        self.assertEqual(lowered.count(r"^{\dagger}"), 4)
        self.assertNotIn("not total", lowered)

    def test_cli_requires_and_honors_an_explicit_output(self) -> None:
        with self.assertRaises(SystemExit):
            table.parse_args([])
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "BHS_Params.tex"
            self.assertEqual(table.main(["--output", str(output)]), 0)
            self.assertEqual(output.read_text(encoding="utf-8"), self.rendered)

    def test_checked_in_manuscript_table_is_exact_build(self) -> None:
        if not table.DEFAULT_MANUSCRIPT_OUTPUT.exists():
            self.skipTest("revision-package manuscript is not part of the Zenodo archive")
        self.assertEqual(
            table.DEFAULT_MANUSCRIPT_OUTPUT.read_text(encoding="utf-8"),
            self.rendered,
        )


if __name__ == "__main__":
    unittest.main()
