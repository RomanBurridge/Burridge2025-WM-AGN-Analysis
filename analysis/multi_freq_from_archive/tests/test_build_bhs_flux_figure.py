from __future__ import annotations

import csv
import hashlib
import math
import sys
import tempfile
import unittest
from pathlib import Path


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))

import build_bhs_flux_figure as figure_builder  # noqa: E402


class ExactLogUncertaintyTests(unittest.TestCase):
    def test_linear_bounds_are_transformed_exactly_and_asymmetrically(self) -> None:
        lower, upper = figure_builder.exact_log10_uncertainties(10.0, 2.0, 5.0)
        self.assertAlmostEqual(lower, math.log10(10.0) - math.log10(8.0))
        self.assertAlmostEqual(upper, math.log10(15.0) - math.log10(10.0))
        self.assertNotAlmostEqual(lower, 2.0 / (10.0 * math.log(10.0)))
        self.assertNotAlmostEqual(upper, 5.0 / (10.0 * math.log(10.0)))

    def test_nonpositive_lower_bound_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "remain positive"):
            figure_builder.exact_log10_uncertainties(1.0, 1.0, 0.1)


class PublicationFigureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = figure_builder.load_plot_data()
        cls.by_source = {point.source_key: point for point in cls.data.points}

    def test_publication_tables_join_to_the_expected_sample(self) -> None:
        self.assertEqual(len(self.data.points), 17)
        self.assertEqual(
            set(self.by_source),
            {
                "Circinus",
                "ESO558-G009",
                "IC2560",
                "NGC1068",
                "NGC1194",
                "NGC2273",
                "NGC2960",
                "NGC3079",
                "NGC3393",
                "NGC4258",
                "NGC4388",
                "NGC4945",
                "NGC5765b",
                "NGC6264",
                "NGC6323",
                "UGC3789",
                "UGC6093",
            },
        )
        # CGCG 074-064 has no canonical 200--400 GHz S230 anchor (its
        # available continuum datum is at 178.3 GHz), so neither its adopted
        # paired result nor the Villaraos conditional fit belongs in this plot.
        self.assertNotIn("CGCG074-064", self.by_source)
        self.assertEqual(set(self.data.skipped_without_size), {"Mrk1029", "NGC1320"})
        self.assertEqual(
            {
                point.source_key
                for point in self.data.points
                if point.uncertainty_basis == "missing"
            },
            {"NGC4945", "NGC3079"},
        )
        self.assertEqual(self.by_source["NGC5765b"].display_name, "NGC5765B")

    def test_total_size_bounds_are_used_in_exact_log_space(self) -> None:
        point = self.by_source["NGC4258"]
        self.assertEqual(point.uncertainty_basis, "total")
        expected_lower, expected_upper = figure_builder.exact_log10_uncertainties(
            point.size_uas,
            0.0054150793401046308,
            0.005420935810321683,
        )
        self.assertAlmostEqual(point.size_log_unc_minus, expected_lower)
        self.assertAlmostEqual(point.size_log_unc_plus, expected_upper)
        self.assertNotEqual(point.size_log_unc_minus, point.size_log_unc_plus)

    def test_ngc5765b_keeps_the_paired_primary_and_total_error_bar(self) -> None:
        point = self.by_source["NGC5765b"]
        self.assertEqual(point.size_uas, 0.037924638072366623)
        self.assertNotEqual(point.size_uas, 0.038231037738514835)
        self.assertEqual(point.flux_mjy, 0.39481163015431925)
        self.assertEqual(point.log_size, -1.4210785551404157)
        self.assertEqual(point.log_flux, -0.4036100625944722)
        self.assertEqual(point.uncertainty_basis, "total")
        self.assertEqual(
            point.size_stat_component_status, "paired_chain_quantile_interval"
        )
        self.assertEqual(
            point.size_systematic_status,
            "numeric_one_percent_non_keplerian_allowance",
        )
        self.assertEqual(
            point.size_total_status, "quadrature_statistical_plus_systematic"
        )
        expected_lower, expected_upper = figure_builder.exact_log10_uncertainties(
            point.size_uas,
            0.00039922401893172424,
            0.00040259137756787664,
        )
        self.assertEqual(point.size_log_unc_minus, expected_lower)
        self.assertEqual(point.size_log_unc_plus, expected_upper)

    def test_flux_bounds_are_used_in_exact_log_space(self) -> None:
        point = self.by_source["NGC4258"]
        expected_lower, expected_upper = figure_builder.exact_log10_uncertainties(
            point.flux_mjy,
            1.8200615335911141,
            1.8200615335911141,
        )
        self.assertAlmostEqual(point.flux_log_unc_minus, expected_lower)
        self.assertAlmostEqual(point.flux_log_unc_plus, expected_upper)
        self.assertNotEqual(point.flux_log_unc_minus, point.flux_log_unc_plus)

    def test_incomplete_total_uses_identified_component_only(self) -> None:
        point = self.by_source["NGC1068"]
        self.assertEqual(point.uncertainty_basis, "component_only")
        self.assertEqual(point.size_systematic_status, "unquantified_model")
        self.assertEqual(point.size_total_status, "incomplete_unquantified_systematic")
        self.assertIsNotNone(point.size_log_unc_minus)
        self.assertIsNotNone(point.size_log_unc_plus)

        ic2560 = self.by_source["IC2560"]
        self.assertEqual(ic2560.uncertainty_basis, "component_only")
        self.assertEqual(ic2560.size_systematic_status, "unquantified_model")
        self.assertEqual(
            ic2560.size_total_status, "incomplete_unquantified_systematic"
        )
        self.assertIsNotNone(ic2560.size_log_unc_minus)
        self.assertIsNotNone(ic2560.size_log_unc_plus)

        circinus = self.by_source["Circinus"]
        self.assertEqual(circinus.uncertainty_basis, "component_only")
        self.assertEqual(
            circinus.size_stat_component_status, "angular_radius_only_component"
        )
        expected_lower, expected_upper = figure_builder.exact_log10_uncertainties(
            circinus.size_uas,
            float("0.00078339438650467655"),
            float("0.00078339438650467655"),
        )
        self.assertEqual(circinus.size_log_unc_minus, expected_lower)
        self.assertEqual(circinus.size_log_unc_plus, expected_upper)

    def test_stale_circinus_formal_component_token_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "BHS_params.csv"
            text = figure_builder.DEFAULT_BHS_PATH.read_text(encoding="utf-8")
            text = text.replace(
                "angular_radius_only_component",
                "angular_radius_formal_term",
                1,
            )
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unsupported size statistical"):
                figure_builder.load_plot_data(bhs_path=path)

    def test_semantic_statuses_remain_in_audited_plot_data(self) -> None:
        self.assertEqual(self.by_source["NGC1194"].uncertainty_basis, "component_only")
        self.assertEqual(self.by_source["NGC2273"].uncertainty_basis, "component_only")
        self.assertEqual(self.by_source["NGC2960"].uncertainty_basis, "component_only")
        self.assertEqual(
            self.by_source["NGC1194"].size_total_status,
            "reported_undecomposed_uncertainty_component_only",
        )
        self.assertEqual(self.by_source["NGC3393"].uncertainty_basis, "model_span")
        self.assertEqual(self.by_source["ESO558-G009"].uncertainty_basis, "total")
        self.assertEqual(
            self.by_source["ESO558-G009"].size_total_status,
            "quadrature_mixed_fit_position_plus_systematic",
        )

    def test_ngc1194_uses_the_adopted_primary_size_not_the_comparison_fit(self) -> None:
        point = self.by_source["NGC1194"]
        self.assertEqual(point.size_uas, 0.12533097499265997)
        self.assertNotEqual(point.size_uas, 0.13437794321828395)
        expected_lower, expected_upper = figure_builder.exact_log10_uncertainties(
            point.size_uas,
            0.0077126753841636909,
            0.0077126753841636909,
        )
        self.assertEqual(point.size_log_unc_minus, expected_lower)
        self.assertEqual(point.size_log_unc_plus, expected_upper)

    def test_ngc2273_uses_the_adopted_primary_size_not_the_comparison_fit(self) -> None:
        point = self.by_source["NGC2273"]
        self.assertEqual(point.size_uas, 0.029935384029973856)
        self.assertNotEqual(point.size_uas, 0.033676548670658495)
        expected_lower, expected_upper = figure_builder.exact_log10_uncertainties(
            point.size_uas,
            0.0019956922686649235,
            0.0019956922686649235,
        )
        self.assertEqual(point.size_log_unc_minus, expected_lower)
        self.assertEqual(point.size_log_unc_plus, expected_upper)

    def test_ngc4388_uses_reviewed_primary_component_not_separate_checks(self) -> None:
        point = self.by_source["NGC4388"]
        self.assertEqual(point.size_uas, 0.04535053125888249)
        self.assertNotEqual(point.size_uas, 0.043421514018549316)
        self.assertNotEqual(point.size_uas, 0.04539588179014137)
        self.assertEqual(point.uncertainty_basis, "component_only")
        self.assertEqual(
            point.size_stat_component_status,
            "fit_plus_absolute_position_mixed_component",
        )
        self.assertEqual(point.size_systematic_status, "unquantified_model")
        self.assertEqual(point.size_total_status, "incomplete_unquantified_systematic")
        expected_lower, expected_upper = figure_builder.exact_log10_uncertainties(
            point.size_uas,
            0.0010797745537829164,
            0.0010797745537829164,
        )
        self.assertEqual(point.size_log_unc_minus, expected_lower)
        self.assertEqual(point.size_log_unc_plus, expected_upper)

    def test_ngc4388_primary_figure_contract_fails_closed(self) -> None:
        with figure_builder.DEFAULT_BHS_PATH.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames
            originals = list(reader)
        self.assertIsNotNone(fields)

        mutations = (
            ("size_uas", "0.043421514018549316"),
            ("size_stat_unc_minus_uas", "0.0010"),
            ("size_stat_unc_plus_uas", "0.0010"),
            ("size_stat_component_status", "fixed_distance_formal_fit_interval"),
            ("size_systematic_status", "none_adopted"),
            ("size_total_status", "formal_component_only"),
        )
        for field, value in mutations:
            with self.subTest(field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(row) for row in originals]
                next(row for row in rows if row["source_key"] == "NGC4388")[field] = value
                path = Path(tmpdir) / "BHS_params.csv"
                with path.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=fields)
                    writer.writeheader()
                    writer.writerows(rows)
                with self.assertRaisesRegex(
                    ValueError, "NGC4388 reviewed primary size/component"
                ):
                    figure_builder.load_plot_data(bhs_path=path)

    def test_ngc3079_model_check_does_not_become_a_figure_error_bar(self) -> None:
        point = self.by_source["NGC3079"]
        self.assertEqual(point.size_uas, 0.011858795677384631)
        self.assertNotEqual(point.size_uas, 0.015890786207695407)
        self.assertEqual(point.uncertainty_basis, "missing")
        self.assertEqual(point.size_stat_component_status, "no_formal_uncertainty")
        self.assertEqual(point.size_systematic_status, "unquantified_model")
        self.assertEqual(
            point.size_total_status, "incomplete_unquantified_systematic"
        )
        self.assertIsNone(point.size_log_unc_minus)
        self.assertIsNone(point.size_log_unc_plus)

    def test_ngc4945_model_check_does_not_become_a_figure_error_bar(self) -> None:
        point = self.by_source["NGC4945"]
        self.assertEqual(point.size_uas, 0.03881351774408862)
        self.assertNotEqual(point.size_uas, 0.04812876200266989)
        self.assertEqual(point.uncertainty_basis, "missing")
        self.assertEqual(point.size_stat_component_status, "no_formal_uncertainty")
        self.assertEqual(point.size_systematic_status, "no_formal_uncertainty")
        self.assertEqual(point.size_total_status, "no_formal_uncertainty")
        self.assertIsNone(point.size_log_unc_minus)
        self.assertIsNone(point.size_log_unc_plus)

    def test_marker_legend_is_compact_and_complete(self) -> None:
        figure, axis = figure_builder.plt.subplots()
        try:
            figure_builder._marker_legend(axis)
            labels = [text.get_text() for text in axis.get_legend().get_texts()]
        finally:
            figure_builder.plt.close(figure)
        self.assertEqual(
            labels,
            ["Detection", "Upper limit", "M87* and Sgr A*"],
        )

    def test_publication_points_omit_all_error_bars(self) -> None:
        class RecordingAxis:
            def __init__(self) -> None:
                self.errorbar_calls: list[dict[str, object]] = []

            def errorbar(self, *args, **kwargs) -> None:
                self.errorbar_calls.append(dict(kwargs))

            def text(self, *args, **kwargs) -> None:
                return None

        axis = RecordingAxis()
        figure_builder._plot_source_point(axis, self.by_source["NGC4258"])
        self.assertEqual(len(axis.errorbar_calls), 1)
        self.assertNotIn("xerr", axis.errorbar_calls[0])
        self.assertIsNone(axis.errorbar_calls[0]["yerr"])

    def test_comparison_points_omit_all_error_bars(self) -> None:
        class RecordingAxis:
            def __init__(self) -> None:
                self.errorbar_calls: list[dict[str, object]] = []

            def errorbar(self, *args, **kwargs) -> None:
                self.errorbar_calls.append(dict(kwargs))

            def text(self, *args, **kwargs) -> None:
                return None

        axis = RecordingAxis()
        figure_builder._plot_reference(axis, figure_builder.reference_points()[0])
        self.assertEqual(len(axis.errorbar_calls), 1)
        self.assertNotIn("xerr", axis.errorbar_calls[0])
        self.assertIsNone(axis.errorbar_calls[0]["yerr"])

    def test_comparison_sizes_match_manuscript_table(self) -> None:
        by_label = {point.label: point for point in figure_builder.reference_points()}
        self.assertEqual(by_label["M87*"].size_uas, 42.0)
        self.assertEqual(by_label["Sgr A*"].size_uas, 51.8)
        self.assertNotEqual(
            by_label["M87*"].flux_log_unc_minus,
            by_label["M87*"].flux_log_unc_plus,
        )

    def test_generated_pdf_and_png_are_byte_reproducible(self) -> None:
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            result_one = figure_builder.build_figure(output_dir=Path(first))
            result_two = figure_builder.build_figure(output_dir=Path(second))

            self.assertTrue(result_one.pdf_path.read_bytes().startswith(b"%PDF"))
            self.assertTrue(result_one.png_path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))
            for first_path, second_path in (
                (result_one.pdf_path, result_two.pdf_path),
                (result_one.png_path, result_two.png_path),
            ):
                first_hash = hashlib.sha256(first_path.read_bytes()).digest()
                second_hash = hashlib.sha256(second_path.read_bytes()).digest()
                self.assertEqual(first_hash, second_hash)

    def test_checked_in_publication_artifacts_match_the_builder(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            generated = figure_builder.build_figure(output_dir=Path(tmpdir))
            for generated_path in (generated.pdf_path, generated.png_path):
                checked_in = ANALYSIS_DIR / generated_path.name
                self.assertTrue(checked_in.exists(), checked_in)
                self.assertEqual(
                    generated_path.read_bytes(),
                    checked_in.read_bytes(),
                    f"stale checked-in figure: {checked_in}",
                )

    def test_manuscript_figures_are_exact_repository_copies(self) -> None:
        revision_dir = ANALYSIS_DIR.parents[2] / "revision_package"
        if not revision_dir.is_dir():
            self.skipTest("adjacent manuscript revision package is not present")
        for variant in ("working_clean", "working_redline"):
            manuscript_dir = revision_dir / variant
            for filename in ("BHS_Flux_Graph.pdf", "BHS_Flux_Graph.png"):
                with self.subTest(variant=variant, filename=filename):
                    repository_copy = ANALYSIS_DIR / filename
                    manuscript_copy = manuscript_dir / filename
                    self.assertTrue(manuscript_copy.exists(), manuscript_copy)
                    self.assertEqual(
                        repository_copy.read_bytes(),
                        manuscript_copy.read_bytes(),
                        f"stale manuscript figure: {manuscript_copy}",
                    )


if __name__ == "__main__":
    unittest.main()
