from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_DIR))

from bhs_params_utils import (
    get_bhs_plot_uncertainties,
    get_bhs_row,
    load_bhs_params,
    require_final_size_uncertainties,
)


class BhsParamsUtilsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.table = load_bhs_params()

    def test_all_21_sources_and_correct_j0437_alias(self) -> None:
        self.assertEqual(
            len({row["source_key"] for row in self.table.values()}),
            21,
        )
        self.assertEqual(get_bhs_row("WISEAJ043703", self.table)["source_label"], "J0437+2456")

    def test_paired_total_errors_are_transformed_exactly_in_log_space(self) -> None:
        result = get_bhs_plot_uncertainties("NGC4258", self.table)
        self.assertIsNotNone(result)
        self.assertTrue(result["complete"])
        self.assertEqual(result["component"], "total")
        size = result["size_uas"]
        minus = result["linear_minus_uas"]
        plus = result["linear_plus_uas"]
        self.assertEqual(result["log10_minus"], math.log10(size) - math.log10(size - minus))
        self.assertEqual(result["log10_plus"], math.log10(size + plus) - math.log10(size))

    def test_ngc4258_model_check_does_not_replace_the_plotted_error(self) -> None:
        row = get_bhs_row("NGC4258", self.table)
        self.assertEqual(row["size_model_check_lower_uas"], 0.53679114913867576)
        self.assertEqual(row["size_model_check_upper_uas"], 0.54888539817298243)
        self.assertEqual(row["size_model_check_probability_status"], "non_probabilistic")
        self.assertEqual(
            row["size_model_check_adoption_status"],
            "reported_separately_not_adopted",
        )
        result = get_bhs_plot_uncertainties("NGC4258", self.table)
        self.assertEqual(result["linear_minus_uas"], 0.0054150793401046308)
        self.assertEqual(result["linear_plus_uas"], 0.005420935810321683)
        self.assertNotEqual(
            result["linear_minus_uas"],
            result["size_uas"] - row["size_model_check_lower_uas"],
        )

        mutated_table = {key: dict(item) for key, item in self.table.items()}
        for mutated_row in mutated_table.values():
            if mutated_row["source_key"] == "NGC4258":
                mutated_row["size_model_check_lower_uas"] = 0.1
                mutated_row["size_model_check_upper_uas"] = 9.9
        mutated = get_bhs_plot_uncertainties("NGC4258", mutated_table)
        for field in (
            "size_uas",
            "linear_minus_uas",
            "linear_plus_uas",
            "log10_minus",
            "log10_plus",
            "component",
            "complete",
        ):
            self.assertEqual(mutated[field], result[field], field)

    def test_unquantified_model_uses_visible_component_only_interval(self) -> None:
        result = get_bhs_plot_uncertainties("NGC1068", self.table)
        self.assertFalse(result["complete"])
        self.assertEqual(result["component"], "component_only")

        circinus = get_bhs_plot_uncertainties("Circinus", self.table)
        self.assertFalse(circinus["complete"])
        self.assertEqual(circinus["component"], "component_only")
        self.assertEqual(
            circinus["row"]["size_stat_component_status"],
            "angular_radius_only_component",
        )
        self.assertEqual(
            circinus["linear_minus_uas"], float("0.00078339438650467655")
        )
        self.assertEqual(
            circinus["row"]["size_model_check_upper_uas"],
            0.04276549955929029,
        )
        self.assertEqual(circinus["row"]["mass_display_unc_minus_msun"], 300_000.0)
        self.assertEqual(
            circinus["row"]["mass_display_interval_status"],
            "source_reported_total",
        )
        self.assertEqual(
            circinus["row"]["distance_stat_component_status"],
            "reported_legacy_multimethod",
        )
        self.assertNotEqual(
            circinus["linear_plus_uas"],
            circinus["row"]["size_model_check_upper_uas"]
            - circinus["size_uas"],
        )

        ngc4388 = get_bhs_plot_uncertainties("NGC4388", self.table)
        self.assertFalse(ngc4388["complete"])
        self.assertEqual(ngc4388["component"], "component_only")
        self.assertEqual(ngc4388["size_uas"], 0.04535053125888249)
        self.assertEqual(ngc4388["linear_minus_uas"], 0.0010797745537829164)
        self.assertEqual(
            ngc4388["row"]["size_stat_component_status"],
            "fit_plus_absolute_position_mixed_component",
        )
        self.assertEqual(
            ngc4388["row"]["distance_stat_component_status"], "not_available"
        )
        self.assertIsNone(ngc4388["row"]["distance_stat_unc_minus_mpc"])
        self.assertEqual(
            ngc4388["row"]["distance_systematic_unc_minus_mpc"], 2.09
        )
        self.assertEqual(
            ngc4388["row"]["size_model_check_upper_uas"],
            0.04539588179014137,
        )
        self.assertNotEqual(
            ngc4388["size_uas"], ngc4388["row"]["nonadopted_fit_size_uas"]
        )
        self.assertNotEqual(
            ngc4388["linear_plus_uas"],
            ngc4388["row"]["size_model_check_upper_uas"] - ngc4388["size_uas"],
        )

    def test_reported_uncertainty_is_not_misreported_as_a_complete_total(self) -> None:
        result = get_bhs_plot_uncertainties("NGC 1194", self.table)
        self.assertFalse(result["complete"])
        self.assertEqual(result["component"], "component_only")
        self.assertEqual(
            result["row"]["size_total_status"],
            "reported_undecomposed_uncertainty_component_only",
        )

    def test_model_envelope_is_not_misreported_as_probabilistic_total(self) -> None:
        result = get_bhs_plot_uncertainties("NGC 3393", self.table)
        self.assertFalse(result["complete"])
        self.assertEqual(result["component"], "model_span")
        self.assertGreater(result["linear_minus_uas"], 0)
        self.assertGreater(result["linear_plus_uas"], 0)

    def test_no_formal_error_and_unquantified_model_are_explicit(self) -> None:
        ngc4945_result = get_bhs_plot_uncertainties("NGC4945", self.table)
        self.assertFalse(ngc4945_result["complete"])
        self.assertEqual(ngc4945_result["component"], "none")
        for field in (
            "linear_minus_uas",
            "linear_plus_uas",
            "log10_minus",
            "log10_plus",
        ):
            self.assertIsNone(ngc4945_result[field], field)

        result = get_bhs_plot_uncertainties("NGC3079", self.table)
        self.assertFalse(result["complete"])
        self.assertEqual(result["component"], "none")
        self.assertEqual(
            result["row"]["size_stat_component_status"],
            "no_formal_uncertainty",
        )
        self.assertEqual(
            result["row"]["size_systematic_status"],
            "unquantified_model",
        )
        self.assertEqual(
            result["row"]["size_total_status"],
            "incomplete_unquantified_systematic",
        )
        self.assertEqual(
            result["row"]["distance_stat_component_status"],
            "reported_legacy_tully_fisher",
        )
        for field in (
            "linear_minus_uas",
            "linear_plus_uas",
            "log10_minus",
            "log10_plus",
        ):
            self.assertIsNone(result[field], field)

        ngc3079 = get_bhs_row("NGC3079", self.table)
        result = get_bhs_plot_uncertainties("NGC3079", self.table)
        self.assertEqual(result["size_uas"], 0.011858795677384631)
        self.assertEqual(ngc3079["size_model_check_lower_uas"], result["size_uas"])
        self.assertEqual(
            ngc3079["size_model_check_upper_uas"], 0.015890786207695407
        )
        self.assertNotEqual(ngc3079["size_model_check_upper_uas"], result["size_uas"])

        ngc4945 = get_bhs_row("NGC4945", self.table)
        result = get_bhs_plot_uncertainties("NGC4945", self.table)
        self.assertEqual(result["size_uas"], 0.03881351774408862)
        self.assertEqual(ngc4945["size_model_check_lower_uas"], result["size_uas"])
        self.assertEqual(
            ngc4945["size_model_check_upper_uas"], 0.04812876200266989
        )
        self.assertNotEqual(ngc4945["size_model_check_upper_uas"], result["size_uas"])

    def test_ic2560_formal_result_is_visible_and_gate_complete(self) -> None:
        row = get_bhs_row("IC 2560", self.table)
        self.assertEqual(row["size_systematic_status"], "unquantified_model")
        self.assertEqual(row["size_total_status"], "incomplete_unquantified_systematic")
        result = get_bhs_plot_uncertainties("IC 2560", self.table)
        self.assertFalse(result["complete"])
        self.assertEqual(result["component"], "component_only")
        self.assertIsNone(require_final_size_uncertainties(self.table))

    def test_nonadopted_ratio_fits_do_not_replace_primary_plot_values(self) -> None:
        expected = {
            "NGC1194": (1_310_000.0, 0.13437794321828395),
            "NGC2273": (328_300.0, 0.033676548670658495),
            "NGC4388": (423_300.0, 0.043421514018549316),
            "CGCG074-064": (274_900.0, 0.02819885235931776),
            "NGC6264": (212_600.0, 0.021808206662753567),
            "UGC3789": (227_700.0, 0.02335714326015516),
            "ESO558-G009": (159_400.0, 0.016351026067934704),
            "NGC2960": (173_800.0, 0.017828157657509735),
            "NGC5765b": (372_700.0, 0.038231037738514835),
            "NGC6323": (91_600.0, 0.00939619816701894),
            "J0437": (42_500.0, 0.004359589760898525),
            "Mrk1029": (14_400.0, 0.0014771315895750299),
            "NGC1320": (149_700.0, 0.015356013816623746),
        }
        retained_primary = {
            "NGC1194": (
                0.12533097499265997,
                0.0077126753841636909,
                0.0077126753841636909,
            ),
            "NGC2273": (
                0.029935384029973856,
                0.0019956922686649235,
                0.0019956922686649235,
            ),
            "NGC4388": (
                0.04535053125888249,
                0.0010797745537829164,
                0.0010797745537829164,
            ),
            "CGCG074-064": (
                0.02834708992683376,
                0.00033735712639871657,
                0.000349670517390297,
            ),
            "NGC5765b": (
                0.037924638072366623,
                0.00039922401893172424,
                0.00040259137756787664,
            ),
            "NGC6264": (
                0.021455952583840227,
                0.00022076184028948517,
                0.0002220423160751122,
            ),
            "UGC3789": (
                0.023705886582872197,
                0.00024954339475930794,
                0.00025156070016598445,
            ),
            "ESO558-G009": (
                0.01620665338623987,
                0.001041900020436918,
                0.001041900020436918,
            ),
            "NGC2960": (
                0.017574751618808785,
                8.524535950917762e-05,
                9.945291942737387e-05,
            ),
            "NGC6323": (
                0.00954948892245025,
                0.00011076350422409036,
                0.00011348459749234575,
            ),
            "J0437": (
                0.004191067803754549,
                0.0002930816645982202,
                0.0002930816645982202,
            ),
        }
        for key, (ratio, implied_size) in expected.items():
            row = get_bhs_row(key, self.table)
            self.assertEqual(row["nonadopted_fit_ratio_msun_per_mpc"], ratio)
            self.assertEqual(row["nonadopted_fit_size_uas"], implied_size)
            self.assertEqual(
                row["nonadopted_fit_adoption_status"],
                "reported_separately_not_adopted",
            )
            if key in retained_primary:
                primary_size, primary_minus, primary_plus = retained_primary[key]
                self.assertEqual(row["size_uas"], primary_size)
                self.assertNotEqual(row["size_uas"], implied_size)
                plot = get_bhs_plot_uncertainties(key, self.table)
                self.assertEqual(plot["size_uas"], primary_size)
                self.assertEqual(plot["linear_minus_uas"], primary_minus)
                self.assertEqual(plot["linear_plus_uas"], primary_plus)
                if key in {
                    "CGCG074-064",
                    "ESO558-G009",
                    "NGC5765b",
                    "NGC6264",
                    "UGC3789",
                    "NGC6323",
                }:
                    self.assertTrue(plot["complete"])
                    self.assertEqual(plot["component"], "total")
                    self.assertEqual(
                        plot["log10_minus"],
                        math.log10(primary_size)
                        - math.log10(primary_size - primary_minus),
                    )
                    self.assertEqual(
                        plot["log10_plus"],
                        math.log10(primary_size + primary_plus)
                        - math.log10(primary_size),
                    )
                else:
                    self.assertFalse(plot["complete"])
                    self.assertEqual(plot["component"], "component_only")
            else:
                self.assertIsNone(row["size_uas"])
                self.assertIsNone(get_bhs_plot_uncertainties(key, self.table))

    def test_mrk1029_standalone_diagnostics_remain_unplotted_without_a_primary(
        self,
    ) -> None:
        row = get_bhs_row("MRK1029", self.table)
        self.assertEqual(row["source_key"], "Mrk1029")
        self.assertEqual(row["source_order"], 20.0)
        self.assertEqual(row["treatment"], "no_robust_mass")
        self.assertEqual(row["distance_mpc"], 119.32882273342354)
        self.assertEqual(
            row["distance_stat_unc_minus_mpc"], 3.4105504657071966
        )
        self.assertEqual(
            row["distance_stat_unc_plus_mpc"], 3.4105504657071966
        )
        self.assertEqual(
            row["distance_systematic_unc_minus_mpc"], 4.844201193508398
        )
        self.assertEqual(row["distance_unc_minus_mpc"], 5.924368293946939)
        self.assertIsNone(row["size_uas"])
        self.assertIsNone(row["size_unc_adopted_uas"])
        self.assertIsNone(row["mass_msun"])
        self.assertIsNone(row["size_stat_unc_minus_uas"])
        self.assertIsNone(row["size_systematic_unc_minus_uas"])
        self.assertIsNone(row["size_total_unc_minus_uas"])
        self.assertEqual(row["adopted_uncertainty_model"], "not_available")
        self.assertEqual(row["size_stat_component_status"], "not_available")
        self.assertEqual(row["size_systematic_status"], "not_applicable")
        self.assertEqual(row["size_total_status"], "not_available")
        self.assertEqual(row["ratio_provenance"], "no_robust_primary_mass")
        self.assertEqual(
            row["uncertainty_status"],
            "no robust primary black-hole mass; separate non-probabilistic "
            "standalone disk-model diagnostic not combined; conditional model-"
            "dependent literature ratio reported separately and not adopted",
        )

        self.assertEqual(
            row["size_model_check_key"],
            "deepshikha_2024_mrk1029_mestel_disk_mass_model",
        )
        self.assertEqual(row["size_model_check_type"], "disk_mass_model_choice")
        self.assertEqual(
            row["size_model_check_data_relation"],
            "reanalysis_of_published_gao2017_vlbi_data_data_selection_"
            "unreported",
        )
        self.assertEqual(
            row["size_model_check_range_basis"],
            "reported_symmetric_mass_error_span",
        )
        self.assertEqual(row["size_model_check_accepted_configuration_count"], 1.0)
        self.assertEqual(row["size_model_check_excluded_configuration_ids"], "NA")
        self.assertEqual(
            row["size_model_check_lower_uas"], 0.0003317155455200028
        )
        self.assertEqual(
            row["size_model_check_upper_uas"], 0.0005906154834868343
        )
        self.assertEqual(
            row["size_model_check_probability_status"], "non_probabilistic"
        )
        self.assertEqual(
            row["size_model_check_adoption_status"],
            "reported_separately_not_adopted",
        )
        self.assertEqual(row["size_model_check_reference_key"], "Deepshikha2024")
        self.assertEqual(
            row["size_model_check_reference_doi"],
            "10.1017/S1743921323002491",
        )

        self.assertEqual(
            row["nonadopted_fit_key"],
            "gonzalez_juarez_2024_mrk1029_gr_ratio_fit",
        )
        self.assertEqual(row["nonadopted_fit_ratio_msun_per_mpc"], 14_400.0)
        self.assertEqual(
            row["nonadopted_fit_ratio_unc_minus_msun_per_mpc"], 1_100.0
        )
        self.assertEqual(
            row["nonadopted_fit_ratio_unc_plus_msun_per_mpc"], 1_100.0
        )
        self.assertEqual(row["nonadopted_fit_size_uas"], 0.0014771315895750299)
        self.assertEqual(
            row["nonadopted_fit_size_unc_minus_uas"], 0.00011283644087031476
        )
        self.assertEqual(
            row["nonadopted_fit_size_unc_plus_uas"], 0.00011283644087031476
        )
        self.assertEqual(
            row["nonadopted_fit_basis"],
            "restricted_same_data_fit_conditional_on_fixed_gao2017_inclination_"
            "vertical_center_and_selected_scatter_angle_does_not_resolve_primary_"
            "complex_geometry_and_dynamics",
        )
        self.assertEqual(
            row["nonadopted_fit_reference_doi"],
            "10.1051/0004-6361/202450098",
        )

        self.assertIsNone(get_bhs_plot_uncertainties("Mrk 1029", self.table))
        self.assertIsNone(require_final_size_uncertainties(self.table))

        mutated_table = {key: dict(item) for key, item in self.table.items()}
        mutated = get_bhs_row("Mrk1029", mutated_table)
        mutated["size_model_check_lower_uas"] = 0.1
        mutated["size_model_check_upper_uas"] = 9.9
        mutated["nonadopted_fit_ratio_msun_per_mpc"] = 999_999.0
        mutated["nonadopted_fit_size_uas"] = 8.8
        self.assertIsNone(get_bhs_plot_uncertainties("Mrk1029", mutated_table))
        self.assertIsNone(require_final_size_uncertainties(mutated_table))
        for field in (
            "size_uas",
            "size_unc_adopted_uas",
            "mass_msun",
            "size_stat_unc_minus_uas",
            "size_systematic_unc_minus_uas",
            "size_total_unc_minus_uas",
            "adopted_uncertainty_model",
            "size_total_status",
        ):
            self.assertEqual(mutated[field], row[field], field)

    def test_ngc1320_standalone_diagnostics_remain_unplotted_without_a_primary(
        self,
    ) -> None:
        row = get_bhs_row("NGC 1320", self.table)
        self.assertEqual(row["source_key"], "NGC1320")
        self.assertEqual(row["source_order"], 21.0)
        self.assertEqual(row["treatment"], "no_robust_mass")
        self.assertEqual(row["distance_mpc"], 33.783491204330176)
        self.assertEqual(
            row["distance_stat_unc_minus_mpc"], 3.389871133734789
        )
        self.assertEqual(
            row["distance_stat_unc_plus_mpc"], 3.389871133734789
        )
        self.assertEqual(
            row["distance_systematic_unc_minus_mpc"], 1.3714543114071789
        )
        self.assertEqual(
            row["distance_systematic_unc_plus_mpc"], 1.3714543114071789
        )
        self.assertEqual(row["distance_unc_minus_mpc"], 3.656790017434105)
        self.assertEqual(row["distance_unc_plus_mpc"], 3.656790017434105)
        self.assertIsNone(row["size_uas"])
        self.assertIsNone(row["size_unc_adopted_uas"])
        self.assertIsNone(row["mass_msun"])
        self.assertIsNone(row["mass_fit_unc_minus_msun"])
        self.assertIsNone(row["mass_shared_distance_unc_minus_msun"])
        self.assertIsNone(row["mass_total_unc_minus_msun"])
        self.assertIsNone(row["size_stat_unc_minus_uas"])
        self.assertIsNone(row["size_systematic_unc_minus_uas"])
        self.assertIsNone(row["size_total_unc_minus_uas"])
        self.assertEqual(row["adopted_uncertainty_model"], "not_available")
        self.assertEqual(row["size_stat_component_status"], "not_available")
        self.assertEqual(row["size_systematic_status"], "not_applicable")
        self.assertEqual(row["size_total_status"], "not_available")
        self.assertEqual(row["ratio_provenance"], "no_robust_primary_mass")
        self.assertEqual(
            row["uncertainty_status"],
            "no robust primary black-hole mass; separate non-probabilistic "
            "standalone disk-model diagnostic not combined; conditional model-"
            "dependent literature ratio reported separately and not adopted",
        )

        self.assertEqual(
            row["size_model_check_key"],
            "deepshikha_2024_ngc1320_mestel_disk_mass_model",
        )
        self.assertEqual(row["size_model_check_type"], "disk_mass_model_choice")
        self.assertEqual(
            row["size_model_check_data_relation"],
            "reanalysis_of_published_gao2017_vlbi_data_data_selection_"
            "unreported",
        )
        self.assertEqual(
            row["size_model_check_range_basis"],
            "reported_symmetric_mass_error_span",
        )
        self.assertEqual(row["size_model_check_accepted_configuration_count"], 1.0)
        self.assertEqual(row["size_model_check_excluded_configuration_ids"], "NA")
        self.assertEqual(
            row["size_model_check_lower_uas"], 0.01103254026501894
        )
        self.assertEqual(
            row["size_model_check_upper_uas"], 0.014131268281343874
        )
        self.assertAlmostEqual(
            row["size_model_check_lower_uas"] * (4_710_000.0 / 4_130_000.0),
            0.012581904273181408,
            places=17,
        )
        self.assertEqual(
            row["size_model_check_probability_status"], "non_probabilistic"
        )
        self.assertEqual(
            row["size_model_check_adoption_status"],
            "reported_separately_not_adopted",
        )
        self.assertEqual(row["size_model_check_reference_key"], "Deepshikha2024")
        self.assertEqual(
            row["size_model_check_reference_doi"],
            "10.1017/S1743921323002491",
        )
        self.assertNotEqual(
            row["size_model_check_key"],
            "kuo_2020_ngc1320_turbulent_support_bias",
        )

        self.assertEqual(
            row["nonadopted_fit_key"],
            "gonzalez_juarez_2024_ngc1320_gr_ratio_fit",
        )
        self.assertEqual(row["nonadopted_fit_ratio_msun_per_mpc"], 149_700.0)
        self.assertEqual(
            row["nonadopted_fit_ratio_unc_minus_msun_per_mpc"], 6_900.0
        )
        self.assertEqual(
            row["nonadopted_fit_ratio_unc_plus_msun_per_mpc"], 7_100.0
        )
        self.assertEqual(row["nonadopted_fit_size_uas"], 0.015356013816623746)
        self.assertEqual(
            row["nonadopted_fit_size_unc_minus_uas"], 0.0007077922200047018
        )
        self.assertEqual(
            row["nonadopted_fit_size_unc_plus_uas"], 0.0007283079365265772
        )
        self.assertEqual(
            row["nonadopted_fit_basis"],
            "restricted_same_data_fit_conditional_on_fixed_edge_on_inclination_"
            "and_selected_scatter_angle_does_not_resolve_primary_complex_"
            "geometry_and_dynamics",
        )
        self.assertEqual(
            row["nonadopted_fit_reference_key"], "GonzalezJuarez2024"
        )
        self.assertEqual(
            row["nonadopted_fit_reference_doi"],
            "10.1051/0004-6361/202450098",
        )

        self.assertIsNone(get_bhs_plot_uncertainties("NGC1320", self.table))
        self.assertIsNone(require_final_size_uncertainties(self.table))

        mutated_table = {key: dict(item) for key, item in self.table.items()}
        mutated = get_bhs_row("NGC1320", mutated_table)
        mutated["size_model_check_lower_uas"] = 0.1
        mutated["size_model_check_upper_uas"] = 9.9
        mutated["nonadopted_fit_ratio_msun_per_mpc"] = 999_999.0
        mutated["nonadopted_fit_size_uas"] = 8.8
        self.assertIsNone(get_bhs_plot_uncertainties("NGC 1320", mutated_table))
        self.assertIsNone(require_final_size_uncertainties(mutated_table))
        for field in (
            "size_uas",
            "size_unc_adopted_uas",
            "mass_msun",
            "size_stat_unc_minus_uas",
            "size_systematic_unc_minus_uas",
            "size_total_unc_minus_uas",
            "adopted_uncertainty_model",
            "size_total_status",
        ):
            self.assertEqual(mutated[field], row[field], field)

    def test_j0437_sidecars_do_not_replace_the_formal_only_primary_plot(self) -> None:
        row = get_bhs_row("J0437+2456", self.table)
        self.assertEqual(row["size_uas"], 0.004191067803754549)
        self.assertEqual(
            row["size_stat_unc_minus_uas"], 0.0002930816645982202
        )
        self.assertEqual(
            row["size_stat_unc_plus_uas"], 0.0002930816645982202
        )
        self.assertEqual(row["size_systematic_status"], "unquantified_model")
        self.assertEqual(
            row["size_total_status"], "incomplete_unquantified_systematic"
        )
        self.assertIsNone(row["size_total_unc_minus_uas"])
        self.assertIsNone(row["size_total_unc_plus_uas"])

        self.assertEqual(
            row["size_model_check_key"],
            "deepshikha_2024_j0437_mestel_disk_mass_model",
        )
        self.assertEqual(row["size_model_check_type"], "disk_mass_model_choice")
        self.assertEqual(
            row["size_model_check_data_relation"],
            "reanalysis_of_published_gao2017_vlbi_data_data_selection_"
            "unreported",
        )
        self.assertEqual(
            row["size_model_check_range_basis"],
            "reported_symmetric_mass_error_span",
        )
        self.assertEqual(row["size_model_check_accepted_configuration_count"], 1.0)
        self.assertEqual(row["size_model_check_excluded_configuration_ids"], "NA")
        self.assertEqual(
            row["size_model_check_lower_uas"], 0.0024341189955645356
        )
        self.assertEqual(
            row["size_model_check_upper_uas"], 0.0031252885868976758
        )
        self.assertAlmostEqual(
            row["size_model_check_lower_uas"] * (1_850_000.0 / 1_620_000.0),
            0.002779703791231106,
            places=18,
        )
        self.assertEqual(
            row["size_model_check_probability_status"], "non_probabilistic"
        )
        self.assertEqual(
            row["size_model_check_adoption_status"],
            "reported_separately_not_adopted",
        )
        self.assertEqual(row["size_model_check_reference_key"], "Deepshikha2024")
        self.assertEqual(
            row["size_model_check_reference_doi"],
            "10.1017/S1743921323002491",
        )

        self.assertEqual(
            row["nonadopted_fit_key"], "villaraos_2022_j0437_gr_ratio_fit"
        )
        self.assertEqual(row["nonadopted_fit_ratio_msun_per_mpc"], 42_500.0)
        self.assertEqual(row["nonadopted_fit_ratio_unc_minus_msun_per_mpc"], 1_200.0)
        self.assertEqual(row["nonadopted_fit_ratio_unc_plus_msun_per_mpc"], 1_000.0)
        self.assertEqual(row["nonadopted_fit_size_uas"], 0.004359589760898525)
        self.assertEqual(
            row["nonadopted_fit_basis"],
            "same_data_reanalysis_conditional_on_fixed_gao2017_inclination_"
            "vertical_center_and_zero_scatter_angle",
        )

        plotted = get_bhs_plot_uncertainties("J0437", self.table)
        self.assertEqual(plotted["size_uas"], row["size_uas"])
        self.assertEqual(
            plotted["linear_minus_uas"], 0.0002930816645982202
        )
        self.assertEqual(
            plotted["linear_plus_uas"], 0.0002930816645982202
        )
        self.assertFalse(plotted["complete"])
        self.assertEqual(plotted["component"], "component_only")

        mutated_table = {key: dict(item) for key, item in self.table.items()}
        for mutated_row in mutated_table.values():
            if mutated_row["source_key"] != "J0437":
                continue
            mutated_row["size_model_check_lower_uas"] = 0.001
            mutated_row["size_model_check_upper_uas"] = 9.9
            mutated_row["nonadopted_fit_ratio_msun_per_mpc"] = 999_999.0
            mutated_row["nonadopted_fit_size_uas"] = 8.8
        mutated = get_bhs_plot_uncertainties("J0437", mutated_table)
        for field in (
            "size_uas",
            "linear_minus_uas",
            "linear_plus_uas",
            "log10_minus",
            "log10_plus",
            "component",
            "complete",
        ):
            self.assertEqual(mutated[field], plotted[field], field)

    def test_ngc2960_sidecars_do_not_replace_the_formal_only_primary_plot(self) -> None:
        row = get_bhs_row("NGC 2960", self.table)
        self.assertEqual(row["size_uas"], 0.017574751618808785)
        self.assertEqual(row["size_systematic_status"], "unquantified_astrometric")
        self.assertEqual(
            row["size_total_status"], "incomplete_unquantified_systematic"
        )
        self.assertIsNone(row["size_total_unc_minus_uas"])
        self.assertIsNone(row["size_total_unc_plus_uas"])

        self.assertEqual(
            row["size_model_check_key"],
            "kuo_2018_ngc2960_physically_allowed_disk_mass_models",
        )
        self.assertEqual(row["size_model_check_type"], "disk_mass_model_choice")
        self.assertEqual(
            row["size_model_check_lower_uas"], 0.017503713819217805
        )
        self.assertEqual(
            row["size_model_check_upper_uas"], 0.017574751618808785
        )
        self.assertEqual(
            row["size_model_check_probability_status"], "non_probabilistic"
        )
        self.assertEqual(
            row["size_model_check_adoption_status"],
            "reported_separately_not_adopted",
        )
        self.assertAlmostEqual(
            100.0
            * (
                row["size_model_check_upper_uas"]
                - row["size_model_check_lower_uas"]
            )
            / row["size_model_check_upper_uas"],
            0.4042037186742026,
            places=15,
        )

        self.assertEqual(
            row["nonadopted_fit_key"], "villaraos_2022_ngc2960_gr_ratio_fit"
        )
        self.assertEqual(row["nonadopted_fit_ratio_msun_per_mpc"], 173_800.0)
        self.assertEqual(row["nonadopted_fit_size_uas"], 0.017828157657509735)
        self.assertEqual(
            row["nonadopted_fit_basis"],
            "older_published_data_reanalysis_conditional_on_fixed_kuo2011_"
            "inclination_vertical_center_and_selected_scatter_angle",
        )
        self.assertNotEqual(row["size_uas"], row["nonadopted_fit_size_uas"])

        plotted = get_bhs_plot_uncertainties("NGC2960", self.table)
        self.assertEqual(plotted["size_uas"], row["size_uas"])
        self.assertEqual(plotted["linear_minus_uas"], 8.524535950917762e-05)
        self.assertEqual(plotted["linear_plus_uas"], 9.945291942737387e-05)
        self.assertFalse(plotted["complete"])
        self.assertEqual(plotted["component"], "component_only")

        mutated_table = {key: dict(item) for key, item in self.table.items()}
        for mutated_row in mutated_table.values():
            if mutated_row["source_key"] != "NGC2960":
                continue
            mutated_row["size_model_check_lower_uas"] = 0.001
            mutated_row["size_model_check_upper_uas"] = 9.9
            mutated_row["nonadopted_fit_ratio_msun_per_mpc"] = 999_999.0
            mutated_row["nonadopted_fit_size_uas"] = 8.8
        mutated = get_bhs_plot_uncertainties("NGC2960", mutated_table)
        for field in (
            "size_uas",
            "linear_minus_uas",
            "linear_plus_uas",
            "log10_minus",
            "log10_plus",
            "component",
            "complete",
        ):
            self.assertEqual(mutated[field], plotted[field], field)

    def test_ngc6323_sidecars_do_not_replace_the_paired_total_plot(self) -> None:
        row = get_bhs_row("NGC 6323", self.table)
        self.assertEqual(row["size_uas"], 0.00954948892245025)
        self.assertEqual(row["size_posterior_q16_uas"], 0.009493370476285018)
        self.assertEqual(row["size_posterior_q50_uas"], 0.00954948892245025)
        self.assertEqual(row["size_posterior_q84_uas"], 0.009610803521038186)
        self.assertEqual(
            row["size_systematic_status"],
            "numeric_one_percent_non_keplerian_allowance",
        )
        self.assertEqual(
            row["size_total_status"], "quadrature_statistical_plus_systematic"
        )

        self.assertEqual(
            row["size_model_check_key"],
            "kuo_2018_ngc6323_physically_allowed_disk_mass_models",
        )
        self.assertEqual(row["size_model_check_type"], "disk_mass_model_choice")
        self.assertEqual(
            row["size_model_check_data_relation"],
            "predecessor_fit_data_alternative_accretion_disk_models",
        )
        self.assertEqual(row["size_model_check_accepted_configuration_count"], 3)
        self.assertEqual(row["size_model_check_excluded_configuration_ids"], "2;4")
        self.assertEqual(
            row["size_model_check_lower_uas"], 0.009459918021333369
        )
        self.assertEqual(
            row["size_model_check_upper_uas"], 0.009536827923945834
        )
        self.assertEqual(
            row["size_model_check_probability_status"], "non_probabilistic"
        )
        self.assertEqual(
            row["size_model_check_adoption_status"],
            "reported_separately_not_adopted",
        )

        self.assertEqual(
            row["nonadopted_fit_key"], "villaraos_2022_ngc6323_gr_ratio_fit"
        )
        self.assertEqual(row["nonadopted_fit_ratio_msun_per_mpc"], 91_600.0)
        self.assertEqual(row["nonadopted_fit_size_uas"], 0.00939619816701894)
        self.assertEqual(
            row["nonadopted_fit_basis"],
            "older_partial_data_reanalysis_conditional_on_fixed_kuo2015_"
            "inclination_and_selected_scatter_angle",
        )
        self.assertAlmostEqual(
            100.0 * (row["nonadopted_fit_size_uas"] / row["size_uas"] - 1.0),
            -1.6052247054911262,
            places=14,
        )

        plotted = get_bhs_plot_uncertainties("NGC6323", self.table)
        self.assertEqual(plotted["size_uas"], row["size_uas"])
        self.assertEqual(
            plotted["linear_minus_uas"], 0.00011076350422409036
        )
        self.assertEqual(
            plotted["linear_plus_uas"], 0.00011348459749234575
        )
        self.assertTrue(plotted["complete"])
        self.assertEqual(plotted["component"], "total")

        mutated_table = {key: dict(item) for key, item in self.table.items()}
        for mutated_row in mutated_table.values():
            if mutated_row["source_key"] != "NGC6323":
                continue
            mutated_row["size_model_check_lower_uas"] = 0.001
            mutated_row["size_model_check_upper_uas"] = 9.9
            mutated_row["nonadopted_fit_ratio_msun_per_mpc"] = 999_999.0
            mutated_row["nonadopted_fit_size_uas"] = 8.8
        mutated = get_bhs_plot_uncertainties("NGC6323", mutated_table)
        for field in (
            "size_uas",
            "linear_minus_uas",
            "linear_plus_uas",
            "log10_minus",
            "log10_plus",
            "component",
            "complete",
        ):
            self.assertEqual(mutated[field], plotted[field], field)

    def test_cgcg_mass_model_systematic_is_typed_but_not_a_size_error(self) -> None:
        row = get_bhs_row("CGCG 074-064", self.table)
        self.assertEqual(row["mass_systematic_unc_minus_msun"], 500_000.0)
        self.assertEqual(row["mass_systematic_unc_plus_msun"], 500_000.0)
        self.assertEqual(
            row["log_mass_systematic_unc_minus"], 0.009067019970327397
        )
        self.assertEqual(
            row["log_mass_systematic_unc_plus"], 0.008881587279234472
        )
        self.assertEqual(
            row["mass_systematic_status"],
            "source_reported_model_selection_systematic",
        )
        self.assertEqual(row["mass_total_unc_minus_msun"], 2_061_552.8128088303)
        self.assertEqual(row["mass_total_unc_plus_msun"], 2_256_102.8345356956)

        plotted = get_bhs_plot_uncertainties("CGCG074-064", self.table)
        self.assertEqual(plotted["size_uas"], 0.02834708992683376)
        self.assertEqual(
            plotted["linear_minus_uas"], 0.00033735712639871657
        )
        self.assertEqual(plotted["linear_plus_uas"], 0.000349670517390297)
        self.assertTrue(plotted["complete"])
        self.assertEqual(plotted["component"], "total")

    def test_eso558_geometry_mass_systematic_and_primary_plot_are_retained(
        self,
    ) -> None:
        row = get_bhs_row("ESO 558-G009", self.table)
        self.assertEqual(row["mass_msun"], 16_792_963.464140728)
        self.assertEqual(
            row["mass_systematic_unc_minus_msun"], 1_008_137.4178931694
        )
        self.assertEqual(
            row["mass_systematic_unc_plus_msun"], 1_008_137.4178931694
        )
        self.assertEqual(
            row["log_mass_systematic_unc_minus"], 0.02688754291121219
        )
        self.assertEqual(
            row["log_mass_systematic_unc_plus"], 0.02531951831781824
        )
        self.assertEqual(
            row["mass_systematic_status"],
            "source_reported_midline_inclination_systematic",
        )
        self.assertEqual(
            row["mass_display_interval_status"],
            "separate_decomposed_and_geometry_systematic_components",
        )
        self.assertEqual(row["mass_display_unc_minus_msun"], 950_221.0191704632)
        self.assertEqual(row["mass_display_unc_plus_msun"], 950_221.0191704632)
        self.assertEqual(row["mass_total_unc_minus_msun"], 1_385_373.9706771097)
        self.assertEqual(row["mass_total_unc_plus_msun"], 1_385_373.9706771097)

        self.assertEqual(
            row["nonadopted_fit_key"],
            "villaraos_2022_eso558g009_gr_ratio_fit",
        )
        self.assertEqual(row["nonadopted_fit_ratio_msun_per_mpc"], 159_400.0)
        self.assertEqual(row["nonadopted_fit_size_uas"], 0.016351026067934704)
        self.assertEqual(
            row["nonadopted_fit_adoption_status"],
            "reported_separately_not_adopted",
        )

        plotted = get_bhs_plot_uncertainties("ESO558", self.table)
        self.assertEqual(plotted["size_uas"], 0.01620665338623987)
        self.assertEqual(plotted["linear_minus_uas"], 0.001041900020436918)
        self.assertEqual(plotted["linear_plus_uas"], 0.001041900020436918)
        self.assertTrue(plotted["complete"])
        self.assertEqual(plotted["component"], "total")
        self.assertEqual(
            plotted["row"]["size_total_status"],
            "quadrature_mixed_fit_position_plus_systematic",
        )
        self.assertNotEqual(plotted["size_uas"], row["nonadopted_fit_size_uas"])

        # Neither the separately recorded Villaraos comparison nor marginal
        # physical-mass bookkeeping is allowed to replace/reconstruct the
        # canonical shadow-size point used downstream.
        mutated_table = {key: dict(item) for key, item in self.table.items()}
        for mutated_row in mutated_table.values():
            if mutated_row["source_key"] != "ESO558-G009":
                continue
            mutated_row["nonadopted_fit_ratio_msun_per_mpc"] = 999_999.0
            mutated_row["nonadopted_fit_size_uas"] = 9.9
            mutated_row["mass_systematic_unc_minus_msun"] = 1.0
            mutated_row["mass_systematic_unc_plus_msun"] = 2.0
            mutated_row["mass_total_unc_minus_msun"] = 3.0
            mutated_row["mass_total_unc_plus_msun"] = 4.0
        mutated = get_bhs_plot_uncertainties("ESO558-G009", mutated_table)
        for field in (
            "size_uas",
            "linear_minus_uas",
            "linear_plus_uas",
            "log10_minus",
            "log10_plus",
            "component",
            "complete",
        ):
            self.assertEqual(mutated[field], plotted[field], field)

    def test_unknown_total_status_fails_closed(self) -> None:
        row = dict(get_bhs_row("NGC 1194", self.table))
        row["size_total_status"] = "unknown_status"
        with self.assertRaisesRegex(ValueError, "unsupported size_total_status"):
            get_bhs_plot_uncertainties("NGC1194", {"NGC1194": row})

    def test_stale_circinus_formal_component_token_fails_closed(self) -> None:
        row = dict(get_bhs_row("Circinus", self.table))
        row["size_stat_component_status"] = "angular_radius_formal_term"
        with self.assertRaisesRegex(ValueError, "unsupported size statistical"):
            get_bhs_plot_uncertainties("Circinus", {"CIRCINUS": row})

    def test_stale_ngc4388_component_or_distance_classification_fails_closed(self) -> None:
        for field, value, message in (
            (
                "size_stat_component_status",
                "fixed_distance_formal_fit_interval",
                "fit plus position",
            ),
            ("distance_stat_component_status", "statistical", "2.09 Mpc systematic"),
            ("distance_stat_unc_minus_mpc", 2.09, "2.09 Mpc systematic"),
            ("distance_systematic_unc_minus_mpc", 2.08, "2.09 Mpc systematic"),
        ):
            with self.subTest(field=field):
                row = dict(get_bhs_row("NGC4388", self.table))
                row[field] = value
                with self.assertRaisesRegex(ValueError, message):
                    get_bhs_plot_uncertainties("NGC4388", {"NGC4388": row})


if __name__ == "__main__":
    unittest.main()
