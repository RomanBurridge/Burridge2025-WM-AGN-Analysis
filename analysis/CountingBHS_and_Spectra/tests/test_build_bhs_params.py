from __future__ import annotations

import csv
import io
import math
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path


MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import build_bhs_params as bhs  # noqa: E402


EXPECTED_CENTRAL_SIZES = {
    "NGC4258": 0.5386030224537178,
    "NGC1194": 0.12533097499265997,
    "NGC1068": 0.11825031050803191,
    "NGC3393": 0.06359872121781378,
    "NGC4388": 0.04535053125888249,
    "Circinus": 0.04151990248474786,
    "NGC5765b": 0.03792463807236662,
    "NGC4945": 0.03881351774408862,
    "NGC2273": 0.029935384029973856,
    "CGCG074-064": 0.02834708992683376,
    "NGC6264": 0.021455952583840227,
    "UGC3789": 0.023705886582872197,
    "UGC6093": 0.017445797174172238,
    "ESO558-G009": 0.01620665338623987,
    "NGC2960": 0.017574751618808785,
    "IC2560": 0.009543265887703845,
    "NGC3079": 0.011858795677384631,
    "NGC6323": 0.00954948892245025,
    "J0437": 0.004191067803754549,
    "Mrk1029": None,
    "NGC1320": None,
}


class BuildBhsParamsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.inputs = bhs.load_source_inputs()
        cls.quantiles = bhs.load_paired_quantiles()
        cls.source_keys = {row["source_key"] for row in cls.inputs}
        cls.source_treatments = {
            row["source_key"]: row["treatment"] for row in cls.inputs
        }
        cls.model_checks = bhs.load_model_checks(known_sources=cls.source_keys)
        cls.audit_only_model_checks = bhs.load_audit_only_model_checks(
            known_sources=cls.source_keys
        )
        cls.nonadopted_fits = bhs.load_nonadopted_ratio_fits(
            source_treatments=cls.source_treatments
        )
        cls.rows = bhs.build_rows()
        cls.by_key = {row["source_key"]: row for row in cls.rows}

    def _write_model_checks(
        self,
        directory: str,
        rows: list[dict[str, str]],
        *,
        fields: list[str] | None = None,
    ) -> Path:
        path = Path(directory) / "model_checks.csv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields or bhs.MODEL_CHECK_FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        return path

    def _write_source_inputs(
        self, directory: str, rows: list[dict[str, str]]
    ) -> Path:
        path = Path(directory) / "source_inputs.csv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=bhs.INPUT_FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        return path

    def _write_audit_only_model_checks(
        self,
        directory: str,
        rows: list[dict[str, str]],
        *,
        fields: list[str] | None = None,
    ) -> Path:
        path = Path(directory) / "audit_only_model_checks.csv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fields or bhs.AUDIT_ONLY_MODEL_CHECK_FIELDS,
            )
            writer.writeheader()
            writer.writerows(rows)
        return path

    def _write_quantiles(
        self, directory: str, rows: list[dict[str, str]]
    ) -> Path:
        path = Path(directory) / "quantiles.csv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=bhs.QUANTILE_FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        return path

    def _write_nonadopted_fits(
        self,
        directory: str,
        rows: list[dict[str, str]],
        *,
        fields: list[str] | None = None,
    ) -> Path:
        path = Path(directory) / "nonadopted_fits.csv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields or bhs.NONADOPTED_FIT_FIELDS
            )
            writer.writeheader()
            writer.writerows(rows)
        return path

    def test_physical_constant_policy(self) -> None:
        self.assertEqual(bhs.NOMINAL_SOLAR_MASS_PARAMETER, 1.3271244e20)
        self.assertEqual(bhs.SPEED_OF_LIGHT, 299_792_458.0)
        self.assertEqual(bhs.ASTRONOMICAL_UNIT_M, 149_597_870_700.0)
        self.assertTrue(
            math.isclose(bhs.PARSEC_M, 3.085677581491367e16, rel_tol=2e-16)
        )

    def test_exact_source_inventory_and_treatments(self) -> None:
        self.assertEqual(len(self.inputs), 21)
        self.assertEqual(len(self.rows), 21)
        self.assertEqual(
            Counter(row["treatment"] for row in self.inputs),
            Counter(bhs.TREATMENT_COUNTS),
        )
        self.assertNotIn("independent", {row["treatment"] for row in self.inputs})
        self.assertEqual([row["source_order"] for row in self.rows], list(range(1, 22)))

    def test_every_central_value_against_frozen_regression_values(self) -> None:
        self.assertEqual(set(self.by_key), set(EXPECTED_CENTRAL_SIZES))
        for key, expected in EXPECTED_CENTRAL_SIZES.items():
            actual = self.by_key[key]["size_uas"]
            if expected is None:
                self.assertIsNone(actual, key)
            else:
                self.assertTrue(
                    math.isclose(actual, expected, rel_tol=2e-15, abs_tol=0.0),
                    f"{key}: {actual} != {expected}",
                )

    def test_paired_theta_g_conversion_and_components(self) -> None:
        paired_keys = {
            row["source_key"]
            for row in self.inputs
            if row["treatment"] == "paired_posterior"
        }
        self.assertEqual(paired_keys, set(self.quantiles))
        for key in paired_keys:
            row = self.by_key[key]
            supplied = self.quantiles[key]
            q16 = bhs.theta_g_nas_to_shadow_uas(supplied["q16_theta_g_nas"])
            q50 = bhs.theta_g_nas_to_shadow_uas(supplied["q50_theta_g_nas"])
            q84 = bhs.theta_g_nas_to_shadow_uas(supplied["q84_theta_g_nas"])
            self.assertEqual(row["size_posterior_q16_uas"], q16)
            self.assertEqual(row["size_posterior_q50_uas"], q50)
            self.assertEqual(row["size_posterior_q84_uas"], q84)
            self.assertEqual(row["size_stat_unc_minus_uas"], q50 - q16)
            self.assertEqual(row["size_stat_unc_plus_uas"], q84 - q50)
            self.assertEqual(row["size_systematic_unc_minus_uas"], 0.01 * q50)
            self.assertEqual(row["size_systematic_unc_plus_uas"], 0.01 * q50)
            self.assertEqual(
                row["size_total_unc_minus_uas"],
                math.hypot(q50 - q16, 0.01 * q50),
            )
            self.assertEqual(
                row["size_total_unc_plus_uas"],
                math.hypot(q84 - q50, 0.01 * q50),
            )

    def test_ngc5765b_primary_and_conditional_fit_are_exact_and_separate(self) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "NGC5765b")
        paired = self.quantiles["NGC5765b"]
        baseline = bhs.build_paired_row(source, paired)
        row = self.by_key["NGC5765b"]

        expected_primary = {
            "size_uas": 0.037924638072366623,
            "size_unc_adopted_uas": 0.00040090769824980047,
            "log_mass": 7.6180480967120925,
            "log_mass_unc": 0.020397035514562505,
            "distance_mpc": 112.2,
            "distance_unc_mpc": 5.25,
            "mass_msun": 41_500_000.0,
            "mass_unc_minus_msun": 1_900_000.0,
            "mass_unc_plus_msun": 2_000_000.0,
            "log_mass_unc_minus": 0.020352910786580414,
            "log_mass_unc_plus": 0.020441160242544592,
            "distance_unc_minus_mpc": 5.1,
            "distance_unc_plus_mpc": 5.4,
            "size_posterior_q16_uas": 0.037799930414221673,
            "size_posterior_q50_uas": 0.037924638072366623,
            "size_posterior_q84_uas": 0.038059738035356998,
            "size_stat_unc_minus_uas": 0.00012470765814494994,
            "size_stat_unc_plus_uas": 0.00013509996299037458,
            "size_systematic_unc_minus_uas": 0.00037924638072366626,
            "size_systematic_unc_plus_uas": 0.00037924638072366626,
            "size_total_unc_minus_uas": 0.00039922401893172424,
            "size_total_unc_plus_uas": 0.00040259137756787664,
        }
        for field, expected in expected_primary.items():
            self.assertEqual(row[field], expected, field)
        self.assertEqual(row["treatment"], "paired_posterior")
        self.assertEqual(
            row["size_stat_component_status"], "paired_chain_quantile_interval"
        )
        self.assertEqual(
            row["size_systematic_status"],
            "numeric_one_percent_non_keplerian_allowance",
        )
        self.assertEqual(
            row["size_total_status"], "quadrature_statistical_plus_systematic"
        )
        self.assertEqual(
            row["legacy_uncertainty_basis"], "symmetrized_quadrature_total"
        )
        self.assertEqual(
            row["ratio_provenance"], "D_Pesce_private_communication_2026-08-17"
        )

        for field in bhs.OUTPUT_FIELDS:
            if field.startswith("nonadopted_fit_"):
                continue
            self.assertEqual(row[field], baseline[field], field)

        self.assertEqual(
            row["nonadopted_fit_key"], "villaraos_2022_ngc5765b_gr_ratio_fit"
        )
        self.assertEqual(row["nonadopted_fit_ratio_msun_per_mpc"], 372_700.0)
        self.assertEqual(
            row["nonadopted_fit_ratio_unc_minus_msun_per_mpc"], 1_300.0
        )
        self.assertEqual(
            row["nonadopted_fit_ratio_unc_plus_msun_per_mpc"], 1_300.0
        )
        self.assertEqual(row["nonadopted_fit_size_uas"], 0.038231037738514835)
        self.assertEqual(
            row["nonadopted_fit_size_unc_minus_uas"], 0.00013335215739219018
        )
        self.assertEqual(
            row["nonadopted_fit_size_unc_plus_uas"], 0.00013335215739219018
        )
        self.assertEqual(
            row["nonadopted_fit_interval_status"],
            "bayesian_formal_conditional_on_model",
        )
        self.assertEqual(
            row["nonadopted_fit_adoption_status"],
            "reported_separately_not_adopted",
        )
        self.assertEqual(
            row["nonadopted_fit_basis"],
            "same_data_reanalysis_conditional_on_fixed_gao2016_inclination_and_"
            "selected_scatter_angle",
        )
        self.assertEqual(row["nonadopted_fit_reference_key"], "Villaraos2022")
        self.assertEqual(
            row["nonadopted_fit_reference_doi"], "10.1093/mnras/stac2973"
        )
        self.assertLess(
            row["nonadopted_fit_size_uas"] - row["size_uas"],
            row["size_systematic_unc_plus_uas"],
        )

    def test_cgcg_mass_systematic_conditional_fit_and_bhs_are_exact_and_separate(
        self,
    ) -> None:
        source = next(
            row for row in self.inputs if row["source_key"] == "CGCG074-064"
        )
        paired = self.quantiles["CGCG074-064"]
        baseline = bhs.build_paired_row(source, paired)
        row = self.by_key["CGCG074-064"]

        self.assertEqual(source["reference_mass_msun"], "24200000")
        self.assertEqual(source["reference_mass_unc_minus_msun"], "2000000")
        self.assertEqual(source["reference_mass_unc_plus_msun"], "2200000")
        self.assertEqual(
            source["reference_mass_systematic_unc_minus_msun"], "500000"
        )
        self.assertEqual(
            source["reference_mass_systematic_unc_plus_msun"], "500000"
        )
        self.assertEqual(
            source["mass_systematic_status"],
            "source_reported_model_selection_systematic",
        )
        self.assertEqual(
            source["mass_display_interval_status"],
            "separate_statistical_and_systematic_components",
        )

        expected_primary = {
            "size_uas": 0.02834708992683376,
            "size_posterior_q16_uas": 0.02816418536155449,
            "size_posterior_q50_uas": 0.02834708992683376,
            "size_posterior_q84_uas": 0.028551818332288398,
            "size_stat_unc_minus_uas": 0.0001829045652792724,
            "size_stat_unc_plus_uas": 0.00020472840545463639,
            "size_systematic_unc_minus_uas": 0.0002834708992683376,
            "size_systematic_unc_plus_uas": 0.0002834708992683376,
            "size_total_unc_minus_uas": 0.00033735712639871657,
            "size_total_unc_plus_uas": 0.000349670517390297,
            "size_unc_adopted_uas": 0.0003435138218945068,
            "mass_msun": 24_200_000.0,
            "mass_unc_minus_msun": 2_000_000.0,
            "mass_unc_plus_msun": 2_200_000.0,
            "log_mass": 7.383815365980431,
            "log_mass_unc_minus": 0.03746239152979263,
            "log_mass_unc_plus": 0.037788560889399754,
            "log_mass_unc": 0.03762547620959619,
            "mass_systematic_unc_minus_msun": 500_000.0,
            "mass_systematic_unc_plus_msun": 500_000.0,
            "log_mass_systematic_unc_minus": 0.009067019970327397,
            "log_mass_systematic_unc_plus": 0.008881587279234472,
            "mass_total_unc_minus_msun": 2_061_552.8128088303,
            "mass_total_unc_plus_msun": 2_256_102.8345356956,
            "mass_display_unc_minus_msun": 2_000_000.0,
            "mass_display_unc_plus_msun": 2_200_000.0,
        }
        for field, expected in expected_primary.items():
            self.assertEqual(row[field], expected, field)
        self.assertEqual(
            row["mass_systematic_status"],
            "source_reported_model_selection_systematic",
        )
        self.assertEqual(
            row["mass_display_interval_status"],
            "separate_statistical_and_systematic_components",
        )
        self.assertEqual(
            row["log_mass_unc"],
            (row["log_mass_unc_minus"] + row["log_mass_unc_plus"]) / 2.0,
        )

        # The source-reported mass model systematic is diagnostic and is not
        # folded into the paired-chain M/D shadow-size fields.
        for field in bhs.OUTPUT_FIELDS:
            if field.startswith("nonadopted_fit_"):
                continue
            self.assertEqual(row[field], baseline[field], field)
        self.assertEqual(
            row["nonadopted_fit_key"],
            "villaraos_2022_cgcg074064_gr_ratio_fit",
        )
        self.assertEqual(row["nonadopted_fit_ratio_msun_per_mpc"], 274_900.0)
        self.assertEqual(
            row["nonadopted_fit_ratio_unc_minus_msun_per_mpc"], 1_500.0
        )
        self.assertEqual(
            row["nonadopted_fit_ratio_unc_plus_msun_per_mpc"], 1_700.0
        )
        self.assertEqual(row["nonadopted_fit_size_uas"], 0.02819885235931776)
        self.assertEqual(
            row["nonadopted_fit_size_unc_minus_uas"],
            0.0001538678739140656,
        )
        self.assertEqual(
            row["nonadopted_fit_size_unc_plus_uas"],
            0.000174383590435941,
        )
        fractional_offset = (
            row["size_uas"] - row["nonadopted_fit_size_uas"]
        ) / row["size_uas"]
        self.assertLess(fractional_offset, 0.01)
        self.assertGreater(fractional_offset, 0.0)

    def test_ugc6093_velocity_terms_position_systematic_and_shadow_are_exact(
        self,
    ) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "UGC6093")
        row = self.by_key["UGC6093"]

        expected_source = {
            "source_order": "13",
            "source_label": "UGC 6093",
            "treatment": "shared_distance",
            "reference_mass_msun": "25800000",
            "reference_mass_unc_minus_msun": "900000",
            "reference_mass_unc_plus_msun": "1100000",
            "reference_mass_systematic_unc_minus_msun": "1290000",
            "reference_mass_systematic_unc_plus_msun": "1290000",
            "reference_distance_mpc": "151.7",
            "adopted_distance_mpc": "143.7235453315291",
            "distance_measurement_unc_minus_mpc": "0.16441136671177264",
            "distance_measurement_unc_plus_mpc": "0.14614343707713126",
            "distance_flow_unc_mpc": "3.3829499323410013",
            "distance_shared_systematic_unc_mpc": "5.8345146954612623",
            "size_stat_fraction_minus_override": "",
            "size_stat_fraction_plus_override": "",
            "size_systematic_fraction_minus": "0.05",
            "size_systematic_fraction_plus": "0.05",
            "systematic_status": "numeric_systematic",
            "mass_reference_doi": "10.3847/1538-4357/aaa95c",
            "mass_uncertainty_reference_doi": "10.3847/1538-4357/aaa95c",
            "distance_reference_doi": "10.3847/1538-4357/aaa95c",
            "uncertainty_note": "absolute_position_systematic",
            "mass_display_interval_status": (
                "separate_decomposed_and_position_systematic_components"
            ),
            "distance_stat_component_status": "statistical",
            "mass_systematic_status": (
                "source_reported_absolute_position_systematic"
            ),
        }
        for field, expected in expected_source.items():
            self.assertEqual(source[field], expected, field)

        # Zhao et al.'s asymmetric v_LSR interval is retained rather than
        # symmetrized before division by the adopted H0.
        self.assertAlmostEqual(
            float(source["adopted_distance_mpc"]) * 73.9, 10_621.17, places=12
        )
        self.assertAlmostEqual(
            float(source["distance_measurement_unc_minus_mpc"]) * 73.9,
            12.15,
            places=13,
        )
        self.assertAlmostEqual(
            float(source["distance_measurement_unc_plus_mpc"]) * 73.9,
            10.80,
            places=13,
        )
        self.assertEqual(float(source["distance_flow_unc_mpc"]) * 73.9, 250.0)
        self.assertEqual(
            float(source["distance_shared_systematic_unc_mpc"])
            * 73.9
            / float(source["adopted_distance_mpc"]),
            3.0,
        )

        expected_output = {
            "distance_mpc": 143.7235453315291,
            "distance_stat_unc_minus_mpc": 3.386942772210658,
            "distance_stat_unc_plus_mpc": 3.386105159165424,
            "distance_systematic_unc_minus_mpc": 5.834514695461262,
            "distance_systematic_unc_plus_mpc": 5.834514695461262,
            "distance_unc_minus_mpc": 6.746328117856665,
            "distance_unc_plus_mpc": 6.7459076394863375,
            "mass_msun": 24_443_424.32138069,
            "mass_fit_unc_minus_msun": 852_677.5926063032,
            "mass_fit_unc_plus_msun": 1_042_161.5020743705,
            "mass_shared_distance_unc_minus_msun": 1_147_364.9666493207,
            "mass_shared_distance_unc_plus_msun": 1_147_293.4548368326,
            "mass_systematic_unc_minus_msun": 1_222_171.2160690345,
            "mass_systematic_unc_plus_msun": 1_222_171.2160690345,
            "log_mass_systematic_unc_minus": 0.02227639471115221,
            "log_mass_systematic_unc_plus": 0.021189299069938092,
            "mass_display_unc_minus_msun": 1_429_512.309715127,
            "mass_display_unc_plus_msun": 1_549_962.214996657,
            "mass_total_unc_minus_msun": 1_880_746.6403039885,
            "mass_total_unc_plus_msun": 1_973_850.386758076,
            "size_uas": 0.017445797174172238,
            "size_stat_unc_minus_uas": 0.0006085743200292641,
            "size_stat_unc_plus_uas": 0.000743813057813545,
            "size_systematic_unc_minus_uas": 0.000872289858708612,
            "size_systematic_unc_plus_uas": 0.000872289858708612,
            "size_total_unc_minus_uas": 0.0010636034508241178,
            "size_total_unc_plus_uas": 0.0011463627098697107,
        }
        for field, expected in expected_output.items():
            self.assertEqual(row[field], expected, field)

        self.assertEqual(
            row["mass_display_interval_status"],
            "separate_decomposed_and_position_systematic_components",
        )
        self.assertEqual(
            row["mass_systematic_status"],
            "source_reported_absolute_position_systematic",
        )
        self.assertEqual(
            row["size_stat_component_status"],
            "fixed_distance_formal_fit_interval",
        )
        self.assertEqual(row["size_systematic_status"], "numeric_systematic")
        self.assertEqual(
            row["size_total_status"], "quadrature_statistical_plus_systematic"
        )

        # The physical-mass display interval is fit+distance, with the
        # source-reported 5% absolute-position term retained separately.  The
        # machine total includes each of those three terms exactly once.
        self.assertEqual(
            row["mass_systematic_unc_minus_msun"], 0.05 * row["mass_msun"]
        )
        self.assertEqual(
            row["mass_systematic_unc_plus_msun"], 0.05 * row["mass_msun"]
        )
        self.assertEqual(
            row["mass_display_unc_minus_msun"],
            bhs.positive_quadrature(
                row["mass_fit_unc_minus_msun"],
                row["mass_shared_distance_unc_minus_msun"],
            ),
        )
        self.assertEqual(
            row["mass_display_unc_plus_msun"],
            bhs.positive_quadrature(
                row["mass_fit_unc_plus_msun"],
                row["mass_shared_distance_unc_plus_msun"],
            ),
        )
        for suffix in ("minus", "plus"):
            self.assertEqual(
                row[f"mass_total_unc_{suffix}_msun"],
                bhs.positive_quadrature(
                    row[f"mass_display_unc_{suffix}_msun"],
                    row[f"mass_systematic_unc_{suffix}_msun"],
                ),
                suffix,
            )

        # Shared H0/distance-scale uncertainty belongs to the physical mass,
        # but cancels from M/D.  The shadow total therefore combines only the
        # formal fit and independent absolute-position terms.
        h0_shadow_term = (
            row["size_uas"]
            * row["distance_systematic_unc_minus_mpc"]
            / row["distance_mpc"]
        )
        for suffix in ("minus", "plus"):
            expected_shadow_total = bhs.positive_quadrature(
                row[f"size_stat_unc_{suffix}_uas"],
                row[f"size_systematic_unc_{suffix}_uas"],
            )
            self.assertEqual(
                row[f"size_total_unc_{suffix}_uas"], expected_shadow_total, suffix
            )
            self.assertNotEqual(
                row[f"size_total_unc_{suffix}_uas"],
                bhs.positive_quadrature(
                    row[f"size_stat_unc_{suffix}_uas"],
                    row[f"size_systematic_unc_{suffix}_uas"],
                    h0_shadow_term,
                ),
                suffix,
            )

        without_h0 = dict(source)
        without_h0["distance_shared_systematic_unc_mpc"] = "0"
        rebuilt = bhs.build_shared_row(without_h0)
        for field in (
            "size_uas",
            "size_stat_unc_minus_uas",
            "size_stat_unc_plus_uas",
            "size_systematic_unc_minus_uas",
            "size_systematic_unc_plus_uas",
            "size_total_unc_minus_uas",
            "size_total_unc_plus_uas",
        ):
            self.assertEqual(rebuilt[field], row[field], field)
        self.assertLess(
            rebuilt["mass_shared_distance_unc_minus_msun"],
            row["mass_shared_distance_unc_minus_msun"],
        )
        self.assertLess(
            rebuilt["mass_total_unc_plus_msun"], row["mass_total_unc_plus_msun"]
        )

    def test_ugc6093_reviewed_inputs_and_outputs_fail_closed(self) -> None:
        source_mutations = (
            ("reference_mass_msun", "25800001"),
            ("reference_mass_unc_minus_msun", "900001"),
            ("reference_mass_unc_plus_msun", "1100001"),
            ("reference_mass_systematic_unc_minus_msun", "1290001"),
            ("reference_mass_systematic_unc_plus_msun", "1290001"),
            ("reference_distance_mpc", "151.8"),
            ("adopted_distance_mpc", "143.8"),
            ("distance_measurement_unc_minus_mpc", "0.14614343707713126"),
            ("distance_measurement_unc_plus_mpc", "0.16441136671177264"),
            ("distance_flow_unc_mpc", "3.4"),
            ("distance_shared_systematic_unc_mpc", "5.9"),
            ("size_systematic_fraction_minus", "0.049"),
            ("size_systematic_fraction_plus", "0.051"),
            ("systematic_status", "none_adopted"),
            ("uncertainty_note", "different_note"),
            ("mass_display_interval_status", "decomposed_components"),
            ("distance_stat_component_status", "not_available"),
            ("mass_systematic_status", "not_available"),
        )
        for field, value in source_mutations:
            with self.subTest(source_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in self.inputs]
                next(item for item in rows if item["source_key"] == "UGC6093")[
                    field
                ] = value
                path = self._write_source_inputs(tmpdir, rows)
                with self.assertRaisesRegex(ValueError, "UGC6093"):
                    bhs.load_source_inputs(path)

        output_mutations = (
            ("distance_stat_unc_minus_mpc", 3.3869),
            ("distance_stat_unc_plus_mpc", 3.3862),
            ("size_stat_unc_minus_uas", 0.00061),
            ("size_systematic_unc_plus_uas", 0.00088),
            ("size_total_unc_plus_uas", 0.00115),
            ("mass_fit_unc_minus_msun", 852_678.0),
            ("mass_shared_distance_unc_plus_msun", 1_147_294.0),
            ("mass_systematic_unc_minus_msun", 1_222_172.0),
            ("log_mass_systematic_unc_plus", 0.0212),
            ("mass_display_unc_minus_msun", 1_429_513.0),
            ("mass_total_unc_plus_msun", 1_973_851.0),
            ("mass_systematic_status", "not_available"),
            ("mass_display_interval_status", "decomposed_components"),
            ("size_stat_component_status", "fit_plus_absolute_position_mixed_component"),
            ("size_systematic_status", "none_adopted"),
            ("size_total_status", "formal_component_only"),
        )
        for field, value in output_mutations:
            with self.subTest(output_field=field):
                rows = [dict(item) for item in self.rows]
                next(item for item in rows if item["source_key"] == "UGC6093")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "UGC6093"):
                    bhs.render_csv(rows)

    def test_eso558_distance_geometry_mass_and_shadow_are_exact(self) -> None:
        source = next(
            row for row in self.inputs if row["source_key"] == "ESO558-G009"
        )
        row = self.by_key["ESO558-G009"]

        expected_source = {
            "source_order": "14",
            "source_label": "ESO 558-G009",
            "treatment": "shared_distance",
            "reference_mass_msun": "17000000",
            "reference_mass_unc_minus_msun": "",
            "reference_mass_unc_plus_msun": "",
            "reference_mass_systematic_unc_minus_msun": "1020566.5093466471",
            "reference_mass_systematic_unc_plus_msun": "1020566.5093466471",
            "reference_distance_mpc": "107.6",
            "adopted_distance_mpc": "106.28958051420837",
            "distance_measurement_unc_minus_mpc": "0.3653585926928281",
            "distance_measurement_unc_plus_mpc": "0.3653585926928281",
            "distance_flow_unc_mpc": "3.3829499323410013",
            "distance_shared_systematic_unc_mpc": "4.314867950509135",
            "size_stat_fraction_minus_override": "0.023",
            "size_stat_fraction_plus_override": "0.023",
            "size_systematic_fraction_minus": "0.06003332407921454",
            "size_systematic_fraction_plus": "0.06003332407921454",
            "systematic_status": "numeric_systematic",
            "mass_reference_doi": "10.3847/1538-4357/834/1/52",
            "mass_uncertainty_reference_doi": "10.3847/1538-4357/834/1/52",
            "distance_reference_doi": "10.3847/1538-4357/834/1/52",
            "uncertainty_note": (
                "fit_plus_absolute_position_separate_from_midline_and_"
                "inclination_systematics"
            ),
            "mass_display_interval_status": (
                "separate_decomposed_and_geometry_systematic_components"
            ),
            "distance_stat_component_status": "statistical",
            "mass_systematic_status": (
                "source_reported_midline_inclination_systematic"
            ),
        }
        for field, expected in expected_source.items():
            self.assertEqual(source[field], expected, field)

        # Gao's Table 1 source-information velocity uncertainty is 27 km/s.
        # It is distinct from the fitted recession velocity and is retained
        # exactly rather than rounded to 0.365 Mpc.
        self.assertEqual(
            float(source["distance_measurement_unc_minus_mpc"]) * 73.9,
            27.0,
        )
        self.assertEqual(
            float(source["distance_measurement_unc_plus_mpc"]) * 73.9,
            27.0,
        )
        self.assertEqual(float(source["distance_flow_unc_mpc"]) * 73.9, 250.0)
        self.assertEqual(
            float(source["distance_shared_systematic_unc_mpc"])
            * 73.9
            / float(source["adopted_distance_mpc"]),
            3.0,
        )
        geometry_fraction = math.hypot(0.052, 0.03)
        self.assertEqual(
            float(source["size_systematic_fraction_minus"]), geometry_fraction
        )
        self.assertEqual(
            float(source["reference_mass_systematic_unc_minus_msun"]),
            17_000_000.0 * geometry_fraction,
        )

        expected_output = {
            "distance_mpc": 106.28958051420837,
            "distance_stat_unc_minus_mpc": 3.4026220986146067,
            "distance_stat_unc_plus_mpc": 3.4026220986146067,
            "distance_systematic_unc_minus_mpc": 4.314867950509135,
            "distance_systematic_unc_plus_mpc": 4.314867950509135,
            "distance_unc_minus_mpc": 5.49508167148691,
            "distance_unc_plus_mpc": 5.49508167148691,
            "mass_msun": 16_792_963.464140728,
            "mass_fit_unc_minus_msun": 386_238.15967523673,
            "mass_fit_unc_plus_msun": 386_238.15967523673,
            "mass_shared_distance_unc_minus_msun": 868_182.0484691216,
            "mass_shared_distance_unc_plus_msun": 868_182.0484691216,
            "mass_systematic_unc_minus_msun": 1_008_137.4178931694,
            "mass_systematic_unc_plus_msun": 1_008_137.4178931694,
            "log_mass_systematic_unc_minus": 0.02688754291121219,
            "log_mass_systematic_unc_plus": 0.02531951831781824,
            "mass_display_unc_minus_msun": 950_221.0191704632,
            "mass_display_unc_plus_msun": 950_221.0191704632,
            "mass_total_unc_minus_msun": 1_385_373.9706771097,
            "mass_total_unc_plus_msun": 1_385_373.9706771097,
            "size_uas": 0.01620665338623987,
            "size_stat_unc_minus_uas": 0.00037275302788351703,
            "size_stat_unc_plus_uas": 0.00037275302788351703,
            "size_systematic_unc_minus_uas": 0.0009729392749756379,
            "size_systematic_unc_plus_uas": 0.0009729392749756379,
            "size_total_unc_minus_uas": 0.001041900020436918,
            "size_total_unc_plus_uas": 0.001041900020436918,
            "nonadopted_fit_ratio_msun_per_mpc": 159_400.0,
            "nonadopted_fit_ratio_unc_minus_msun_per_mpc": 2_600.0,
            "nonadopted_fit_ratio_unc_plus_msun_per_mpc": 2_500.0,
            "nonadopted_fit_size_uas": 0.016351026067934704,
            "nonadopted_fit_size_unc_minus_uas": 0.00026670431478438035,
            "nonadopted_fit_size_unc_plus_uas": 0.00025644645652344266,
        }
        for field, expected in expected_output.items():
            self.assertEqual(row[field], expected, field)

        self.assertEqual(
            row["mass_display_interval_status"],
            "separate_decomposed_and_geometry_systematic_components",
        )
        self.assertEqual(
            row["mass_systematic_status"],
            "source_reported_midline_inclination_systematic",
        )
        self.assertEqual(
            row["size_stat_component_status"],
            "fit_plus_absolute_position_mixed_component",
        )
        self.assertEqual(row["size_systematic_status"], "numeric_systematic")
        self.assertEqual(
            row["size_total_status"],
            "quadrature_mixed_fit_position_plus_systematic",
        )

        # The physical-mass display interval contains fit plus all adopted
        # distance terms; the geometry allowance remains a separately rendered
        # component.  The diagnostic full total contains each term once.
        self.assertTrue(
            math.isclose(
                row["mass_systematic_unc_minus_msun"],
                row["mass_msun"] * geometry_fraction,
                rel_tol=2e-16,
                abs_tol=0.0,
            )
        )
        for suffix in ("minus", "plus"):
            self.assertEqual(
                row[f"mass_display_unc_{suffix}_msun"],
                bhs.positive_quadrature(
                    row[f"mass_fit_unc_{suffix}_msun"],
                    row[f"mass_shared_distance_unc_{suffix}_msun"],
                ),
                suffix,
            )
            self.assertTrue(
                math.isclose(
                    row[f"mass_total_unc_{suffix}_msun"],
                    bhs.positive_quadrature(
                        row[f"mass_display_unc_{suffix}_msun"],
                        row[f"mass_systematic_unc_{suffix}_msun"],
                    ),
                    rel_tol=3e-16,
                    abs_tol=0.0,
                ),
                suffix,
            )
        self.assertEqual(
            math.log10(
                row["mass_msun"]
                / (row["mass_msun"] - row["mass_display_unc_minus_msun"])
            ),
            0.025296981151847565,
        )
        self.assertEqual(
            math.log10(
                (row["mass_msun"] + row["mass_display_unc_plus_msun"])
                / row["mass_msun"]
            ),
            0.02390422500815494,
        )
        self.assertEqual(
            math.log10(
                row["mass_msun"]
                / (row["mass_msun"] - row["mass_total_unc_minus_msun"])
            ),
            0.03739264413720991,
        )
        self.assertEqual(
            math.log10(
                (row["mass_msun"] + row["mass_total_unc_plus_msun"])
                / row["mass_msun"]
            ),
            0.0344268176356242,
        )

        # Distance and H0 terms change the physical-mass budget but cancel
        # from M/D.  Only fit-plus-position and geometry enter the shadow total.
        h0_shadow_term = (
            row["size_uas"]
            * row["distance_systematic_unc_minus_mpc"]
            / row["distance_mpc"]
        )
        for suffix in ("minus", "plus"):
            self.assertEqual(
                row[f"size_total_unc_{suffix}_uas"],
                bhs.positive_quadrature(
                    row[f"size_stat_unc_{suffix}_uas"],
                    row[f"size_systematic_unc_{suffix}_uas"],
                ),
                suffix,
            )
            self.assertNotEqual(
                row[f"size_total_unc_{suffix}_uas"],
                bhs.positive_quadrature(
                    row[f"size_stat_unc_{suffix}_uas"],
                    row[f"size_systematic_unc_{suffix}_uas"],
                    h0_shadow_term,
                ),
                suffix,
            )

        without_h0 = dict(source)
        without_h0["distance_shared_systematic_unc_mpc"] = "0"
        rebuilt = bhs.build_shared_row(without_h0)
        for field in (
            "size_uas",
            "size_stat_unc_minus_uas",
            "size_stat_unc_plus_uas",
            "size_systematic_unc_minus_uas",
            "size_systematic_unc_plus_uas",
            "size_total_unc_minus_uas",
            "size_total_unc_plus_uas",
        ):
            self.assertEqual(rebuilt[field], row[field], field)
        self.assertLess(
            rebuilt["mass_shared_distance_unc_minus_msun"],
            row["mass_shared_distance_unc_minus_msun"],
        )
        self.assertLess(
            rebuilt["mass_total_unc_plus_msun"], row["mass_total_unc_plus_msun"]
        )

    def test_eso558_reviewed_inputs_outputs_and_conditional_fit_fail_closed(
        self,
    ) -> None:
        source_mutations = (
            ("source_label", "ESO 558 G009"),
            ("reference_mass_msun", "17000001"),
            ("reference_mass_unc_minus_msun", "391000"),
            ("reference_mass_unc_plus_msun", "391000"),
            ("reference_mass_systematic_unc_minus_msun", "1020567"),
            ("reference_mass_systematic_unc_plus_msun", "1020567"),
            ("reference_distance_mpc", "107.7"),
            ("adopted_distance_mpc", "106.3"),
            ("distance_measurement_unc_minus_mpc", "0.365"),
            ("distance_measurement_unc_plus_mpc", "0.365"),
            ("distance_flow_unc_mpc", "3.4"),
            ("distance_shared_systematic_unc_mpc", "4.3"),
            ("size_stat_fraction_minus_override", "0.022"),
            ("size_stat_fraction_plus_override", "0.024"),
            ("size_systematic_fraction_minus", "0.06"),
            ("size_systematic_fraction_plus", "0.06"),
            ("systematic_status", "none_adopted"),
            ("mass_reference_doi", "10.0000/wrong"),
            ("mass_uncertainty_reference_doi", "10.0000/wrong"),
            ("distance_reference_doi", "10.0000/wrong"),
            ("uncertainty_note", "different_note"),
            ("mass_display_interval_status", "decomposed_components"),
            ("distance_stat_component_status", "not_available"),
            ("mass_systematic_status", "not_available"),
        )
        for field, value in source_mutations:
            with self.subTest(source_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in self.inputs]
                next(
                    item for item in rows if item["source_key"] == "ESO558-G009"
                )[field] = value
                path = self._write_source_inputs(tmpdir, rows)
                with self.assertRaisesRegex(ValueError, "ESO558-G009"):
                    bhs.load_source_inputs(path)

        with bhs.DEFAULT_NONADOPTED_FITS.open(newline="", encoding="utf-8") as handle:
            original_fits = list(csv.DictReader(handle))
        fit_mutations = (
            ("fit_key", "villaraos_2022_eso558g009_other_fit"),
            ("data_relation", "reanalysis_of_unpublished_data"),
            ("model_basis", "unknown_model"),
            ("ratio_msun_per_mpc", "159401"),
            ("ratio_unc_minus_msun_per_mpc", "2599"),
            ("ratio_unc_plus_msun_per_mpc", "2501"),
            ("interval_status", "unconditional"),
            ("adoption_status", "adopted"),
            ("nonadoption_basis", "none"),
            ("reference_key", "GonzalezJuarez2024"),
            ("reference_doi", "10.0000/wrong"),
            ("provenance_note", "published_table_1_ratio"),
        )
        for field, value in fit_mutations:
            with self.subTest(fit_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in original_fits]
                next(
                    item for item in rows if item["source_key"] == "ESO558-G009"
                )[field] = value
                path = self._write_nonadopted_fits(tmpdir, rows)
                with self.assertRaisesRegex(ValueError, "ESO558-G009"):
                    bhs.load_nonadopted_ratio_fits(
                        path, source_treatments=self.source_treatments
                    )

        output_mutations = (
            ("distance_stat_unc_minus_mpc", 3.4026),
            ("distance_unc_plus_mpc", 5.495),
            ("mass_fit_unc_minus_msun", 386_239.0),
            ("mass_shared_distance_unc_plus_msun", 868_183.0),
            ("mass_systematic_unc_minus_msun", 1_008_138.0),
            ("log_mass_systematic_unc_plus", 0.0253),
            ("mass_display_unc_minus_msun", 950_222.0),
            ("mass_total_unc_plus_msun", 1_385_375.0),
            ("size_stat_unc_minus_uas", 0.000373),
            ("size_systematic_unc_plus_uas", 0.000973),
            ("size_total_unc_plus_uas", 0.001042),
            ("mass_systematic_status", "not_available"),
            ("mass_display_interval_status", "decomposed_components"),
            ("size_stat_component_status", "fixed_distance_formal_fit_interval"),
            ("size_systematic_status", "none_adopted"),
            ("size_total_status", "formal_component_only"),
            ("nonadopted_fit_ratio_msun_per_mpc", 159_401.0),
            ("nonadopted_fit_reference_doi", "10.0000/wrong"),
        )
        for field, value in output_mutations:
            with self.subTest(output_field=field):
                rows = [dict(item) for item in self.rows]
                next(
                    item for item in rows if item["source_key"] == "ESO558-G009"
                )[field] = value
                with self.assertRaisesRegex(ValueError, "ESO558-G009"):
                    bhs.render_csv(rows)

    def test_cgcg_source_quantiles_and_conditional_fit_fail_closed(self) -> None:
        source_mutations = (
            ("source_label", "CGCG 074-64"),
            ("reference_mass_msun", "24200001"),
            ("reference_mass_unc_minus_msun", "2000001"),
            ("reference_mass_unc_plus_msun", "2200001"),
            ("reference_mass_systematic_unc_minus_msun", "500001"),
            ("reference_mass_systematic_unc_plus_msun", "500001"),
            ("mass_systematic_status", "not_available"),
            ("reference_distance_mpc", "87.7"),
            ("adopted_distance_mpc", "87.7"),
            ("distance_measurement_unc_minus_mpc", "7.3"),
            ("distance_measurement_unc_plus_mpc", "8.0"),
            ("distance_flow_unc_mpc", "0.1"),
            ("distance_shared_systematic_unc_mpc", "1.6"),
            ("size_stat_fraction_minus_override", "0.01"),
            ("size_stat_fraction_plus_override", "0.01"),
            ("size_systematic_fraction_minus", "0.01"),
            ("size_systematic_fraction_plus", "0.01"),
            ("systematic_status", "none_adopted"),
            ("mass_reference_doi", "10.0000/wrong"),
            ("mass_uncertainty_reference_doi", "10.0000/wrong"),
            ("distance_reference_doi", "10.0000/wrong"),
            ("uncertainty_note", "different_note"),
            ("mass_display_interval_status", "decomposed_components"),
            ("distance_stat_component_status", "not_available"),
        )
        for field, value in source_mutations:
            with self.subTest(source_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in self.inputs]
                next(item for item in rows if item["source_key"] == "CGCG074-064")[
                    field
                ] = value
                path = self._write_source_inputs(tmpdir, rows)
                with self.assertRaises(ValueError):
                    bhs.load_source_inputs(path)

        with bhs.DEFAULT_QUANTILES.open(newline="", encoding="utf-8") as handle:
            quantiles = list(csv.DictReader(handle))
        quantile_mutations = (
            ("q16_theta_g_nas", "2.7102"),
            ("q50_theta_g_nas", "2.7278"),
            ("q84_theta_g_nas", "2.7475"),
            ("systematic_fraction", "0.011"),
            ("provenance", "different_provenance"),
        )
        for field, value in quantile_mutations:
            with self.subTest(quantile_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in quantiles]
                next(item for item in rows if item["source_key"] == "CGCG074-064")[
                    field
                ] = value
                path = self._write_quantiles(tmpdir, rows)
                with self.assertRaisesRegex(ValueError, "CGCG074-064"):
                    bhs.load_paired_quantiles(path)

        with bhs.DEFAULT_NONADOPTED_FITS.open(
            newline="", encoding="utf-8"
        ) as handle:
            conditional_fits = list(csv.DictReader(handle))
        fit_mutations = (
            ("fit_key", "different_fit"),
            ("data_relation", "reanalysis_of_published_kuo2011_vlbi_spot_data"),
            ("model_basis", "unknown_model"),
            ("ratio_msun_per_mpc", "274901"),
            ("ratio_unc_minus_msun_per_mpc", "1501"),
            ("ratio_unc_plus_msun_per_mpc", "1701"),
            ("interval_status", "unconditional"),
            ("adoption_status", "adopted"),
            (
                "nonadoption_basis",
                "same_data_reanalysis_conditional_on_fixed_vertical_center_and_"
                "selected_scatter_angle",
            ),
            ("reference_key", "GonzalezJuarez2024"),
            ("reference_doi", "10.0000/wrong"),
            ("provenance_note", "published_table_2_ratio_and_formal_interval"),
        )
        for field, value in fit_mutations:
            with self.subTest(fit_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in conditional_fits]
                next(item for item in rows if item["source_key"] == "CGCG074-064")[
                    field
                ] = value
                path = self._write_nonadopted_fits(tmpdir, rows)
                with self.assertRaisesRegex(ValueError, "CGCG074-064"):
                    bhs.load_nonadopted_ratio_fits(
                        path, source_treatments=self.source_treatments
                    )

    def test_ngc5765b_source_and_paired_quantile_provenance_fail_closed(self) -> None:
        source_mutations = (
            ("source_label", "NGC 5765B"),
            ("reference_mass_msun", "41500001"),
            ("reference_mass_unc_minus_msun", "1900001"),
            ("reference_mass_unc_plus_msun", "2000001"),
            ("reference_distance_mpc", "112.3"),
            ("adopted_distance_mpc", "112.3"),
            ("distance_measurement_unc_minus_mpc", "5.2"),
            ("distance_measurement_unc_plus_mpc", "5.5"),
            ("distance_flow_unc_mpc", "0.1"),
            ("distance_shared_systematic_unc_mpc", "0.1"),
            ("size_stat_fraction_minus_override", "0.01"),
            ("size_stat_fraction_plus_override", "0.01"),
            ("size_systematic_fraction_minus", "0.01"),
            ("size_systematic_fraction_plus", "0.01"),
            ("systematic_status", "none_adopted"),
            ("mass_reference_doi", "10.0000/wrong"),
            ("mass_uncertainty_reference_doi", "10.0000/wrong"),
            ("distance_reference_doi", "10.0000/wrong"),
            ("uncertainty_note", "different_note"),
        )
        for field, value in source_mutations:
            with self.subTest(source_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in self.inputs]
                next(item for item in rows if item["source_key"] == "NGC5765b")[
                    field
                ] = value
                path = self._write_source_inputs(tmpdir, rows)
                with self.assertRaisesRegex(ValueError, "NGC5765b"):
                    bhs.load_source_inputs(path)

        with bhs.DEFAULT_QUANTILES.open(newline="", encoding="utf-8") as handle:
            originals = list(csv.DictReader(handle))
        quantile_mutations = (
            ("q16_theta_g_nas", "3.6374"),
            ("q50_theta_g_nas", "3.6494"),
            ("q84_theta_g_nas", "3.6624"),
            ("systematic_fraction", "0.011"),
            ("provenance", "different_provenance"),
        )
        for field, value in quantile_mutations:
            with self.subTest(quantile_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in originals]
                next(item for item in rows if item["source_key"] == "NGC5765b")[
                    field
                ] = value
                path = self._write_quantiles(tmpdir, rows)
                with self.assertRaisesRegex(ValueError, "NGC5765b"):
                    bhs.load_paired_quantiles(path)

        with bhs.DEFAULT_NONADOPTED_FITS.open(
            newline="", encoding="utf-8"
        ) as handle:
            originals = list(csv.DictReader(handle))
        comparison_mutations = (
            ("fit_key", "different_fit"),
            ("data_relation", "reanalysis_of_published_kuo2011_vlbi_spot_data"),
            ("model_basis", "unknown_model"),
            ("ratio_msun_per_mpc", "372701"),
            ("ratio_unc_minus_msun_per_mpc", "1301"),
            ("ratio_unc_plus_msun_per_mpc", "1301"),
            ("interval_status", "unconditional"),
            ("adoption_status", "adopted"),
            (
                "nonadoption_basis",
                "same_data_reanalysis_conditional_on_fixed_vertical_center_and_"
                "selected_scatter_angle",
            ),
            ("reference_key", "GonzalezJuarez2024"),
            ("reference_doi", "10.0000/wrong"),
            ("provenance_note", "different_provenance"),
        )
        for field, value in comparison_mutations:
            with self.subTest(comparison_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in originals]
                next(item for item in rows if item["source_key"] == "NGC5765b")[
                    field
                ] = value
                path = self._write_nonadopted_fits(tmpdir, rows)
                with self.assertRaisesRegex(ValueError, "NGC5765b"):
                    bhs.load_nonadopted_ratio_fits(
                        path, source_treatments=self.source_treatments
                    )

    def test_ngc4258_model_check_is_separate_and_reproducible(self) -> None:
        self.assertEqual(len(self.model_checks), 10)
        self.assertEqual(
            set(self.model_checks),
            {
                "NGC4258",
                "NGC4388",
                "Circinus",
                "NGC3079",
                "NGC4945",
                "NGC2960",
                "NGC6323",
                "J0437",
                "Mrk1029",
                "NGC1320",
            },
        )
        check = self.model_checks["NGC4258"]
        row = self.by_key["NGC4258"]

        expected_lower = bhs.shadow_size_uas(4.234e7, 8.091)
        expected_upper = bhs.shadow_size_uas(3.843e7, 7.182)
        self.assertEqual(check["size_lower_uas"], expected_lower)
        self.assertEqual(check["size_upper_uas"], expected_upper)
        self.assertEqual(check["accepted_configuration_count"], 15)
        self.assertEqual(check["excluded_configuration_ids"], "2;17")
        self.assertEqual(check["range_basis"], "accepted_configuration_medians")
        self.assertEqual(check["lower_configuration_id"], 1)
        self.assertEqual(check["upper_configuration_id"], 8)
        self.assertEqual(check["probability_status"], "non_probabilistic")
        self.assertEqual(check["adoption_status"], "reported_separately_not_adopted")
        self.assertEqual(check["reference_doi"], "10.1093/mnras/stag1327")

        self.assertEqual(row["size_model_check_lower_uas"], expected_lower)
        self.assertEqual(row["size_model_check_upper_uas"], expected_upper)
        self.assertEqual(
            row["size_model_check_data_relation"],
            "subset_reanalysis_of_argon2007_vlbi_epochs_used_by_primary",
        )
        self.assertEqual(row["size_model_check_probability_status"], "non_probabilistic")
        self.assertEqual(
            row["size_model_check_adoption_status"],
            "reported_separately_not_adopted",
        )

        # Freeze the adopted result: the separate check must never enter these fields.
        self.assertEqual(row["size_posterior_q16_uas"], 0.53804287722255006)
        self.assertEqual(row["size_posterior_q50_uas"], 0.53860302245371783)
        self.assertEqual(row["size_posterior_q84_uas"], 0.53921720767008174)
        self.assertEqual(row["size_stat_unc_minus_uas"], 0.00056014523116776349)
        self.assertEqual(row["size_stat_unc_plus_uas"], 0.0006141852163639161)
        self.assertEqual(row["size_systematic_unc_minus_uas"], 0.0053860302245371788)
        self.assertEqual(row["size_systematic_unc_plus_uas"], 0.0053860302245371788)
        self.assertEqual(row["size_total_unc_minus_uas"], 0.0054150793401046308)
        self.assertEqual(row["size_total_unc_plus_uas"], 0.005420935810321683)
        self.assertEqual(row["size_unc_adopted_uas"], 0.0054180075752131573)
        self.assertEqual(row["size_total_status"], "quadrature_statistical_plus_systematic")
        self.assertIn("not combined", row["uncertainty_status"])

        for key, other in self.by_key.items():
            if key in {
                "NGC4258",
                "NGC4388",
                "Circinus",
                "NGC3079",
                "NGC4945",
                "NGC2960",
                "NGC6323",
                "J0437",
                "Mrk1029",
                "NGC1320",
            }:
                continue
            self.assertIsNone(other["size_model_check_key"], key)
            self.assertIsNone(other["size_model_check_lower_uas"], key)
            self.assertIsNone(other["size_model_check_upper_uas"], key)

    def test_ngc2960_formal_astrometric_model_and_villaraos_results_are_exact_and_separate(
        self,
    ) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "NGC2960")
        baseline = bhs.build_shared_row(source)
        row = self.by_key["NGC2960"]
        check = self.model_checks["NGC2960"]
        alternate = self.nonadopted_fits["NGC2960"]

        expected_source = {
            "source_order": "15",
            "source_label": "NGC 2960",
            "treatment": "shared_distance",
            "reference_mass_msun": "12370000",
            "reference_mass_unc_minus_msun": "60000",
            "reference_mass_unc_plus_msun": "70000",
            "reference_distance_mpc": "72.2",
            "adopted_distance_mpc": "71.32070365358592",
            "distance_measurement_unc_minus_mpc": "0.20297699594046006",
            "distance_measurement_unc_plus_mpc": "0.20297699594046006",
            "distance_flow_unc_mpc": "3.3829499323410013",
            "distance_shared_systematic_unc_mpc": "2.8952924351929332",
            "size_stat_fraction_minus_override": "",
            "size_stat_fraction_plus_override": "",
            "size_systematic_fraction_minus": "",
            "size_systematic_fraction_plus": "",
            "systematic_status": "unquantified_astrometric",
            "mass_reference_doi": "10.3847/1538-4357/aabff1",
            "mass_uncertainty_reference_doi": "10.3847/1538-4357/aabff1",
            "distance_reference_doi": "10.1088/0004-637X/727/1/20",
            "uncertainty_note": (
                "later_three_dimensional_disk_fit_formal_interval_only_absolute_"
                "coordinate_error_floor_completeness_unreported_legacy_three_to_"
                "five_percent_not_transplanted_disk_model_checks_not_adopted"
            ),
            "mass_display_interval_status": "decomposed_components",
            "distance_stat_component_status": "statistical",
            "reference_mass_systematic_unc_minus_msun": "",
            "reference_mass_systematic_unc_plus_msun": "",
            "mass_systematic_status": "not_available",
        }
        for field, expected in expected_source.items():
            self.assertEqual(source[field], expected, field)

        expected_primary = {
            "distance_mpc": 71.32070365358592,
            "distance_stat_unc_minus_mpc": 3.389033771682867,
            "distance_stat_unc_plus_mpc": 3.389033771682867,
            "distance_systematic_unc_minus_mpc": 2.8952924351929332,
            "distance_systematic_unc_plus_mpc": 2.8952924351929332,
            "distance_unc_minus_mpc": 4.457383558870879,
            "distance_unc_plus_mpc": 4.457383558870879,
            "mass_msun": 12_219_350.47361299,
            "mass_fit_unc_minus_msun": 59_269.282814614344,
            "mass_fit_unc_plus_msun": 69_147.49661705006,
            "mass_shared_distance_unc_minus_msun": 763_681.9199893735,
            "mass_shared_distance_unc_plus_msun": 763_681.9199893735,
            "mass_total_unc_minus_msun": 765_978.4088367078,
            "mass_total_unc_plus_msun": 766_806.006240862,
            "mass_display_unc_minus_msun": 765_978.4088367078,
            "mass_display_unc_plus_msun": 766_806.006240862,
            "size_uas": 0.017574751618808785,
            "size_unc_adopted_uas": 9.234913946827575e-05,
            "size_stat_unc_minus_uas": 8.524535950917762e-05,
            "size_stat_unc_plus_uas": 9.945291942737387e-05,
        }
        for field, expected in expected_primary.items():
            self.assertEqual(row[field], expected, field)

        self.assertEqual(
            row["size_stat_component_status"], "fixed_distance_formal_fit_interval"
        )
        self.assertEqual(row["size_systematic_status"], "unquantified_astrometric")
        self.assertEqual(
            row["size_total_status"], "incomplete_unquantified_systematic"
        )
        self.assertIsNone(row["size_systematic_unc_minus_uas"])
        self.assertIsNone(row["size_systematic_unc_plus_uas"])
        self.assertIsNone(row["size_total_unc_minus_uas"])
        self.assertIsNone(row["size_total_unc_plus_uas"])
        self.assertEqual(
            row["legacy_uncertainty_basis"],
            "symmetrized_formal_statistical_term_total_incomplete",
        )
        self.assertIn("absolute-coordinate", row["uncertainty_status"])
        self.assertIn("disk-model sensitivity not combined", row["uncertainty_status"])

        self.assertEqual(
            check["check_key"],
            "kuo_2018_ngc2960_physically_allowed_disk_mass_models",
        )
        self.assertEqual(check["check_type"], "disk_mass_model_choice")
        self.assertEqual(
            check["data_relation"],
            "same_fit_data_alternative_accretion_disk_models",
        )
        self.assertEqual(
            check["range_basis"], "published_physically_allowed_model_envelope"
        )
        self.assertEqual(check["endpoint_derivation"], "published_mass_distance_endpoints")
        self.assertEqual(check["sensitivity_direction"], "downward_only")
        self.assertEqual(check["accepted_configuration_count"], 3)
        self.assertEqual(check["excluded_configuration_ids"], "2;4")
        self.assertEqual(check["lower_configuration_id"], 3)
        self.assertEqual(check["lower_mass_msun"], 12_320_000.0)
        self.assertEqual(check["lower_distance_mpc"], 72.2)
        self.assertEqual(check["upper_configuration_id"], 1)
        self.assertEqual(check["upper_mass_msun"], 12_370_000.0)
        self.assertEqual(check["upper_distance_mpc"], 72.2)
        self.assertEqual(check["size_lower_uas"], 0.017503713819217805)
        self.assertEqual(check["size_upper_uas"], 0.017574751618808785)
        self.assertEqual(check["probability_status"], "non_probabilistic")
        self.assertEqual(check["adoption_status"], "reported_separately_not_adopted")
        self.assertEqual(check["reference_key"], "Kuo2018Disk")
        self.assertEqual(check["reference_doi"], "10.3847/1538-4357/aabff1")
        self.assertEqual(row["size_model_check_lower_uas"], check["size_lower_uas"])
        self.assertEqual(row["size_model_check_upper_uas"], check["size_upper_uas"])
        self.assertAlmostEqual(
            100.0
            * (check["size_upper_uas"] - check["size_lower_uas"])
            / check["size_upper_uas"],
            0.4042037186742026,
            places=15,
        )

        self.assertEqual(
            alternate["fit_key"], "villaraos_2022_ngc2960_gr_ratio_fit"
        )
        self.assertEqual(alternate["ratio_msun_per_mpc"], 173_800.0)
        self.assertEqual(alternate["ratio_unc_minus_msun_per_mpc"], 1_600.0)
        self.assertEqual(alternate["ratio_unc_plus_msun_per_mpc"], 1_500.0)
        self.assertEqual(alternate["size_uas"], 0.017828157657509735)
        self.assertEqual(alternate["size_unc_minus_uas"], 0.0001641257321750033)
        self.assertEqual(alternate["size_unc_plus_uas"], 0.0001538678739140656)
        self.assertEqual(
            alternate["nonadoption_basis"],
            "older_published_data_reanalysis_conditional_on_fixed_kuo2011_"
            "inclination_vertical_center_and_selected_scatter_angle",
        )
        self.assertEqual(alternate["reference_key"], "Villaraos2022")
        self.assertEqual(alternate["reference_doi"], "10.1093/mnras/stac2973")
        self.assertIn("y0_variation_up_to_4_percent", alternate["provenance_note"])
        self.assertIn("exact_spot_overlap_with_kuo2018_unverifiable", alternate["provenance_note"])
        self.assertEqual(row["nonadopted_fit_size_uas"], alternate["size_uas"])
        self.assertNotEqual(row["size_uas"], row["nonadopted_fit_size_uas"])
        self.assertAlmostEqual(
            100.0 * (row["nonadopted_fit_size_uas"] / row["size_uas"] - 1.0),
            1.4418755052546528,
            places=14,
        )

        # Neither sidecar may alter the adopted formal-only primary result.
        for field in bhs.OUTPUT_FIELDS:
            if (
                field.startswith("size_model_check_")
                or field.startswith("nonadopted_fit_")
                or field == "uncertainty_status"
            ):
                continue
            self.assertEqual(row[field], baseline[field], field)

    def test_ngc2960_reviewed_inputs_and_sidecars_fail_closed(self) -> None:
        source_mutations = (
            ("reference_mass_msun", "12370001"),
            ("reference_mass_unc_minus_msun", "60001"),
            ("reference_mass_unc_plus_msun", "70001"),
            ("reference_distance_mpc", "72.3"),
            ("adopted_distance_mpc", "71.4"),
            ("distance_measurement_unc_minus_mpc", "0.203"),
            ("distance_measurement_unc_plus_mpc", "0.203"),
            ("distance_flow_unc_mpc", "3.4"),
            ("distance_shared_systematic_unc_mpc", "2.9"),
            ("size_stat_fraction_minus_override", "0.01"),
            ("size_systematic_fraction_minus", "0.03"),
            ("systematic_status", "unquantified_model"),
            ("mass_reference_doi", "10.0000/wrong"),
            ("mass_uncertainty_reference_doi", "10.0000/wrong"),
            ("distance_reference_doi", "10.0000/wrong"),
            ("uncertainty_note", "different_note"),
            ("mass_display_interval_status", "source_reported_total"),
            ("distance_stat_component_status", "not_available"),
            ("reference_mass_systematic_unc_minus_msun", "300000"),
            ("mass_systematic_status", "source_reported_absolute_position_systematic"),
        )
        for field, value in source_mutations:
            with self.subTest(source_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in self.inputs]
                next(item for item in rows if item["source_key"] == "NGC2960")[
                    field
                ] = value
                path = self._write_source_inputs(tmpdir, rows)
                with self.assertRaises(ValueError):
                    bhs.load_source_inputs(path)

        with bhs.DEFAULT_MODEL_CHECKS.open(newline="", encoding="utf-8") as handle:
            model_rows = list(csv.DictReader(handle))
        reviewed_model = next(
            row for row in model_rows if row["source_key"] == "NGC2960"
        )
        model_mutations = (
            ("check_key", "different_check"),
            ("check_type", "data_selection_component_tracking"),
            ("data_relation", "subset_reanalysis_of_argon2007_vlbi_epochs_used_by_primary"),
            ("range_basis", "accepted_configuration_medians"),
            ("endpoint_derivation", "adopted_primary_ratio_factors"),
            ("sensitivity_direction", "two_sided_range"),
            ("accepted_configuration_count", "2"),
            ("excluded_configuration_ids", "2;5"),
            ("lower_configuration_id", "5"),
            ("lower_mass_msun", "12320001"),
            ("lower_distance_mpc", "72.3"),
            ("upper_configuration_id", "5"),
            ("upper_mass_msun", "12370001"),
            ("upper_distance_mpc", "72.3"),
            ("lower_primary_ratio_factor", "1"),
            ("probability_status", "probabilistic"),
            ("adoption_status", "adopted"),
            ("reference_key", "Kuo2020"),
            ("reference_doi", "10.0000/wrong"),
            ("provenance_note", "different_provenance"),
        )
        for field, value in model_mutations:
            with self.subTest(model_field=field), tempfile.TemporaryDirectory() as tmpdir:
                mutated = dict(reviewed_model)
                mutated[field] = value
                path = self._write_model_checks(tmpdir, [mutated])
                with self.assertRaises(ValueError):
                    bhs.load_model_checks(path, known_sources=self.source_keys)

        with bhs.DEFAULT_NONADOPTED_FITS.open(
            newline="", encoding="utf-8"
        ) as handle:
            comparison_rows = list(csv.DictReader(handle))
        comparison_mutations = (
            ("fit_key", "different_fit"),
            ("data_relation", "reanalysis_of_published_gao2017_vlbi_spot_data"),
            ("model_basis", "unknown_model"),
            ("ratio_msun_per_mpc", "173801"),
            ("ratio_unc_minus_msun_per_mpc", "1601"),
            ("ratio_unc_plus_msun_per_mpc", "1501"),
            ("interval_status", "unconditional"),
            ("adoption_status", "adopted"),
            (
                "nonadoption_basis",
                "same_data_reanalysis_conditional_on_fixed_vertical_center_and_"
                "selected_scatter_angle",
            ),
            ("reference_key", "GonzalezJuarez2024"),
            ("reference_doi", "10.0000/wrong"),
            ("provenance_note", "different_provenance"),
        )
        for field, value in comparison_mutations:
            with self.subTest(comparison_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in comparison_rows]
                next(item for item in rows if item["source_key"] == "NGC2960")[
                    field
                ] = value
                path = self._write_nonadopted_fits(tmpdir, rows)
                with self.assertRaisesRegex(ValueError, "NGC2960"):
                    bhs.load_nonadopted_ratio_fits(
                        path, source_treatments=self.source_treatments
                    )

    def test_ngc2960_sidecar_endpoint_changes_cannot_change_the_primary(self) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "NGC2960")
        baseline = bhs.build_shared_row(source)

        model = dict(self.model_checks["NGC2960"])
        model["size_lower_uas"] = 0.001
        model["size_upper_uas"] = 9.9
        candidate = bhs.attach_model_check(dict(baseline), model)

        alternate = dict(self.nonadopted_fits["NGC2960"])
        alternate["ratio_msun_per_mpc"] = 999_999.0
        alternate["size_uas"] = bhs.shadow_size_uas(999_999.0, 1.0)
        candidate = bhs.attach_nonadopted_ratio_fit(candidate, alternate)

        self.assertEqual(candidate["size_model_check_lower_uas"], 0.001)
        self.assertEqual(candidate["size_model_check_upper_uas"], 9.9)
        self.assertEqual(
            candidate["nonadopted_fit_size_uas"],
            bhs.shadow_size_uas(999_999.0, 1.0),
        )
        for field in bhs.OUTPUT_FIELDS:
            if (
                field.startswith("size_model_check_")
                or field.startswith("nonadopted_fit_")
                or field == "uncertainty_status"
            ):
                continue
            self.assertEqual(candidate[field], baseline[field], field)

    def test_ngc6323_paired_primary_and_two_checks_are_exact_and_separate(
        self,
    ) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "NGC6323")
        paired = self.quantiles["NGC6323"]
        baseline = bhs.build_paired_row(source, paired)
        row = self.by_key["NGC6323"]
        check = self.model_checks["NGC6323"]
        alternate = self.nonadopted_fits["NGC6323"]

        expected_source = {
            "source_order": "18",
            "source_label": "NGC 6323",
            "treatment": "paired_posterior",
            "reference_mass_msun": "10200000",
            "reference_mass_unc_minus_msun": "2200000",
            "reference_mass_unc_plus_msun": "3200000",
            "reference_distance_mpc": "109.4",
            "adopted_distance_mpc": "109.4",
            "distance_measurement_unc_minus_mpc": "23.4",
            "distance_measurement_unc_plus_mpc": "34.2",
            "distance_flow_unc_mpc": "0",
            "distance_shared_systematic_unc_mpc": "0",
            "size_stat_fraction_minus_override": "",
            "size_stat_fraction_plus_override": "",
            "size_systematic_fraction_minus": "",
            "size_systematic_fraction_plus": "",
            "systematic_status": "paired_quantile_file",
            "mass_reference_doi": "10.3847/2041-8213/ab75f0",
            "mass_uncertainty_reference_doi": "10.3847/2041-8213/ab75f0",
            "distance_reference_doi": "10.3847/2041-8213/ab75f0",
            "uncertainty_note": (
                "paired_theta_g_quantiles_plus_one_percent_systematic"
            ),
            "mass_display_interval_status": "decomposed_components",
            "distance_stat_component_status": "statistical",
            "reference_mass_systematic_unc_minus_msun": "",
            "reference_mass_systematic_unc_plus_msun": "",
            "mass_systematic_status": "not_available",
        }
        for field, expected in expected_source.items():
            self.assertEqual(source[field], expected, field)

        expected_quantiles = {
            "q16_theta_g_nas": 0.9135,
            "q50_theta_g_nas": 0.9189,
            "q84_theta_g_nas": 0.9248,
            "systematic_fraction": 0.01,
            "provenance": "D_Pesce_private_communication_2026-08-17",
        }
        for field, expected in expected_quantiles.items():
            self.assertEqual(paired[field], expected, field)

        expected_primary = {
            "size_uas": 0.00954948892245025,
            "size_posterior_q16_uas": 0.009493370476285018,
            "size_posterior_q50_uas": 0.00954948892245025,
            "size_posterior_q84_uas": 0.009610803521038186,
            "size_stat_unc_minus_uas": 5.611844616523233e-05,
            "size_stat_unc_plus_uas": 6.131459858793598e-05,
            "size_systematic_unc_minus_uas": 9.549488922450251e-05,
            "size_systematic_unc_plus_uas": 9.549488922450251e-05,
            "size_total_unc_minus_uas": 0.00011076350422409036,
            "size_total_unc_plus_uas": 0.00011348459749234575,
            "size_unc_adopted_uas": 0.00011212405085821805,
            "mass_msun": 10_200_000.0,
            "mass_unc_minus_msun": 2_200_000.0,
            "mass_unc_plus_msun": 3_200_000.0,
            "distance_mpc": 109.4,
            "distance_unc_minus_mpc": 23.4,
            "distance_unc_plus_mpc": 34.2,
        }
        for field, expected in expected_primary.items():
            self.assertEqual(row[field], expected, field)
        self.assertEqual(row["treatment"], "paired_posterior")
        self.assertEqual(
            row["size_stat_component_status"], "paired_chain_quantile_interval"
        )
        self.assertEqual(
            row["size_systematic_status"],
            "numeric_one_percent_non_keplerian_allowance",
        )
        self.assertEqual(
            row["size_total_status"], "quadrature_statistical_plus_systematic"
        )

        self.assertEqual(
            check["check_key"],
            "kuo_2018_ngc6323_physically_allowed_disk_mass_models",
        )
        self.assertEqual(check["check_type"], "disk_mass_model_choice")
        self.assertEqual(
            check["data_relation"],
            "predecessor_fit_data_alternative_accretion_disk_models",
        )
        self.assertEqual(
            check["range_basis"], "published_physically_allowed_model_envelope"
        )
        self.assertEqual(check["endpoint_derivation"], "published_mass_distance_endpoints")
        self.assertEqual(check["sensitivity_direction"], "downward_only")
        self.assertEqual(check["accepted_configuration_count"], 3)
        self.assertEqual(check["excluded_configuration_ids"], "2;4")
        self.assertEqual(check["lower_configuration_id"], 3)
        self.assertEqual(check["lower_mass_msun"], 9_840_000.0)
        self.assertEqual(check["lower_distance_mpc"], 106.7)
        self.assertEqual(check["upper_configuration_id"], 1)
        self.assertEqual(check["upper_mass_msun"], 9_920_000.0)
        self.assertEqual(check["upper_distance_mpc"], 106.7)
        self.assertEqual(check["size_lower_uas"], 0.009459918021333369)
        self.assertEqual(check["size_upper_uas"], 0.009536827923945834)
        self.assertEqual(check["probability_status"], "non_probabilistic")
        self.assertEqual(check["adoption_status"], "reported_separately_not_adopted")
        self.assertEqual(check["reference_key"], "Kuo2018Disk")
        self.assertEqual(check["reference_doi"], "10.3847/1538-4357/aabff1")
        self.assertEqual(
            check["provenance_note"],
            "table_2_configuration_ids_1_fiducial_2_linear_no_density_excluded_"
            "3_linear_physical_lower_4_herrnstein_no_density_excluded_5_"
            "herrnstein_physical_intermediate_upper_endpoint_fiducial",
        )
        self.assertEqual(row["size_model_check_lower_uas"], check["size_lower_uas"])
        self.assertEqual(row["size_model_check_upper_uas"], check["size_upper_uas"])
        self.assertAlmostEqual(
            100.0
            * (check["size_upper_uas"] - check["size_lower_uas"])
            / check["size_upper_uas"],
            0.806451612903214,
            places=14,
        )

        self.assertEqual(
            alternate["fit_key"], "villaraos_2022_ngc6323_gr_ratio_fit"
        )
        self.assertEqual(
            alternate["data_relation"],
            "reanalysis_of_published_kuo2011_vlbi_spot_data",
        )
        self.assertEqual(
            alternate["model_basis"], "circular_equatorial_schwarzschild_disk"
        )
        self.assertEqual(alternate["ratio_msun_per_mpc"], 91_600.0)
        self.assertEqual(alternate["ratio_unc_minus_msun_per_mpc"], 500.0)
        self.assertEqual(alternate["ratio_unc_plus_msun_per_mpc"], 500.0)
        self.assertEqual(alternate["size_uas"], 0.00939619816701894)
        self.assertEqual(alternate["size_unc_minus_uas"], 5.128929130468853e-05)
        self.assertEqual(alternate["size_unc_plus_uas"], 5.128929130468853e-05)
        self.assertEqual(
            alternate["nonadoption_basis"],
            "older_partial_data_reanalysis_conditional_on_fixed_kuo2015_"
            "inclination_and_selected_scatter_angle",
        )
        self.assertEqual(
            alternate["provenance_note"],
            "published_table_2_ratio_and_formal_interval_four_parameters_"
            "y0_fitted_table_2_scatter_11deg_chi2red_1p276_conflicts_with_"
            "section_4p1p2_chi2red_1p251",
        )
        self.assertEqual(row["nonadopted_fit_size_uas"], alternate["size_uas"])
        self.assertAlmostEqual(
            100.0 * (row["nonadopted_fit_size_uas"] / row["size_uas"] - 1.0),
            -1.6052247054911262,
            places=14,
        )

        # Both sidecars are attached after the paired primary is constructed.
        for field in bhs.OUTPUT_FIELDS:
            if (
                field.startswith("size_model_check_")
                or field.startswith("nonadopted_fit_")
                or field == "uncertainty_status"
            ):
                continue
            self.assertEqual(row[field], baseline[field], field)

    def test_ngc6323_reviewed_inputs_sidecars_and_output_fail_closed(self) -> None:
        source_mutations = (
            ("reference_mass_msun", "10200001"),
            ("reference_mass_unc_plus_msun", "3200001"),
            ("adopted_distance_mpc", "109.5"),
            ("distance_measurement_unc_minus_mpc", "23.5"),
            ("systematic_status", "none_adopted"),
            ("uncertainty_note", "different_note"),
            ("reference_mass_systematic_unc_minus_msun", "1"),
        )
        for field, value in source_mutations:
            with self.subTest(source_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in self.inputs]
                next(item for item in rows if item["source_key"] == "NGC6323")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "NGC6323"):
                    bhs.load_source_inputs(self._write_source_inputs(tmpdir, rows))

        with bhs.DEFAULT_QUANTILES.open(newline="", encoding="utf-8") as handle:
            quantile_rows = list(csv.DictReader(handle))
        quantile_mutations = (
            ("q16_theta_g_nas", "0.9136"),
            ("q50_theta_g_nas", "0.9190"),
            ("q84_theta_g_nas", "0.9249"),
            ("systematic_fraction", "0.02"),
            ("provenance", "different_provenance"),
        )
        for field, value in quantile_mutations:
            with self.subTest(quantile_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in quantile_rows]
                next(item for item in rows if item["source_key"] == "NGC6323")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "NGC6323"):
                    bhs.load_paired_quantiles(self._write_quantiles(tmpdir, rows))

        with bhs.DEFAULT_MODEL_CHECKS.open(newline="", encoding="utf-8") as handle:
            model_rows = list(csv.DictReader(handle))
        model_mutations = (
            ("data_relation", "same_fit_data_alternative_accretion_disk_models"),
            ("accepted_configuration_count", "2"),
            ("excluded_configuration_ids", "2;5"),
            ("lower_configuration_id", "5"),
            ("lower_mass_msun", "9840001"),
            ("upper_mass_msun", "9920001"),
            ("reference_doi", "10.0000/wrong"),
            ("provenance_note", "different_provenance"),
        )
        for field, value in model_mutations:
            with self.subTest(model_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in model_rows]
                next(item for item in rows if item["source_key"] == "NGC6323")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "NGC6323"):
                    bhs.load_model_checks(
                        self._write_model_checks(tmpdir, rows),
                        known_sources=self.source_keys,
                    )

        with bhs.DEFAULT_NONADOPTED_FITS.open(
            newline="", encoding="utf-8"
        ) as handle:
            comparison_rows = list(csv.DictReader(handle))
        comparison_mutations = (
            ("ratio_msun_per_mpc", "91601"),
            ("ratio_unc_minus_msun_per_mpc", "501"),
            (
                "nonadoption_basis",
                "older_partial_data_reanalysis_conditional_on_fixed_vertical_"
                "center_and_selected_scatter_angle",
            ),
            ("reference_doi", "10.0000/wrong"),
            ("provenance_note", "different_provenance"),
        )
        for field, value in comparison_mutations:
            with self.subTest(comparison_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in comparison_rows]
                next(item for item in rows if item["source_key"] == "NGC6323")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "NGC6323"):
                    bhs.load_nonadopted_ratio_fits(
                        self._write_nonadopted_fits(tmpdir, rows),
                        source_treatments=self.source_treatments,
                    )

        output_mutations = (
            ("size_uas", 0.00954948892245025 + 1e-12),
            ("size_posterior_q16_uas", 0.009493370476285018 + 1e-12),
            ("size_total_unc_plus_uas", 0.00011348459749234575 + 1e-12),
            ("mass_msun", 10_200_001.0),
            ("size_model_check_lower_uas", 0.009459918021333369 + 1e-12),
            ("nonadopted_fit_ratio_msun_per_mpc", 91_601.0),
            ("nonadopted_fit_size_uas", 0.00939619816701894 + 1e-12),
            ("nonadopted_fit_key", "different_fit"),
        )
        for field, value in output_mutations:
            with self.subTest(output_field=field):
                rows = [dict(item) for item in self.rows]
                next(item for item in rows if item["source_key"] == "NGC6323")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "NGC6323"):
                    bhs.render_csv(rows)

    def test_j0437_fuller_primary_and_restricted_fit_are_exact_and_separate(
        self,
    ) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "J0437")
        baseline = bhs.build_shared_row(source)
        row = self.by_key["J0437"]
        check = self.model_checks["J0437"]
        alternate = self.nonadopted_fits["J0437"]

        expected_source = {
            "source_order": "19",
            "source_label": "J0437+2456",
            "treatment": "shared_distance",
            "reference_mass_msun": "2860000",
            "reference_mass_unc_minus_msun": "200000",
            "reference_mass_unc_plus_msun": "200000",
            "reference_distance_mpc": "70",
            "adopted_distance_mpc": "64.50473612990527",
            "distance_measurement_unc_minus_mpc": "0.5412719891745602",
            "distance_measurement_unc_plus_mpc": "0.5412719891745602",
            "distance_flow_unc_mpc": "3.3829499323410013",
            "distance_shared_systematic_unc_mpc": "2.6185955127160461",
            "size_stat_fraction_minus_override": "",
            "size_stat_fraction_plus_override": "",
            "size_systematic_fraction_minus": "",
            "size_systematic_fraction_plus": "",
            "systematic_status": "unquantified_model",
            "mass_reference_doi": "10.3847/1538-4357/abde3d",
            "mass_uncertainty_reference_doi": "10.3847/1538-4357/abde3d",
            "distance_reference_doi": "10.3847/1538-4357/834/1/52",
            "uncertainty_note": (
                "updated_nine_parameter_high_velocity_only_thin_disk_fit_formal_"
                "interval_shared_external_distance_inclination_unconstrained_and_"
                "model_completeness_unquantified"
            ),
            "mass_display_interval_status": "decomposed_components",
            "distance_stat_component_status": "statistical",
            "reference_mass_systematic_unc_minus_msun": "",
            "reference_mass_systematic_unc_plus_msun": "",
            "mass_systematic_status": "not_available",
        }
        for field, expected in expected_source.items():
            self.assertEqual(source[field], expected, field)

        measurement = 40.0 / 73.9
        flow = 250.0 / 73.9
        distance = 65.3 * 73.0 / 73.9
        distance_stat = math.hypot(measurement, flow)
        distance_systematic = distance * 3.0 / 73.9
        distance_total = math.hypot(distance_stat, distance_systematic)
        mass = 2_860_000.0 * (distance / 70.0)
        mass_fit = mass * (200_000.0 / 2_860_000.0)
        mass_shared_distance = mass * distance_total / distance
        mass_total = math.hypot(mass_fit, mass_shared_distance)

        self.assertEqual(float(source["distance_measurement_unc_minus_mpc"]), measurement)
        self.assertEqual(float(source["distance_flow_unc_mpc"]), flow)
        self.assertEqual(float(source["adopted_distance_mpc"]), distance)
        self.assertEqual(
            float(source["distance_shared_systematic_unc_mpc"]),
            distance_systematic,
        )

        expected_primary = {
            "size_uas": 0.004191067803754549,
            "size_stat_unc_minus_uas": 0.0002930816645982202,
            "size_stat_unc_plus_uas": 0.0002930816645982202,
            "size_unc_adopted_uas": 0.0002930816645982202,
            "mass_msun": 2_635_479.2190218437,
            "mass_fit_unc_minus_msun": 184_299.24608544362,
            "mass_fit_unc_plus_msun": 184_299.24608544362,
            "mass_shared_distance_unc_minus_msun": 176_180.85157776746,
            "mass_shared_distance_unc_plus_msun": 176_180.85157776746,
            "mass_total_unc_minus_msun": 254_962.55523180307,
            "mass_total_unc_plus_msun": 254_962.55523180307,
            "mass_display_unc_minus_msun": 254_962.55523180307,
            "mass_display_unc_plus_msun": 254_962.55523180307,
            "distance_mpc": 64.50473612990527,
            "distance_stat_unc_minus_mpc": 3.4259780517380682,
            "distance_stat_unc_plus_mpc": 3.4259780517380682,
            "distance_systematic_unc_minus_mpc": 2.618595512716046,
            "distance_systematic_unc_plus_mpc": 2.618595512716046,
            "distance_unc_minus_mpc": 4.312118744910393,
            "distance_unc_plus_mpc": 4.312118744910393,
        }
        for field, expected in expected_primary.items():
            self.assertEqual(row[field], expected, field)
        self.assertEqual(row["mass_msun"], mass)
        self.assertEqual(row["mass_fit_unc_minus_msun"], mass_fit)
        self.assertEqual(
            row["mass_shared_distance_unc_minus_msun"], mass_shared_distance
        )
        self.assertEqual(row["mass_total_unc_minus_msun"], mass_total)
        self.assertEqual(row["distance_stat_unc_minus_mpc"], distance_stat)
        self.assertEqual(row["distance_unc_minus_mpc"], distance_total)
        self.assertEqual(row["treatment"], "shared_distance")
        self.assertEqual(
            row["size_stat_component_status"],
            "fixed_distance_formal_fit_interval",
        )
        self.assertEqual(row["size_systematic_status"], "unquantified_model")
        self.assertEqual(
            row["size_total_status"], "incomplete_unquantified_systematic"
        )
        self.assertEqual(
            row["legacy_uncertainty_basis"],
            "symmetrized_formal_statistical_term_total_incomplete",
        )
        for field in (
            "size_systematic_unc_minus_uas",
            "size_systematic_unc_plus_uas",
            "size_total_unc_minus_uas",
            "size_total_unc_plus_uas",
            "mass_systematic_unc_minus_msun",
            "mass_systematic_unc_plus_msun",
        ):
            self.assertIsNone(row[field], field)

        expected_lower = bhs.shadow_size_uas(1_620_000.0, 68.27)
        expected_upper = bhs.shadow_size_uas(2_080_000.0, 68.27)
        expected_model = {
            "check_key": "deepshikha_2024_j0437_mestel_disk_mass_model",
            "check_type": "disk_mass_model_choice",
            "data_relation": (
                "reanalysis_of_published_gao2017_vlbi_data_data_selection_"
                "unreported"
            ),
            "range_basis": "reported_symmetric_mass_error_span",
            "endpoint_derivation": (
                "reported_mbh_plusminus_error_at_fixed_reported_distance"
            ),
            "sensitivity_direction": "downward_only",
            "accepted_configuration_count": 1,
            "excluded_configuration_ids": None,
            "lower_configuration_id": None,
            "lower_mass_msun": 1_620_000.0,
            "lower_distance_mpc": 68.27,
            "upper_configuration_id": None,
            "upper_mass_msun": 2_080_000.0,
            "upper_distance_mpc": 68.27,
            "lower_primary_ratio_factor": None,
            "upper_primary_ratio_factor": None,
            "assumed_inner_radius_pc": None,
            "assumed_turbulent_velocity_kms": None,
            "assumed_surface_density_power_law_index": -1.0,
            "assumed_rotation_profile": None,
            "mass_bias_eta": None,
            "assumption_profile": (
                "edge_on_coplanar_black_hole_plus_mestel_disk"
            ),
            "size_lower_uas": expected_lower,
            "size_upper_uas": expected_upper,
            "probability_status": "non_probabilistic",
            "adoption_status": "reported_separately_not_adopted",
            "reference_key": "Deepshikha2024",
            "reference_doi": "10.1017/S1743921323002491",
            "provenance_note": (
                "official_supplement_table_1_abs_alpha_0p303_plusminus0p037_"
                "mbh_1p85_plusminus0p23e6_msun_mdisk_0p82_plusminus0p26e6_"
                "msun_distance_68p27_plusminus0p45_mpc_main_text_edge_on_"
                "mestel_assumption_distance_imported_from_literature_exact_"
                "spot_selection_and_error_confidence_semantics_unreported_"
                "ratio_endpoints_hold_reported_distance_fixed"
            ),
        }
        for field, expected in expected_model.items():
            self.assertEqual(check[field], expected, field)
        self.assertAlmostEqual(expected_lower, 0.002434118995564536, places=18)
        self.assertAlmostEqual(expected_upper, 0.0031252885868976762, places=18)
        self.assertAlmostEqual(
            bhs.shadow_size_uas(1_850_000.0, 68.27),
            0.002779703791231106,
            places=18,
        )
        self.assertEqual(row["size_model_check_key"], check["check_key"])
        self.assertEqual(row["size_model_check_type"], check["check_type"])
        self.assertEqual(
            row["size_model_check_data_relation"], check["data_relation"]
        )
        self.assertEqual(
            row["size_model_check_range_basis"], check["range_basis"]
        )
        self.assertEqual(
            row["size_model_check_accepted_configuration_count"], 1
        )
        self.assertIsNone(row["size_model_check_excluded_configuration_ids"])
        self.assertEqual(row["size_model_check_lower_uas"], expected_lower)
        self.assertEqual(row["size_model_check_upper_uas"], expected_upper)
        self.assertEqual(
            row["size_model_check_probability_status"], "non_probabilistic"
        )
        self.assertEqual(
            row["size_model_check_adoption_status"],
            "reported_separately_not_adopted",
        )
        self.assertEqual(
            row["size_model_check_reference_key"], "Deepshikha2024"
        )
        self.assertEqual(
            row["size_model_check_reference_doi"],
            "10.1017/S1743921323002491",
        )

        self.assertEqual(alternate["fit_key"], "villaraos_2022_j0437_gr_ratio_fit")
        self.assertEqual(
            alternate["data_relation"],
            "reanalysis_of_published_gao2017_vlbi_spot_data",
        )
        self.assertEqual(
            alternate["model_basis"], "circular_equatorial_schwarzschild_disk"
        )
        self.assertEqual(alternate["ratio_msun_per_mpc"], 42_500.0)
        self.assertEqual(alternate["ratio_unc_minus_msun_per_mpc"], 1_200.0)
        self.assertEqual(alternate["ratio_unc_plus_msun_per_mpc"], 1_000.0)
        self.assertEqual(alternate["size_uas"], 0.004359589760898525)
        self.assertEqual(
            alternate["size_unc_minus_uas"], 0.00012309429913125248
        )
        self.assertEqual(
            alternate["size_unc_plus_uas"], 0.00010257858260937706
        )
        self.assertEqual(
            alternate["nonadoption_basis"],
            "same_data_reanalysis_conditional_on_fixed_gao2017_inclination_"
            "vertical_center_and_zero_scatter_angle",
        )
        self.assertEqual(
            alternate["provenance_note"],
            "published_table_2_ratio_and_formal_interval_three_parameters_"
            "fixed_i81_y0_minus0p040_zero_scatter_chi2red_1p038_section_"
            "4p1p8_y0_variation_up_to_2_percent_within_formal_interval",
        )
        self.assertEqual(alternate["primary_treatment"], "shared_distance")
        self.assertEqual(alternate["primary_result_status"], "retained")
        self.assertEqual(row["nonadopted_fit_size_uas"], alternate["size_uas"])
        self.assertAlmostEqual(
            100.0 * (row["nonadopted_fit_size_uas"] / row["size_uas"] - 1.0),
            4.020979020979021,
            places=14,
        )

        # Both sidecars are attached after the adopted primary calculation.
        for field in bhs.OUTPUT_FIELDS:
            if (
                field.startswith("size_model_check_")
                or field.startswith("nonadopted_fit_")
                or field == "uncertainty_status"
            ):
                continue
            self.assertEqual(row[field], baseline[field], field)

        mutated_check = dict(check)
        mutated_check["size_lower_uas"] = 0.001
        mutated_check["size_upper_uas"] = 9.9
        candidate = bhs.attach_model_check(dict(baseline), mutated_check)
        mutated_alternate = dict(alternate)
        mutated_alternate["ratio_msun_per_mpc"] = 999_999.0
        mutated_alternate["size_uas"] = bhs.shadow_size_uas(999_999.0, 1.0)
        candidate = bhs.attach_nonadopted_ratio_fit(
            candidate, mutated_alternate
        )
        for field in bhs.OUTPUT_FIELDS:
            if (
                field.startswith("size_model_check_")
                or field.startswith("nonadopted_fit_")
                or field == "uncertainty_status"
            ):
                continue
            self.assertEqual(candidate[field], baseline[field], field)

    def test_j0437_reviewed_inputs_comparison_and_output_fail_closed(self) -> None:
        source_mutations = (
            ("reference_mass_msun", "2860001"),
            ("reference_mass_unc_minus_msun", "200001"),
            ("reference_distance_mpc", "71"),
            ("adopted_distance_mpc", "64.6"),
            ("distance_measurement_unc_minus_mpc", "0.541"),
            ("distance_flow_unc_mpc", "3.38"),
            ("distance_shared_systematic_unc_mpc", "2.62"),
            ("systematic_status", "none_adopted"),
            ("mass_reference_doi", "10.0000/wrong"),
            ("distance_reference_doi", "10.0000/wrong"),
            ("uncertainty_note", "different_note"),
            ("reference_mass_systematic_unc_minus_msun", "1"),
        )
        for field, value in source_mutations:
            with self.subTest(source_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in self.inputs]
                next(item for item in rows if item["source_key"] == "J0437")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "J0437"):
                    bhs.load_source_inputs(self._write_source_inputs(tmpdir, rows))

        with bhs.DEFAULT_MODEL_CHECKS.open(newline="", encoding="utf-8") as handle:
            model_rows = list(csv.DictReader(handle))
        model_mutations = (
            ("check_key", "different_check"),
            ("check_type", "data_selection_component_tracking"),
            ("data_relation", "reanalysis_of_other_data"),
            ("range_basis", "accepted_configuration_medians"),
            ("endpoint_derivation", "adopted_primary_ratio_factors"),
            ("sensitivity_direction", "two_sided_range"),
            ("accepted_configuration_count", "2"),
            ("excluded_configuration_ids", "1"),
            ("lower_configuration_id", "1"),
            ("lower_mass_msun", "1620001"),
            ("lower_distance_mpc", "68.28"),
            ("upper_configuration_id", "1"),
            ("upper_mass_msun", "2080001"),
            ("upper_distance_mpc", "68.28"),
            ("lower_primary_ratio_factor", "1"),
            ("upper_primary_ratio_factor", "1.1"),
            ("assumed_inner_radius_pc", "0.1"),
            ("assumed_turbulent_velocity_kms", "10"),
            ("assumed_surface_density_power_law_index", "-0.9"),
            ("assumed_rotation_profile", "flat_rotation"),
            ("mass_bias_eta", "2"),
            ("assumption_profile", "point_mass_only"),
            ("probability_status", "probabilistic"),
            ("adoption_status", "adopted"),
            ("reference_key", "Gao2016b"),
            ("reference_doi", "10.0000/wrong"),
            ("provenance_note", "different_provenance"),
        )
        for field, value in model_mutations:
            with self.subTest(model_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in model_rows]
                next(item for item in rows if item["source_key"] == "J0437")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "J0437"):
                    bhs.load_model_checks(
                        self._write_model_checks(tmpdir, rows),
                        known_sources=self.source_keys,
                    )

        with bhs.DEFAULT_NONADOPTED_FITS.open(
            newline="", encoding="utf-8"
        ) as handle:
            comparison_rows = list(csv.DictReader(handle))
        comparison_mutations = (
            ("fit_key", "villaraos_2022_j0437_other_fit"),
            ("data_relation", "reanalysis_of_other_data"),
            ("ratio_msun_per_mpc", "42501"),
            ("ratio_unc_minus_msun_per_mpc", "1201"),
            ("ratio_unc_plus_msun_per_mpc", "1001"),
            ("nonadoption_basis", "same_data_reanalysis_conditional_on_fixed_vertical_center"),
            ("reference_doi", "10.0000/wrong"),
            ("provenance_note", "different_provenance"),
        )
        for field, value in comparison_mutations:
            with self.subTest(comparison_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in comparison_rows]
                next(item for item in rows if item["source_key"] == "J0437")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "J0437"):
                    bhs.load_nonadopted_ratio_fits(
                        self._write_nonadopted_fits(tmpdir, rows),
                        source_treatments=self.source_treatments,
                    )

        output_mutations = (
            ("size_uas", 0.004191067803754549 + 1e-12),
            ("size_stat_unc_minus_uas", 0.0002930816645982202 + 1e-12),
            ("mass_msun", 2_635_480.0),
            ("mass_shared_distance_unc_minus_msun", 176_181.0),
            ("mass_total_unc_minus_msun", 254_963.0),
            ("distance_stat_unc_minus_mpc", 3.4261),
            ("distance_unc_minus_mpc", 4.3122),
            ("size_systematic_status", "numeric_systematic"),
            ("size_total_status", "formal_component_only"),
            ("size_total_unc_minus_uas", 0.001),
            ("size_model_check_key", "different_check"),
            ("size_model_check_type", "data_selection_component_tracking"),
            ("size_model_check_lower_uas", 0.0025),
            ("size_model_check_upper_uas", 0.0032),
            ("size_model_check_reference_key", "Gao2016b"),
            ("size_model_check_reference_doi", "10.0000/wrong"),
            ("nonadopted_fit_key", "different_fit"),
            ("nonadopted_fit_ratio_msun_per_mpc", 42_501.0),
            ("nonadopted_fit_size_uas", 0.00436),
            ("nonadopted_fit_basis", "none"),
            ("nonadopted_fit_reference_doi", "10.0000/wrong"),
        )
        for field, value in output_mutations:
            with self.subTest(output_field=field):
                rows = [dict(item) for item in self.rows]
                next(item for item in rows if item["source_key"] == "J0437")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "J0437"):
                    bhs.render_csv(rows)

    def test_mrk1029_no_primary_diagnostics_are_exact_and_separate(self) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "Mrk1029")
        baseline = bhs.build_no_mass_row(source)
        row = self.by_key["Mrk1029"]
        check = self.model_checks["Mrk1029"]
        alternate = self.nonadopted_fits["Mrk1029"]

        expected_source = {
            "source_order": "20",
            "source_label": "Mrk 1029",
            "treatment": "no_robust_mass",
            "reference_mass_msun": "",
            "reference_mass_unc_minus_msun": "",
            "reference_mass_unc_plus_msun": "",
            "reference_distance_mpc": "",
            "adopted_distance_mpc": "119.32882273342354",
            "distance_measurement_unc_minus_mpc": "0.4330175913396481",
            "distance_measurement_unc_plus_mpc": "0.4330175913396481",
            "distance_flow_unc_mpc": "3.3829499323410013",
            "distance_shared_systematic_unc_mpc": "4.8442011935083977",
            "size_stat_fraction_minus_override": "",
            "size_stat_fraction_plus_override": "",
            "size_systematic_fraction_minus": "",
            "size_systematic_fraction_plus": "",
            "systematic_status": "not_applicable",
            "mass_reference_doi": "10.3847/1538-4357/834/1/52",
            "mass_uncertainty_reference_doi": "NA",
            "distance_reference_doi": "10.3847/1538-4357/834/1/52",
            "uncertainty_note": (
                "no_robust_primary_mass_huchra1999_9076_plusminus32_optical_"
                "velocity_component_selected_amid_two_velocity_ambiguity"
            ),
            "mass_display_interval_status": "not_available",
            "distance_stat_component_status": "statistical",
            "reference_mass_systematic_unc_minus_msun": "",
            "reference_mass_systematic_unc_plus_msun": "",
            "mass_systematic_status": "not_available",
        }
        for field, expected in expected_source.items():
            self.assertEqual(source[field], expected, field)

        measurement = 32.0 / 73.9
        flow = 250.0 / 73.9
        distance = 120.8 * 73.0 / 73.9
        distance_stat = math.hypot(measurement, flow)
        distance_systematic = 4.844201193508398
        distance_total = 5.924368293946939
        self.assertEqual(float(source["distance_measurement_unc_minus_mpc"]), measurement)
        self.assertEqual(float(source["distance_measurement_unc_plus_mpc"]), measurement)
        self.assertEqual(float(source["distance_flow_unc_mpc"]), flow)
        self.assertEqual(float(source["adopted_distance_mpc"]), distance)
        self.assertEqual(row["distance_mpc"], distance)
        self.assertEqual(row["distance_stat_unc_minus_mpc"], distance_stat)
        self.assertEqual(row["distance_stat_unc_plus_mpc"], distance_stat)
        self.assertEqual(
            row["distance_systematic_unc_minus_mpc"], distance_systematic
        )
        self.assertEqual(
            row["distance_systematic_unc_plus_mpc"], distance_systematic
        )
        self.assertEqual(row["distance_unc_minus_mpc"], distance_total)
        self.assertEqual(row["distance_unc_plus_mpc"], distance_total)

        expected_model = {
            "check_key": "deepshikha_2024_mrk1029_mestel_disk_mass_model",
            "check_type": "disk_mass_model_choice",
            "data_relation": (
                "reanalysis_of_published_gao2017_vlbi_data_data_selection_"
                "unreported"
            ),
            "range_basis": "reported_symmetric_mass_error_span",
            "endpoint_derivation": (
                "reported_mbh_plusminus_error_at_fixed_reported_distance"
            ),
            "sensitivity_direction": "two_sided_range",
            "accepted_configuration_count": 1,
            "excluded_configuration_ids": None,
            "lower_configuration_id": None,
            "lower_mass_msun": 410_000.0,
            "lower_distance_mpc": 126.787,
            "upper_configuration_id": None,
            "upper_mass_msun": 730_000.0,
            "upper_distance_mpc": 126.787,
            "lower_primary_ratio_factor": None,
            "upper_primary_ratio_factor": None,
            "assumed_inner_radius_pc": None,
            "assumed_turbulent_velocity_kms": None,
            "assumed_surface_density_power_law_index": -1.0,
            "assumed_rotation_profile": None,
            "mass_bias_eta": None,
            "assumption_profile": (
                "edge_on_coplanar_black_hole_plus_mestel_disk"
            ),
            "size_lower_uas": 0.0003317155455200028,
            "size_upper_uas": 0.0005906154834868343,
            "probability_status": "non_probabilistic",
            "adoption_status": "reported_separately_not_adopted",
            "reference_key": "Deepshikha2024",
            "reference_doi": "10.1017/S1743921323002491",
            "provenance_note": (
                "official_supplement_table_1_abs_alpha_0p310_plusminus0p065_"
                "mbh_0p57_plusminus0p16e6_msun_mdisk_0p93_plusminus0p65e6_"
                "msun_distance_126p787_plusminus0p059_mpc_main_text_edge_on_"
                "mestel_assumption_distance_imported_from_literature_exact_"
                "spot_selection_and_error_confidence_semantics_unreported_"
                "ratio_endpoints_hold_reported_distance_fixed_no_adopted_"
                "primary_for_directional_comparison"
            ),
        }
        for field, expected in expected_model.items():
            self.assertEqual(check[field], expected, field)
        self.assertEqual(
            bhs.shadow_size_uas(410_000.0, 126.787),
            0.0003317155455200028,
        )
        self.assertEqual(
            bhs.shadow_size_uas(570_000.0, 126.787),
            0.0004611655145034186,
        )
        self.assertEqual(
            bhs.shadow_size_uas(730_000.0, 126.787),
            0.0005906154834868343,
        )
        self.assertEqual(row["size_model_check_key"], check["check_key"])
        self.assertEqual(row["size_model_check_type"], check["check_type"])
        self.assertEqual(
            row["size_model_check_data_relation"], check["data_relation"]
        )
        self.assertEqual(
            row["size_model_check_range_basis"], check["range_basis"]
        )
        self.assertEqual(row["size_model_check_accepted_configuration_count"], 1)
        self.assertIsNone(row["size_model_check_excluded_configuration_ids"])
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
            alternate["fit_key"], "gonzalez_juarez_2024_mrk1029_gr_ratio_fit"
        )
        self.assertEqual(
            alternate["data_relation"],
            "reanalysis_of_published_gao2017_vlbi_spot_data",
        )
        self.assertEqual(
            alternate["model_basis"], "circular_equatorial_schwarzschild_disk"
        )
        self.assertEqual(alternate["ratio_msun_per_mpc"], 14_400.0)
        self.assertEqual(alternate["ratio_unc_minus_msun_per_mpc"], 1_100.0)
        self.assertEqual(alternate["ratio_unc_plus_msun_per_mpc"], 1_100.0)
        self.assertEqual(alternate["size_uas"], 0.0014771315895750299)
        self.assertEqual(
            alternate["size_unc_minus_uas"], 0.00011283644087031476
        )
        self.assertEqual(
            alternate["size_unc_plus_uas"], 0.00011283644087031476
        )
        self.assertEqual(
            alternate["nonadoption_basis"],
            "restricted_same_data_fit_conditional_on_fixed_gao2017_inclination_"
            "vertical_center_and_selected_scatter_angle_does_not_resolve_primary_"
            "complex_geometry_and_dynamics",
        )
        self.assertEqual(
            alternate["provenance_note"],
            "published_table_2_ratio_and_formal_interval_three_parameters_m_over_"
            "d_zrec_x0_fitted_y0_fixed_i79_fixed_scatter14deg_chi2red1p385_pa_"
            "plusminus10deg_sensitivity_1p19percent_neglected_within_formal_interval",
        )
        self.assertEqual(alternate["primary_treatment"], "no_robust_mass")
        self.assertEqual(alternate["primary_result_status"], "absent")

        for field in (
            "size_uas",
            "size_unc_adopted_uas",
            "mass_msun",
            "mass_fit_unc_minus_msun",
            "mass_fit_unc_plus_msun",
            "mass_total_unc_minus_msun",
            "mass_total_unc_plus_msun",
            "size_stat_unc_minus_uas",
            "size_stat_unc_plus_uas",
            "size_systematic_unc_minus_uas",
            "size_systematic_unc_plus_uas",
            "size_total_unc_minus_uas",
            "size_total_unc_plus_uas",
        ):
            self.assertIsNone(row[field], field)
        self.assertEqual(row["treatment"], "no_robust_mass")
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

        # Both standalone diagnostics attach after the absent-primary row is built.
        for field in bhs.OUTPUT_FIELDS:
            if (
                field.startswith("size_model_check_")
                or field.startswith("nonadopted_fit_")
                or field == "uncertainty_status"
            ):
                continue
            self.assertEqual(row[field], baseline[field], field)

        mutated_check = dict(check)
        mutated_check["size_lower_uas"] = 0.1
        mutated_check["size_upper_uas"] = 0.2
        candidate = bhs.attach_model_check(dict(baseline), mutated_check)
        mutated_alternate = dict(alternate)
        mutated_alternate["ratio_msun_per_mpc"] = 999_999.0
        mutated_alternate["size_uas"] = bhs.shadow_size_uas(999_999.0, 1.0)
        candidate = bhs.attach_nonadopted_ratio_fit(candidate, mutated_alternate)
        for field in bhs.OUTPUT_FIELDS:
            if (
                field.startswith("size_model_check_")
                or field.startswith("nonadopted_fit_")
                or field == "uncertainty_status"
            ):
                continue
            self.assertEqual(candidate[field], baseline[field], field)

    def test_mrk1029_reviewed_inputs_sidecars_and_output_fail_closed(self) -> None:
        source_mutations = (
            ("source_label", "Mrk1029"),
            ("reference_mass_msun", "570000"),
            ("reference_distance_mpc", "126.787"),
            ("adopted_distance_mpc", "119.4"),
            ("distance_measurement_unc_minus_mpc", "0.433"),
            ("distance_measurement_unc_plus_mpc", "0.433"),
            ("distance_flow_unc_mpc", "3.38"),
            ("distance_shared_systematic_unc_mpc", "4.85"),
            ("systematic_status", "none_adopted"),
            ("mass_reference_doi", "10.0000/wrong"),
            ("mass_uncertainty_reference_doi", "10.0000/wrong"),
            ("distance_reference_doi", "10.0000/wrong"),
            ("uncertainty_note", "different_note"),
            ("mass_display_interval_status", "decomposed_components"),
            ("distance_stat_component_status", "not_available"),
            ("mass_systematic_status", "numeric_systematic"),
        )
        for field, value in source_mutations:
            with self.subTest(source_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in self.inputs]
                next(item for item in rows if item["source_key"] == "Mrk1029")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "Mrk1029"):
                    bhs.load_source_inputs(self._write_source_inputs(tmpdir, rows))

        with bhs.DEFAULT_MODEL_CHECKS.open(newline="", encoding="utf-8") as handle:
            model_rows = list(csv.DictReader(handle))
        model_mutations = (
            ("check_key", "different_check"),
            ("check_type", "data_selection_component_tracking"),
            ("data_relation", "reanalysis_of_other_data"),
            ("range_basis", "accepted_configuration_medians"),
            ("endpoint_derivation", "adopted_primary_ratio_factors"),
            ("sensitivity_direction", "downward_only"),
            ("accepted_configuration_count", "2"),
            ("excluded_configuration_ids", "1"),
            ("lower_configuration_id", "1"),
            ("lower_mass_msun", "410001"),
            ("lower_distance_mpc", "126.788"),
            ("upper_configuration_id", "1"),
            ("upper_mass_msun", "730001"),
            ("upper_distance_mpc", "126.788"),
            ("lower_primary_ratio_factor", "1"),
            ("upper_primary_ratio_factor", "1.1"),
            ("assumed_inner_radius_pc", "0.1"),
            ("assumed_turbulent_velocity_kms", "10"),
            ("assumed_surface_density_power_law_index", "-0.9"),
            ("assumed_rotation_profile", "flat_rotation"),
            ("mass_bias_eta", "2"),
            ("assumption_profile", "point_mass_only"),
            ("probability_status", "probabilistic"),
            ("adoption_status", "adopted"),
            ("reference_key", "Gao2016b"),
            ("reference_doi", "10.0000/wrong"),
            ("provenance_note", "different_provenance"),
        )
        for field, value in model_mutations:
            with self.subTest(model_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in model_rows]
                next(item for item in rows if item["source_key"] == "Mrk1029")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "Mrk1029"):
                    bhs.load_model_checks(
                        self._write_model_checks(tmpdir, rows),
                        known_sources=self.source_keys,
                    )

        with bhs.DEFAULT_NONADOPTED_FITS.open(
            newline="", encoding="utf-8"
        ) as handle:
            comparison_rows = list(csv.DictReader(handle))
        comparison_mutations = (
            ("fit_key", "gonzalez_juarez_2024_mrk1029_other_fit"),
            ("data_relation", "reanalysis_of_other_data"),
            ("model_basis", "unknown_model"),
            ("ratio_msun_per_mpc", "14401"),
            ("ratio_unc_minus_msun_per_mpc", "1101"),
            ("ratio_unc_plus_msun_per_mpc", "1101"),
            ("interval_status", "unconditional"),
            ("adoption_status", "adopted"),
            ("nonadoption_basis", "does_not_resolve_primary_complex_geometry_and_dynamics"),
            ("reference_key", "Villaraos2022"),
            ("reference_doi", "10.0000/wrong"),
            ("provenance_note", "published_table_2_ratio_and_formal_interval"),
        )
        for field, value in comparison_mutations:
            with self.subTest(comparison_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in comparison_rows]
                next(item for item in rows if item["source_key"] == "Mrk1029")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "Mrk1029"):
                    bhs.load_nonadopted_ratio_fits(
                        self._write_nonadopted_fits(tmpdir, rows),
                        source_treatments=self.source_treatments,
                    )

        output_mutations = (
            ("distance_mpc", 119.4),
            ("distance_stat_unc_minus_mpc", 3.4106),
            ("distance_stat_unc_plus_mpc", 3.4106),
            ("distance_systematic_unc_minus_mpc", 4.8443),
            ("distance_unc_minus_mpc", 5.9244),
            ("size_uas", 0.0004),
            ("mass_msun", 570_000.0),
            ("size_stat_component_status", "fixed_distance_formal_fit_interval"),
            ("size_systematic_status", "unquantified_model"),
            ("size_total_status", "incomplete_unquantified_systematic"),
            ("uncertainty_status", "different_status"),
            ("size_model_check_key", "different_check"),
            ("size_model_check_type", "data_selection_component_tracking"),
            ("size_model_check_data_relation", "reanalysis_of_other_data"),
            ("size_model_check_range_basis", "accepted_configuration_medians"),
            ("size_model_check_accepted_configuration_count", 2.0),
            ("size_model_check_excluded_configuration_ids", "1"),
            ("size_model_check_lower_uas", 0.00034),
            ("size_model_check_upper_uas", 0.00060),
            ("size_model_check_probability_status", "probabilistic"),
            ("size_model_check_adoption_status", "adopted"),
            ("size_model_check_reference_key", "Gao2016b"),
            ("size_model_check_reference_doi", "10.0000/wrong"),
            ("nonadopted_fit_key", "different_fit"),
            ("nonadopted_fit_ratio_msun_per_mpc", 14_401.0),
            ("nonadopted_fit_ratio_unc_minus_msun_per_mpc", 1_101.0),
            ("nonadopted_fit_size_uas", 0.00148),
            ("nonadopted_fit_basis", "none"),
            ("nonadopted_fit_reference_key", "Villaraos2022"),
            ("nonadopted_fit_reference_doi", "10.0000/wrong"),
        )
        for field, value in output_mutations:
            with self.subTest(output_field=field):
                rows = [dict(item) for item in self.rows]
                next(item for item in rows if item["source_key"] == "Mrk1029")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "Mrk1029"):
                    bhs.render_csv(rows)

    def test_ngc1320_no_primary_diagnostics_and_audit_only_check_are_exact(
        self,
    ) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "NGC1320")
        baseline = bhs.build_no_mass_row(source)
        row = self.by_key["NGC1320"]
        check = self.model_checks["NGC1320"]
        alternate = self.nonadopted_fits["NGC1320"]
        audit_only = self.audit_only_model_checks[
            "kuo_2020_ngc1320_turbulent_support_bias"
        ]

        expected_source = {
            "source_order": "21",
            "source_key": "NGC1320",
            "source_label": "NGC 1320",
            "treatment": "no_robust_mass",
            "reference_mass_msun": "",
            "reference_mass_unc_minus_msun": "",
            "reference_mass_unc_plus_msun": "",
            "reference_distance_mpc": "",
            "adopted_distance_mpc": "33.783491204330176",
            "distance_measurement_unc_minus_mpc": "0.21650879566982406",
            "distance_measurement_unc_plus_mpc": "0.21650879566982406",
            "distance_flow_unc_mpc": "3.3829499323410013",
            "distance_shared_systematic_unc_mpc": "1.3714543114071789",
            "size_stat_fraction_minus_override": "",
            "size_stat_fraction_plus_override": "",
            "size_systematic_fraction_minus": "",
            "size_systematic_fraction_plus": "",
            "systematic_status": "not_applicable",
            "mass_reference_doi": "10.3847/1538-4357/834/1/52",
            "mass_uncertainty_reference_doi": "NA",
            "distance_reference_doi": "10.3847/1538-4357/834/1/52",
            "uncertainty_note": (
                "no_robust_primary_mass_gao2017_2649_plusminus16_optical_"
                "velocity_component"
            ),
            "mass_display_interval_status": "not_available",
            "distance_stat_component_status": "statistical",
            "reference_mass_systematic_unc_minus_msun": "",
            "reference_mass_systematic_unc_plus_msun": "",
            "mass_systematic_status": "not_available",
        }
        self.assertEqual(source, expected_source)

        measurement = 16.0 / 73.9
        flow = 250.0 / 73.9
        distance = 34.2 * 73.0 / 73.9
        distance_stat = math.hypot(measurement, flow)
        distance_systematic = distance * 3.0 / 73.9
        distance_total = math.hypot(distance_stat, distance_systematic)
        self.assertEqual(float(source["distance_measurement_unc_minus_mpc"]), measurement)
        self.assertEqual(float(source["distance_measurement_unc_plus_mpc"]), measurement)
        self.assertEqual(float(source["distance_flow_unc_mpc"]), flow)
        self.assertEqual(float(source["adopted_distance_mpc"]), distance)
        self.assertEqual(
            float(source["distance_shared_systematic_unc_mpc"]),
            distance_systematic,
        )
        self.assertEqual(row["distance_mpc"], distance)
        self.assertEqual(row["distance_stat_unc_minus_mpc"], distance_stat)
        self.assertEqual(row["distance_stat_unc_plus_mpc"], distance_stat)
        self.assertEqual(
            row["distance_systematic_unc_minus_mpc"], distance_systematic
        )
        self.assertEqual(
            row["distance_systematic_unc_plus_mpc"], distance_systematic
        )
        self.assertEqual(row["distance_unc_minus_mpc"], distance_total)
        self.assertEqual(row["distance_unc_plus_mpc"], distance_total)
        self.assertEqual(distance_stat, 3.389871133734789)
        self.assertEqual(distance_total, 3.656790017434105)

        expected_model = {
            "source_key": "NGC1320",
            "check_key": "deepshikha_2024_ngc1320_mestel_disk_mass_model",
            "check_type": "disk_mass_model_choice",
            "data_relation": (
                "reanalysis_of_published_gao2017_vlbi_data_data_selection_"
                "unreported"
            ),
            "range_basis": "reported_symmetric_mass_error_span",
            "endpoint_derivation": (
                "reported_mbh_plusminus_error_at_fixed_reported_distance"
            ),
            "sensitivity_direction": "two_sided_range",
            "accepted_configuration_count": 1,
            "excluded_configuration_ids": None,
            "lower_configuration_id": None,
            "lower_mass_msun": 4_130_000.0,
            "lower_distance_mpc": 38.4,
            "upper_configuration_id": None,
            "upper_mass_msun": 5_290_000.0,
            "upper_distance_mpc": 38.4,
            "lower_primary_ratio_factor": None,
            "upper_primary_ratio_factor": None,
            "assumed_inner_radius_pc": None,
            "assumed_turbulent_velocity_kms": None,
            "assumed_surface_density_power_law_index": -1.0,
            "assumed_rotation_profile": None,
            "mass_bias_eta": None,
            "assumption_profile": "edge_on_coplanar_black_hole_plus_mestel_disk",
            "size_lower_uas": 0.01103254026501894,
            "size_upper_uas": 0.014131268281343874,
            "probability_status": "non_probabilistic",
            "adoption_status": "reported_separately_not_adopted",
            "reference_key": "Deepshikha2024",
            "reference_doi": "10.1017/S1743921323002491",
            "provenance_note": (
                "official_supplement_table_1_abs_alpha_0p479_plusminus0p078_"
                "mbh_4p71_plusminus0p58e6_msun_mdisk_0p10_plusminus1p14e6_"
                "msun_distance_38p4_plusminus0p44_mpc_main_text_edge_on_mestel_"
                "assumption_distance_imported_from_literature_exact_spot_"
                "selection_and_error_confidence_semantics_unreported_ratio_"
                "endpoints_hold_reported_distance_fixed_no_adopted_primary_for_"
                "directional_comparison"
            ),
        }
        self.assertEqual(check, expected_model)
        self.assertEqual(
            bhs.shadow_size_uas(4_130_000.0, 38.4),
            0.01103254026501894,
        )
        self.assertEqual(
            bhs.shadow_size_uas(4_710_000.0, 38.4),
            0.012581904273181408,
        )
        self.assertEqual(
            bhs.shadow_size_uas(5_290_000.0, 38.4),
            0.014131268281343874,
        )

        expected_alternate = {
            "source_key": "NGC1320",
            "fit_key": "gonzalez_juarez_2024_ngc1320_gr_ratio_fit",
            "data_relation": "reanalysis_of_published_gao2017_vlbi_spot_data",
            "model_basis": "circular_equatorial_schwarzschild_disk",
            "ratio_msun_per_mpc": 149_700.0,
            "ratio_unc_minus_msun_per_mpc": 6_900.0,
            "ratio_unc_plus_msun_per_mpc": 7_100.0,
            "size_uas": 0.015356013816623746,
            "size_unc_minus_uas": 0.0007077922200047018,
            "size_unc_plus_uas": 0.0007283079365265772,
            "interval_status": "bayesian_formal_conditional_on_model",
            "adoption_status": "reported_separately_not_adopted",
            "nonadoption_basis": (
                "restricted_same_data_fit_conditional_on_fixed_edge_on_"
                "inclination_and_selected_scatter_angle_does_not_resolve_"
                "primary_complex_geometry_and_dynamics"
            ),
            "reference_key": "GonzalezJuarez2024",
            "reference_doi": "10.1051/0004-6361/202450098",
            "provenance_note": (
                "published_table_2_ratio_and_formal_interval_four_parameters_"
                "m_over_d_zrec_x0_y0_fitted_reference_midpoint_between_"
                "redshifted_and_blueshifted_group_geometric_centers_i90_fixed_"
                "scatter17deg_chi2red1p127_pa_plusminus10deg_sensitivity_"
                "2p64percent_within_formal_interval_bh_position_near_"
                "blueshifted_group_new_vlbi_requested"
            ),
            "primary_treatment": "no_robust_mass",
            "primary_result_status": "absent",
        }
        self.assertEqual(alternate, expected_alternate)

        expected_audit_only = {
            "check_key": "kuo_2020_ngc1320_turbulent_support_bias",
            "source_key": "NGC1320",
            "check_type": "turbulent_support_mass_bias",
            "data_relation": (
                "published_analytic_model_applied_to_ngc1320_indicative_disk_"
                "parameters"
            ),
            "range_basis": "one_sided_fractional_mass_underestimate",
            "endpoint_derivation": "conditional_reference_ratio_factors",
            "sensitivity_direction": "upward_only",
            "assumed_rotation_profile": "keplerian_rotation",
            "assumption_profile": (
                "constant_turbulent_velocity_and_turbulent_pressure_dominates_"
                "thermal"
            ),
            "probability_status": "non_probabilistic",
            "adoption_status": "reported_separately_not_adopted",
            "attachment_status": "audit_only_unattached_no_robust_primary",
            "reference_key": "Kuo2020",
            "reference_doi": "10.1093/mnras/staa2260",
            "provenance_note": (
                "appendix_a_ngc1320_rin_0p07pc_vturb_40kms_sigma_r_minus1_"
                "keplerian_eta_2p5_delta_m_turb_over_m_enc_0p01_conditional_"
                "on_gao2017_indicative_enclosed_mass_no_adopted_primary"
            ),
            "lower_reference_ratio_factor": 1.0,
            "upper_reference_ratio_factor": 1.01,
            "assumed_inner_radius_pc": 0.07,
            "assumed_turbulent_velocity_kms": 40.0,
            "assumed_surface_density_power_law_index": -1.0,
            "mass_bias_eta": 2.5,
        }
        self.assertEqual(len(self.audit_only_model_checks), 1)
        self.assertEqual(audit_only, expected_audit_only)

        expected_attached = {
            "size_model_check_key": check["check_key"],
            "size_model_check_type": check["check_type"],
            "size_model_check_data_relation": check["data_relation"],
            "size_model_check_range_basis": check["range_basis"],
            "size_model_check_accepted_configuration_count": 1,
            "size_model_check_excluded_configuration_ids": None,
            "size_model_check_lower_uas": check["size_lower_uas"],
            "size_model_check_upper_uas": check["size_upper_uas"],
            "size_model_check_probability_status": "non_probabilistic",
            "size_model_check_adoption_status": "reported_separately_not_adopted",
            "size_model_check_reference_key": "Deepshikha2024",
            "size_model_check_reference_doi": "10.1017/S1743921323002491",
            "nonadopted_fit_key": alternate["fit_key"],
            "nonadopted_fit_ratio_msun_per_mpc": 149_700.0,
            "nonadopted_fit_ratio_unc_minus_msun_per_mpc": 6_900.0,
            "nonadopted_fit_ratio_unc_plus_msun_per_mpc": 7_100.0,
            "nonadopted_fit_size_uas": 0.015356013816623746,
            "nonadopted_fit_size_unc_minus_uas": 0.0007077922200047018,
            "nonadopted_fit_size_unc_plus_uas": 0.0007283079365265772,
            "nonadopted_fit_interval_status": "bayesian_formal_conditional_on_model",
            "nonadopted_fit_adoption_status": "reported_separately_not_adopted",
            "nonadopted_fit_basis": alternate["nonadoption_basis"],
            "nonadopted_fit_reference_key": "GonzalezJuarez2024",
            "nonadopted_fit_reference_doi": "10.1051/0004-6361/202450098",
        }
        for field, expected in expected_attached.items():
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
            self.assertIsNone(row[field], field)
        self.assertEqual(row["treatment"], "no_robust_mass")
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
        self.assertNotEqual(row["size_model_check_key"], audit_only["check_key"])
        self.assertNotIn(audit_only["check_key"], bhs.render_csv(self.rows))

        # Attached and audit-only diagnostics cannot manufacture a primary result.
        for field in bhs.OUTPUT_FIELDS:
            if (
                field.startswith("size_model_check_")
                or field.startswith("nonadopted_fit_")
                or field == "uncertainty_status"
            ):
                continue
            self.assertEqual(row[field], baseline[field], field)

        mutated_check = dict(check)
        mutated_check["size_lower_uas"] = 0.1
        mutated_check["size_upper_uas"] = 0.2
        candidate = bhs.attach_model_check(dict(baseline), mutated_check)
        mutated_alternate = dict(alternate)
        mutated_alternate["ratio_msun_per_mpc"] = 999_999.0
        mutated_alternate["size_uas"] = bhs.shadow_size_uas(999_999.0, 1.0)
        candidate = bhs.attach_nonadopted_ratio_fit(candidate, mutated_alternate)
        for field in bhs.OUTPUT_FIELDS:
            if (
                field.startswith("size_model_check_")
                or field.startswith("nonadopted_fit_")
                or field == "uncertainty_status"
            ):
                continue
            self.assertEqual(candidate[field], baseline[field], field)

    def test_ngc1320_reviewed_inputs_sidecars_and_output_fail_closed(self) -> None:
        source_mutations = (
            ("source_order", "20"),
            ("source_label", "NGC1320"),
            ("treatment", "shared_distance"),
            ("reference_mass_msun", "4710000"),
            ("reference_mass_unc_minus_msun", "580000"),
            ("reference_mass_unc_plus_msun", "580000"),
            ("reference_distance_mpc", "38.4"),
            ("adopted_distance_mpc", "33.8"),
            ("distance_measurement_unc_minus_mpc", "0.217"),
            ("distance_measurement_unc_plus_mpc", "0.217"),
            ("distance_flow_unc_mpc", "3.38"),
            ("distance_shared_systematic_unc_mpc", "1.37"),
            ("size_stat_fraction_minus_override", "0.1"),
            ("size_stat_fraction_plus_override", "0.1"),
            ("size_systematic_fraction_minus", "0.1"),
            ("size_systematic_fraction_plus", "0.1"),
            ("systematic_status", "none_adopted"),
            ("mass_reference_doi", "10.0000/wrong"),
            ("mass_uncertainty_reference_doi", "10.0000/wrong"),
            ("distance_reference_doi", "10.0000/wrong"),
            ("uncertainty_note", "different_note"),
            ("mass_display_interval_status", "decomposed_components"),
            ("distance_stat_component_status", "not_available"),
            ("reference_mass_systematic_unc_minus_msun", "1"),
            ("reference_mass_systematic_unc_plus_msun", "1"),
            ("mass_systematic_status", "numeric_systematic"),
        )
        for field, value in source_mutations:
            with self.subTest(source_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in self.inputs]
                next(item for item in rows if item["source_key"] == "NGC1320")[
                    field
                ] = value
                with self.assertRaises((ValueError, StopIteration)):
                    bhs.load_source_inputs(self._write_source_inputs(tmpdir, rows))

        with bhs.DEFAULT_MODEL_CHECKS.open(newline="", encoding="utf-8") as handle:
            model_rows = list(csv.DictReader(handle))
        model_mutations = {
            "source_key": "UNKNOWN",
            "check_key": "different_check",
            "check_type": "turbulent_support_mass_bias",
            "data_relation": "reanalysis_of_other_data",
            "range_basis": "accepted_configuration_medians",
            "endpoint_derivation": "adopted_primary_ratio_factors",
            "sensitivity_direction": "downward_only",
            "accepted_configuration_count": "2",
            "excluded_configuration_ids": "1",
            "lower_configuration_id": "1",
            "lower_mass_msun": "4130001",
            "lower_distance_mpc": "38.5",
            "upper_configuration_id": "1",
            "upper_mass_msun": "5290001",
            "upper_distance_mpc": "38.5",
            "lower_primary_ratio_factor": "1",
            "upper_primary_ratio_factor": "1.1",
            "assumed_inner_radius_pc": "0.1",
            "assumed_turbulent_velocity_kms": "10",
            "assumed_surface_density_power_law_index": "-0.9",
            "assumed_rotation_profile": "flat_rotation",
            "mass_bias_eta": "2",
            "assumption_profile": (
                "constant_turbulent_velocity_and_turbulent_pressure_dominates_"
                "thermal"
            ),
            "probability_status": "probabilistic",
            "adoption_status": "adopted",
            "reference_key": "Kuo2020",
            "reference_doi": "10.0000/wrong",
            "provenance_note": "different_provenance",
        }
        self.assertEqual(set(model_mutations), set(bhs.MODEL_CHECK_FIELDS))
        for field, value in model_mutations.items():
            with self.subTest(model_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in model_rows]
                next(item for item in rows if item["source_key"] == "NGC1320")[
                    field
                ] = value
                with self.assertRaises(ValueError):
                    bhs.load_model_checks(
                        self._write_model_checks(tmpdir, rows),
                        known_sources=self.source_keys,
                    )

        with bhs.DEFAULT_NONADOPTED_FITS.open(
            newline="", encoding="utf-8"
        ) as handle:
            comparison_rows = list(csv.DictReader(handle))
        comparison_mutations = {
            "source_key": "UNKNOWN",
            "fit_key": "different_fit",
            "data_relation": "reanalysis_of_other_data",
            "model_basis": "unknown_model",
            "ratio_msun_per_mpc": "149701",
            "ratio_unc_minus_msun_per_mpc": "6901",
            "ratio_unc_plus_msun_per_mpc": "7101",
            "interval_status": "unconditional",
            "adoption_status": "adopted",
            "nonadoption_basis": "none",
            "reference_key": "Villaraos2022",
            "reference_doi": "10.0000/wrong",
            "provenance_note": "different_provenance",
        }
        self.assertEqual(set(comparison_mutations), set(bhs.NONADOPTED_FIT_FIELDS))
        for field, value in comparison_mutations.items():
            with self.subTest(comparison_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in comparison_rows]
                next(item for item in rows if item["source_key"] == "NGC1320")[
                    field
                ] = value
                with self.assertRaises(ValueError):
                    bhs.load_nonadopted_ratio_fits(
                        self._write_nonadopted_fits(tmpdir, rows),
                        source_treatments=self.source_treatments,
                    )

        with bhs.DEFAULT_AUDIT_ONLY_MODEL_CHECKS.open(
            newline="", encoding="utf-8"
        ) as handle:
            audit_only_rows = list(csv.DictReader(handle))
        audit_only_mutations = {
            "source_key": "UNKNOWN",
            "check_key": "different_check",
            "check_type": "disk_mass_model_choice",
            "data_relation": "published_analytic_model_applied_to_other_parameters",
            "range_basis": "reported_symmetric_mass_error_span",
            "endpoint_derivation": "reported_mbh_plusminus_error",
            "sensitivity_direction": "two_sided_range",
            "lower_reference_ratio_factor": "1.001",
            "upper_reference_ratio_factor": "1.02",
            "assumed_inner_radius_pc": "0.08",
            "assumed_turbulent_velocity_kms": "41",
            "assumed_surface_density_power_law_index": "-0.9",
            "assumed_rotation_profile": "flat_rotation",
            "mass_bias_eta": "2.0",
            "assumption_profile": "different_assumptions",
            "probability_status": "probabilistic",
            "adoption_status": "adopted",
            "attachment_status": "attached",
            "reference_key": "Deepshikha2024",
            "reference_doi": "10.0000/wrong",
            "provenance_note": "different_provenance",
        }
        self.assertEqual(
            set(audit_only_mutations), set(bhs.AUDIT_ONLY_MODEL_CHECK_FIELDS)
        )
        for field, value in audit_only_mutations.items():
            with self.subTest(audit_only_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in audit_only_rows]
                rows[0][field] = value
                with self.assertRaises(ValueError):
                    bhs.load_audit_only_model_checks(
                        self._write_audit_only_model_checks(tmpdir, rows),
                        known_sources=self.source_keys,
                    )

        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ValueError, "exactly match"):
                bhs.load_audit_only_model_checks(
                    self._write_audit_only_model_checks(tmpdir, []),
                    known_sources=self.source_keys,
                )
        with tempfile.TemporaryDirectory() as tmpdir:
            duplicate = [dict(audit_only_rows[0]), dict(audit_only_rows[0])]
            with self.assertRaisesRegex(ValueError, "multiple audit-only"):
                bhs.load_audit_only_model_checks(
                    self._write_audit_only_model_checks(tmpdir, duplicate),
                    known_sources=self.source_keys,
                )
        with tempfile.TemporaryDirectory() as tmpdir:
            fields = bhs.AUDIT_ONLY_MODEL_CHECK_FIELDS[:-1]
            row = {field: audit_only_rows[0][field] for field in fields}
            with self.assertRaisesRegex(ValueError, "unexpected columns"):
                bhs.load_audit_only_model_checks(
                    self._write_audit_only_model_checks(
                        tmpdir, [row], fields=fields
                    ),
                    known_sources=self.source_keys,
                )

        output_mutations = (
            ("distance_mpc", 33.8),
            ("distance_stat_unc_minus_mpc", 3.39),
            ("distance_stat_unc_plus_mpc", 3.39),
            ("distance_systematic_unc_minus_mpc", 1.38),
            ("distance_systematic_unc_plus_mpc", 1.38),
            ("distance_unc_minus_mpc", 3.66),
            ("distance_unc_plus_mpc", 3.66),
            ("size_uas", 0.012),
            ("mass_msun", 4_710_000.0),
            ("treatment", "shared_distance"),
            ("size_stat_component_status", "fixed_distance_formal_fit_interval"),
            ("size_systematic_status", "unquantified_model"),
            ("size_total_status", "incomplete_unquantified_systematic"),
            ("uncertainty_status", "different_status"),
            ("size_model_check_key", "different_check"),
            ("size_model_check_type", "turbulent_support_mass_bias"),
            ("size_model_check_data_relation", "reanalysis_of_other_data"),
            ("size_model_check_range_basis", "accepted_configuration_medians"),
            ("size_model_check_accepted_configuration_count", 2.0),
            ("size_model_check_excluded_configuration_ids", "1"),
            ("size_model_check_lower_uas", 0.0111),
            ("size_model_check_upper_uas", 0.0142),
            ("size_model_check_probability_status", "probabilistic"),
            ("size_model_check_adoption_status", "adopted"),
            ("size_model_check_reference_key", "Kuo2020"),
            ("size_model_check_reference_doi", "10.0000/wrong"),
            ("nonadopted_fit_key", "different_fit"),
            ("nonadopted_fit_ratio_msun_per_mpc", 149_701.0),
            ("nonadopted_fit_ratio_unc_minus_msun_per_mpc", 6_901.0),
            ("nonadopted_fit_ratio_unc_plus_msun_per_mpc", 7_101.0),
            ("nonadopted_fit_size_uas", 0.0154),
            ("nonadopted_fit_size_unc_minus_uas", 0.0008),
            ("nonadopted_fit_size_unc_plus_uas", 0.0008),
            ("nonadopted_fit_basis", "none"),
            ("nonadopted_fit_reference_key", "Villaraos2022"),
            ("nonadopted_fit_reference_doi", "10.0000/wrong"),
        )
        for field, value in output_mutations:
            with self.subTest(output_field=field):
                rows = [dict(item) for item in self.rows]
                next(item for item in rows if item["source_key"] == "NGC1320")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "NGC1320"):
                    bhs.render_csv(rows)

    def test_ngc3079_reviewed_primary_contract_is_exact_and_fails_closed(self) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "NGC3079")
        row = self.by_key["NGC3079"]

        expected_source = {
            "source_order": "17",
            "source_label": "NGC 3079",
            "treatment": "shared_distance",
            "reference_mass_msun": "2000000",
            "reference_mass_unc_minus_msun": "",
            "reference_mass_unc_plus_msun": "",
            "reference_distance_mpc": "17.3",
            "adopted_distance_mpc": "17.3",
            "distance_measurement_unc_minus_mpc": "1.5",
            "distance_measurement_unc_plus_mpc": "1.5",
            "distance_flow_unc_mpc": "0",
            "distance_shared_systematic_unc_mpc": "0",
            "size_stat_fraction_minus_override": "",
            "size_stat_fraction_plus_override": "",
            "size_systematic_fraction_minus": "",
            "size_systematic_fraction_plus": "",
            "systematic_status": "unquantified_model",
            "mass_reference_doi": "10.1086/426101",
            "mass_uncertainty_reference_doi": "NA",
            "distance_reference_doi": "10.1086/426101",
            "uncertainty_note": (
                "approximate_enclosed_mass_black_hole_proxy_lower_limit_factor_of_"
                "order_unity_flat_rotation_disk_mass_significant_no_published_"
                "formal_interval"
            ),
            "mass_display_interval_status": "not_available",
            "distance_stat_component_status": "reported_legacy_tully_fisher",
            "reference_mass_systematic_unc_minus_msun": "",
            "reference_mass_systematic_unc_plus_msun": "",
            "mass_systematic_status": "not_available",
        }
        for field, expected in expected_source.items():
            self.assertEqual(source[field], expected, field)

        self.assertEqual(row["size_uas"], 0.011858795677384631)
        self.assertEqual(row["mass_msun"], 2_000_000.0)
        self.assertEqual(row["log_mass"], 6.301029995663981)
        self.assertEqual(row["distance_mpc"], 17.3)
        self.assertEqual(row["distance_stat_unc_minus_mpc"], 1.5)
        self.assertEqual(row["distance_stat_unc_plus_mpc"], 1.5)
        self.assertEqual(row["distance_unc_minus_mpc"], 1.5)
        self.assertEqual(row["distance_unc_plus_mpc"], 1.5)
        self.assertEqual(
            row["distance_stat_component_status"],
            "reported_legacy_tully_fisher",
        )
        self.assertIsNone(row["mass_fit_unc_minus_msun"])
        self.assertIsNone(row["mass_fit_unc_plus_msun"])
        self.assertEqual(
            row["mass_shared_distance_unc_minus_msun"],
            173_410.40462427746,
        )
        self.assertEqual(
            row["mass_shared_distance_unc_plus_msun"],
            173_410.40462427746,
        )
        for field in (
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
            self.assertIsNone(row[field], field)
        self.assertEqual(row["mass_display_interval_status"], "not_available")
        self.assertEqual(row["size_stat_component_status"], "no_formal_uncertainty")
        self.assertEqual(row["size_systematic_status"], "unquantified_model")
        self.assertEqual(
            row["size_total_status"], "incomplete_unquantified_systematic"
        )
        self.assertEqual(row["legacy_uncertainty_basis"], "not_available")
        self.assertIn("no formal fit interval", row["uncertainty_status"])
        self.assertIn("model systematic unquantified", row["uncertainty_status"])

        source_mutations = (
            ("reference_mass_msun", "2000001"),
            ("reference_mass_unc_minus_msun", "1"),
            ("reference_mass_unc_plus_msun", "1"),
            ("reference_distance_mpc", "17.4"),
            ("adopted_distance_mpc", "17.4"),
            ("distance_measurement_unc_minus_mpc", "1.6"),
            ("distance_measurement_unc_plus_mpc", "1.6"),
            ("distance_flow_unc_mpc", "1"),
            ("distance_shared_systematic_unc_mpc", "1"),
            ("size_stat_fraction_minus_override", "0.1"),
            ("size_systematic_fraction_plus", "0.1"),
            ("systematic_status", "no_formal_uncertainty"),
            ("mass_reference_doi", "10.0000/wrong"),
            ("mass_uncertainty_reference_doi", "10.0000/wrong"),
            ("distance_reference_doi", "10.0000/wrong"),
            ("uncertainty_note", "different_note"),
            ("mass_display_interval_status", "decomposed_components"),
            ("distance_stat_component_status", "statistical"),
            ("reference_mass_systematic_unc_minus_msun", "1"),
            ("mass_systematic_status", "source_reported_model_selection_systematic"),
        )
        for field, value in source_mutations:
            with self.subTest(source_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in self.inputs]
                next(item for item in rows if item["source_key"] == "NGC3079")[
                    field
                ] = value
                path = self._write_source_inputs(tmpdir, rows)
                with self.assertRaisesRegex(ValueError, "NGC3079"):
                    bhs.load_source_inputs(path)

        output_mutations = (
            ("size_uas", 0.011858795677384631 + 1e-12),
            ("mass_msun", 2_000_001.0),
            ("distance_stat_component_status", "statistical"),
            ("mass_fit_unc_minus_msun", 1.0),
            ("mass_shared_distance_unc_plus_msun", 173_411.0),
            ("mass_total_unc_minus_msun", 173_410.40462427746),
            ("mass_display_unc_plus_msun", 173_410.40462427746),
            ("size_stat_unc_minus_uas", 1e-4),
            ("size_systematic_unc_plus_uas", 1e-4),
            ("size_total_unc_plus_uas", 1e-4),
            ("size_stat_component_status", "fixed_distance_formal_fit_interval"),
            ("size_systematic_status", "no_formal_uncertainty"),
            ("size_total_status", "no_formal_uncertainty"),
        )
        for field, value in output_mutations:
            with self.subTest(output_field=field):
                rows = [dict(item) for item in self.rows]
                next(item for item in rows if item["source_key"] == "NGC3079")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "NGC3079"):
                    bhs.render_csv(rows)

    def test_ngc3079_model_check_is_one_sided_separate_and_reproducible(self) -> None:
        check = self.model_checks["NGC3079"]
        row = self.by_key["NGC3079"]
        expected_lower = 0.011858795677384631
        expected_upper = expected_lower * 1.34

        self.assertEqual(check["check_type"], "turbulent_support_mass_bias")
        self.assertEqual(
            check["data_relation"],
            "published_analytic_model_applied_to_ngc3079_disk_parameters",
        )
        self.assertEqual(
            check["range_basis"], "one_sided_fractional_mass_underestimate"
        )
        self.assertEqual(check["endpoint_derivation"], "adopted_primary_ratio_factors")
        self.assertEqual(check["sensitivity_direction"], "upward_only")
        self.assertIsNone(check["accepted_configuration_count"])
        self.assertIsNone(check["excluded_configuration_ids"])
        self.assertEqual(check["lower_primary_ratio_factor"], 1.0)
        self.assertEqual(check["upper_primary_ratio_factor"], 1.34)
        self.assertEqual(check["assumed_inner_radius_pc"], 0.40)
        self.assertEqual(check["assumed_turbulent_velocity_kms"], 60.0)
        self.assertEqual(check["assumed_surface_density_power_law_index"], -1.0)
        self.assertEqual(check["assumed_rotation_profile"], "flat_rotation")
        self.assertEqual(check["mass_bias_eta"], 2.0)
        self.assertEqual(
            check["assumption_profile"],
            "constant_turbulent_velocity_and_turbulent_pressure_dominates_thermal",
        )
        self.assertEqual(check["probability_status"], "non_probabilistic")
        self.assertEqual(check["adoption_status"], "reported_separately_not_adopted")
        self.assertEqual(check["reference_key"], "Kuo2020")
        self.assertEqual(check["reference_doi"], "10.1093/mnras/staa2260")

        self.assertEqual(row["size_model_check_lower_uas"], expected_lower)
        self.assertEqual(row["size_model_check_upper_uas"], expected_upper)
        self.assertEqual(expected_upper, 0.015890786207695407)
        self.assertEqual(expected_upper - expected_lower, 0.004031990530310775)

        # The model check is not an adopted error, systematic, or quadrature term.
        self.assertEqual(row["size_uas"], expected_lower)
        self.assertEqual(row["mass_msun"], 2_000_000.0)
        self.assertIsNone(row["size_unc_adopted_uas"])
        self.assertIsNone(row["size_stat_unc_minus_uas"])
        self.assertIsNone(row["size_stat_unc_plus_uas"])
        self.assertIsNone(row["size_systematic_unc_minus_uas"])
        self.assertIsNone(row["size_systematic_unc_plus_uas"])
        self.assertIsNone(row["size_total_unc_minus_uas"])
        self.assertIsNone(row["size_total_unc_plus_uas"])
        self.assertEqual(row["size_stat_component_status"], "no_formal_uncertainty")
        self.assertEqual(row["size_systematic_status"], "unquantified_model")
        self.assertEqual(
            row["size_total_status"], "incomplete_unquantified_systematic"
        )
        self.assertEqual(row["legacy_uncertainty_basis"], "not_available")
        self.assertIn("one-sided model sensitivity not combined", row["uncertainty_status"])

    def test_ngc4388_reviewed_components_and_separate_checks_are_exact(self) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "NGC4388")
        row = self.by_key["NGC4388"]
        check = self.model_checks["NGC4388"]
        alternate = self.nonadopted_fits["NGC4388"]

        self.assertEqual(source["reference_mass_msun"], "8400000")
        self.assertEqual(source["reference_distance_mpc"], "19")
        self.assertEqual(source["distance_measurement_unc_minus_mpc"], "")
        self.assertEqual(source["distance_measurement_unc_plus_mpc"], "")
        self.assertEqual(source["distance_shared_systematic_unc_mpc"], "2.09")
        self.assertEqual(source["distance_stat_component_status"], "not_available")
        self.assertEqual(source["uncertainty_note"], bhs.NGC4388_UNCERTAINTY_NOTE)

        self.assertEqual(row["size_uas"], 0.04535053125888249)
        self.assertEqual(row["mass_msun"], 8_400_000.0)
        self.assertEqual(row["distance_mpc"], 19.0)
        self.assertEqual(row["size_stat_unc_minus_uas"], 0.0010797745537829164)
        self.assertEqual(row["size_stat_unc_plus_uas"], 0.0010797745537829164)
        self.assertEqual(
            row["size_stat_component_status"],
            "fit_plus_absolute_position_mixed_component",
        )
        self.assertIsNone(row["size_total_unc_minus_uas"])
        self.assertIsNone(row["size_total_unc_plus_uas"])
        self.assertEqual(row["size_total_status"], "incomplete_unquantified_systematic")
        self.assertEqual(
            row["legacy_uncertainty_basis"],
            "symmetrized_fit_plus_absolute_position_component_total_incomplete",
        )
        self.assertIn("fitting plus absolute-position", row["uncertainty_status"])
        self.assertNotIn("numerical formal term only", row["uncertainty_status"])

        self.assertIsNone(row["distance_stat_unc_minus_mpc"])
        self.assertIsNone(row["distance_stat_unc_plus_mpc"])
        self.assertEqual(row["distance_stat_component_status"], "not_available")
        self.assertEqual(row["distance_systematic_unc_minus_mpc"], 2.09)
        self.assertEqual(row["distance_systematic_unc_plus_mpc"], 2.09)
        self.assertEqual(row["distance_unc_minus_mpc"], 2.09)
        self.assertEqual(row["distance_unc_plus_mpc"], 2.09)
        self.assertEqual(row["mass_shared_distance_unc_minus_msun"], 924_000.0)
        self.assertEqual(row["mass_shared_distance_unc_plus_msun"], 924_000.0)
        self.assertEqual(row["mass_total_unc_minus_msun"], 945_397.2709924649)
        self.assertEqual(row["mass_total_unc_plus_msun"], 945_397.2709924649)

        self.assertEqual(check["check_key"], "kuo_2020_ngc4388_turbulent_support_bias")
        self.assertEqual(
            check["data_relation"],
            "published_analytic_model_applied_to_ngc4388_disk_parameters",
        )
        self.assertEqual(check["lower_primary_ratio_factor"], 1.0)
        self.assertEqual(check["upper_primary_ratio_factor"], 1.001)
        self.assertEqual(check["assumed_inner_radius_pc"], 0.24)
        self.assertEqual(check["assumed_turbulent_velocity_kms"], 10.0)
        self.assertEqual(check["assumed_surface_density_power_law_index"], -1.0)
        self.assertEqual(check["assumed_rotation_profile"], "keplerian_rotation")
        self.assertEqual(check["mass_bias_eta"], 2.5)
        self.assertEqual(check["probability_status"], "non_probabilistic")
        self.assertEqual(check["adoption_status"], "reported_separately_not_adopted")
        self.assertEqual(check["reference_doi"], "10.1093/mnras/staa2260")
        self.assertEqual(row["size_model_check_lower_uas"], row["size_uas"])
        self.assertEqual(row["size_model_check_upper_uas"], row["size_uas"] * 1.001)
        self.assertEqual(row["size_model_check_upper_uas"], 0.04539588179014137)

        self.assertEqual(alternate["ratio_msun_per_mpc"], 423_300.0)
        self.assertEqual(alternate["ratio_unc_minus_msun_per_mpc"], 26_800.0)
        self.assertEqual(alternate["ratio_unc_plus_msun_per_mpc"], 18_500.0)
        self.assertEqual(alternate["size_uas"], 0.043421514018549316)
        self.assertEqual(alternate["size_unc_minus_uas"], 0.0027491060139313056)
        self.assertEqual(alternate["size_unc_plus_uas"], 0.0018977037782734758)
        self.assertEqual(row["nonadopted_fit_size_uas"], alternate["size_uas"])
        self.assertNotEqual(row["nonadopted_fit_size_uas"], row["size_uas"])

    def test_ngc4388_source_provenance_and_component_split_fail_closed(self) -> None:
        mutations = (
            ("reference_mass_msun", "8500000"),
            ("reference_mass_unc_minus_msun", "190000"),
            ("reference_distance_mpc", "19.1"),
            ("adopted_distance_mpc", "19.1"),
            ("distance_measurement_unc_minus_mpc", "2.09"),
            ("distance_measurement_unc_plus_mpc", "2.09"),
            ("distance_shared_systematic_unc_mpc", "2.08"),
            ("distance_stat_component_status", "statistical"),
            ("uncertainty_note", "fixed_distance_fit_error_only_sparse_disk_model_caveat"),
            ("mass_reference_doi", "10.0000/wrong"),
            ("mass_uncertainty_reference_doi", "10.0000/wrong"),
            ("distance_reference_doi", "10.0000/wrong"),
        )
        for field, value in mutations:
            with self.subTest(field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in self.inputs]
                next(item for item in rows if item["source_key"] == "NGC4388")[field] = value
                path = self._write_source_inputs(tmpdir, rows)
                with self.assertRaisesRegex(ValueError, "NGC4388"):
                    bhs.load_source_inputs(path)

        stale_outputs = (
            ("size_stat_component_status", "fixed_distance_formal_fit_interval"),
            (
                "legacy_uncertainty_basis",
                "symmetrized_formal_statistical_term_total_incomplete",
            ),
            ("distance_stat_component_status", "statistical"),
            ("distance_stat_unc_minus_mpc", 2.09),
            ("distance_systematic_unc_minus_mpc", None),
        )
        for field, value in stale_outputs:
            with self.subTest(output_field=field):
                rows = [dict(item) for item in self.rows]
                next(item for item in rows if item["source_key"] == "NGC4388")[field] = value
                with self.assertRaisesRegex(ValueError, "NGC4388"):
                    bhs.render_csv(rows)

    def test_circinus_components_display_interval_and_model_check_are_separate(self) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "Circinus")
        check = self.model_checks["Circinus"]
        row = self.by_key["Circinus"]
        expected_size = 0.04151990248474786
        radius_fraction = float(source["size_stat_fraction_minus_override"])
        expected_radius_component = expected_size * radius_fraction
        expected_mass_radius = 1_700_000.0 * radius_fraction
        expected_mass_distance = 1_700_000.0 * 0.8 / 4.2
        expected_mass_diagnostic = math.hypot(
            expected_mass_radius, expected_mass_distance
        )

        self.assertEqual(source["mass_display_interval_status"], "source_reported_total")
        self.assertEqual(
            source["distance_stat_component_status"],
            "reported_legacy_multimethod",
        )
        self.assertEqual(row["size_uas"], expected_size)
        self.assertEqual(row["size_stat_unc_minus_uas"], expected_radius_component)
        self.assertEqual(row["size_stat_unc_plus_uas"], expected_radius_component)
        self.assertEqual(row["size_unc_adopted_uas"], expected_radius_component)
        self.assertEqual(
            row["size_stat_component_status"], "angular_radius_only_component"
        )
        self.assertEqual(row["size_systematic_status"], "unquantified_model")
        self.assertEqual(
            row["size_total_status"], "incomplete_unquantified_systematic"
        )
        self.assertIsNone(row["size_total_unc_minus_uas"])
        self.assertIsNone(row["size_total_unc_plus_uas"])
        self.assertEqual(
            row["legacy_uncertainty_basis"],
            "symmetrized_angular_radius_only_component_total_incomplete",
        )
        self.assertNotIn("formal", row["size_stat_component_status"])
        self.assertNotIn("formal", row["legacy_uncertainty_basis"])
        self.assertNotIn("formal", row["uncertainty_status"])

        self.assertEqual(row["mass_fit_unc_minus_msun"], expected_mass_radius)
        self.assertEqual(row["mass_shared_distance_unc_minus_msun"], expected_mass_distance)
        self.assertEqual(row["mass_total_unc_minus_msun"], expected_mass_diagnostic)
        self.assertEqual(row["mass_total_unc_plus_msun"], expected_mass_diagnostic)
        self.assertEqual(row["mass_display_unc_minus_msun"], 300_000.0)
        self.assertEqual(row["mass_display_unc_plus_msun"], 300_000.0)
        self.assertEqual(row["mass_display_interval_status"], "source_reported_total")
        self.assertEqual(round(expected_mass_diagnostic, -5), 300_000.0)
        self.assertNotEqual(expected_mass_diagnostic, 300_000.0)
        self.assertEqual(
            row["distance_stat_component_status"],
            "reported_legacy_multimethod",
        )

        self.assertEqual(check["check_type"], "turbulent_support_mass_bias")
        self.assertEqual(
            check["data_relation"],
            "published_analytic_model_applied_to_circinus_disk_parameters",
        )
        self.assertEqual(check["range_basis"], "one_sided_fractional_mass_underestimate")
        self.assertEqual(check["endpoint_derivation"], "adopted_primary_ratio_factors")
        self.assertEqual(check["sensitivity_direction"], "upward_only")
        self.assertEqual(check["lower_primary_ratio_factor"], 1.0)
        self.assertEqual(check["upper_primary_ratio_factor"], 1.03)
        self.assertEqual(check["assumed_inner_radius_pc"], 0.11)
        self.assertEqual(check["assumed_turbulent_velocity_kms"], 30.0)
        self.assertEqual(check["assumed_surface_density_power_law_index"], -1.0)
        self.assertEqual(check["assumed_rotation_profile"], "keplerian_rotation")
        self.assertEqual(check["mass_bias_eta"], 2.5)
        self.assertEqual(
            check["assumption_profile"],
            "constant_turbulent_velocity_and_turbulent_pressure_dominates_thermal",
        )
        self.assertEqual(check["probability_status"], "non_probabilistic")
        self.assertEqual(check["adoption_status"], "reported_separately_not_adopted")
        self.assertEqual(check["reference_key"], "Kuo2020")
        self.assertEqual(check["reference_doi"], "10.1093/mnras/staa2260")
        self.assertEqual(row["size_model_check_lower_uas"], expected_size)
        self.assertEqual(row["size_model_check_upper_uas"], expected_size * 1.03)
        self.assertEqual(row["size_model_check_upper_uas"], 0.04276549955929029)
        self.assertNotEqual(
            row["size_model_check_upper_uas"], row["size_unc_adopted_uas"]
        )

    def test_ngc4945_model_check_is_one_sided_separate_and_reproducible(self) -> None:
        check = self.model_checks["NGC4945"]
        row = self.by_key["NGC4945"]
        expected_lower = 0.03881351774408862
        expected_upper = expected_lower * 1.24

        self.assertEqual(check["check_type"], "turbulent_support_mass_bias")
        self.assertEqual(
            check["data_relation"],
            "published_analytic_model_applied_to_ngc4945_disk_parameters",
        )
        self.assertEqual(
            check["range_basis"], "one_sided_fractional_mass_underestimate"
        )
        self.assertEqual(check["endpoint_derivation"], "adopted_primary_ratio_factors")
        self.assertEqual(check["sensitivity_direction"], "upward_only")
        self.assertIsNone(check["accepted_configuration_count"])
        self.assertIsNone(check["excluded_configuration_ids"])
        self.assertEqual(check["lower_primary_ratio_factor"], 1.0)
        self.assertEqual(check["upper_primary_ratio_factor"], 1.24)
        self.assertEqual(check["assumed_inner_radius_pc"], 0.16)
        self.assertEqual(check["assumed_turbulent_velocity_kms"], 60.0)
        self.assertEqual(check["assumed_surface_density_power_law_index"], -1.0)
        self.assertEqual(check["assumed_rotation_profile"], "keplerian_rotation")
        self.assertEqual(check["mass_bias_eta"], 2.5)
        self.assertEqual(
            check["assumption_profile"],
            "constant_turbulent_velocity_and_turbulent_pressure_dominates_thermal",
        )
        self.assertEqual(check["probability_status"], "non_probabilistic")
        self.assertEqual(check["adoption_status"], "reported_separately_not_adopted")
        self.assertEqual(check["reference_key"], "Kuo2020")
        self.assertEqual(check["reference_doi"], "10.1093/mnras/staa2260")

        self.assertEqual(row["size_model_check_lower_uas"], expected_lower)
        self.assertEqual(row["size_model_check_upper_uas"], expected_upper)
        self.assertEqual(expected_upper, 0.04812876200266989)
        self.assertEqual(expected_upper - expected_lower, 0.009315244258581268)

        # The model check is not an adopted error, systematic, or quadrature term.
        self.assertEqual(row["size_uas"], expected_lower)
        self.assertEqual(row["mass_msun"], 1_400_000.0)
        self.assertEqual(row["distance_mpc"], 3.7)
        self.assertEqual(row["log_mass"], 6.146128035678238)
        self.assertEqual(row["parameter_reference_doi"], "10.1086/310643")
        self.assertIsNone(row["mass_uncertainty_reference_doi"])
        self.assertIsNone(row["size_unc_adopted_uas"])
        self.assertIsNone(row["size_stat_unc_minus_uas"])
        self.assertIsNone(row["size_stat_unc_plus_uas"])
        self.assertIsNone(row["size_systematic_unc_minus_uas"])
        self.assertIsNone(row["size_systematic_unc_plus_uas"])
        self.assertIsNone(row["size_total_unc_minus_uas"])
        self.assertIsNone(row["size_total_unc_plus_uas"])
        self.assertEqual(row["size_stat_component_status"], "no_formal_uncertainty")
        self.assertEqual(row["size_systematic_status"], "no_formal_uncertainty")
        self.assertEqual(row["size_total_status"], "no_formal_uncertainty")
        self.assertEqual(row["legacy_uncertainty_basis"], "not_available")
        self.assertIn("one-sided model sensitivity not combined", row["uncertainty_status"])

    def test_model_check_loader_fails_closed(self) -> None:
        with bhs.DEFAULT_MODEL_CHECKS.open(newline="", encoding="utf-8") as handle:
            originals = list(csv.DictReader(handle))
        by_source = {row["source_key"]: row for row in originals}
        original = by_source["NGC4258"]

        mutations = (
            ("check_key", "different_check", "check_key does not match"),
            ("check_type", "unknown", "unsupported model-check type"),
            ("data_relation", "unknown", "unsupported model-check data relation"),
            ("range_basis", "unknown", "unsupported model-check range basis"),
            ("endpoint_derivation", "unknown", "endpoint derivation"),
            ("sensitivity_direction", "unknown", "sensitivity direction"),
            ("probability_status", "probabilistic", "probability status"),
            ("adoption_status", "combine_in_quadrature", "adoption status"),
            ("reference_key", "unknown", "unsupported model-check reference_key"),
            ("reference_doi", "10.0000/wrong", "does not match its DOI"),
            ("accepted_configuration_count", "0", "must be positive"),
            ("accepted_configuration_count", "14", "does not match"),
            ("excluded_configuration_ids", "2;two", "malformed"),
            ("excluded_configuration_ids", "2;16", "does not match"),
            ("lower_configuration_id", "3", "does not match"),
            ("upper_configuration_id", "16", "does not match"),
            ("lower_mass_msun", "-1", "must be positive"),
            ("lower_mass_msun", "42000000", "does not match"),
            ("upper_mass_msun", "1", "must be ordered"),
            (
                "lower_primary_ratio_factor",
                "1",
                "must be literal NA",
            ),
            ("assumed_rotation_profile", "flat_rotation", "must be literal NA"),
            ("mass_bias_eta", "2", "must be literal NA"),
        )
        for field, value, message in mutations:
            with self.subTest(field=field), tempfile.TemporaryDirectory() as tmpdir:
                mutated = dict(original)
                mutated[field] = value
                path = self._write_model_checks(tmpdir, [mutated])
                with self.assertRaisesRegex(ValueError, message):
                    bhs.load_model_checks(path, known_sources=self.source_keys)

        ngc3079 = by_source["NGC3079"]
        ngc3079_mutations = (
            ("check_key", "different_check", "check_key does not match"),
            ("check_type", "unknown", "unsupported model-check type"),
            ("data_relation", "unknown", "unsupported model-check data relation"),
            ("range_basis", "unknown", "unsupported model-check range basis"),
            ("endpoint_derivation", "unknown", "endpoint derivation"),
            ("sensitivity_direction", "two_sided_range", "does not match"),
            ("accepted_configuration_count", "1", "must be literal NA"),
            ("lower_mass_msun", "2000000", "must be literal NA"),
            ("lower_primary_ratio_factor", "0", "positive and ordered"),
            ("lower_primary_ratio_factor", "0.999", "does not match"),
            ("upper_primary_ratio_factor", "1", "positive and ordered"),
            ("upper_primary_ratio_factor", "1.35", "does not match"),
            ("assumed_inner_radius_pc", "0", "must be positive"),
            ("assumed_inner_radius_pc", "0.41", "does not match"),
            ("assumed_turbulent_velocity_kms", "0", "must be positive"),
            ("assumed_turbulent_velocity_kms", "61", "does not match"),
            (
                "assumed_surface_density_power_law_index",
                "NA",
                "power-law index is required",
            ),
            (
                "assumed_surface_density_power_law_index",
                "-0.9",
                "does not match",
            ),
            ("assumed_rotation_profile", "NA", "rotation_profile is required"),
            (
                "assumed_rotation_profile",
                "unknown_rotation",
                "unsupported model-check rotation profile",
            ),
            (
                "assumed_rotation_profile",
                "keplerian_rotation",
                "does not match",
            ),
            ("mass_bias_eta", "NA", "mass_bias_eta is required"),
            ("mass_bias_eta", "0", "mass_bias_eta must be positive"),
            ("mass_bias_eta", "2.5", "does not match"),
            ("assumption_profile", "unknown", "assumption profile"),
            ("probability_status", "probabilistic", "probability status"),
            ("adoption_status", "combine_in_quadrature", "adoption status"),
            ("reference_key", "vanDerBoom2026", "does not match"),
            ("reference_doi", "10.0000/wrong", "does not match its DOI"),
            ("provenance_note", "different", "does not match"),
        )
        for field, value, message in ngc3079_mutations:
            with self.subTest(source="NGC3079", field=field), tempfile.TemporaryDirectory() as tmpdir:
                mutated = dict(ngc3079)
                mutated[field] = value
                path = self._write_model_checks(tmpdir, [mutated])
                with self.assertRaisesRegex(ValueError, message):
                    bhs.load_model_checks(path, known_sources=self.source_keys)

        ngc4945 = by_source["NGC4945"]
        ngc4945_mutations = (
            ("check_key", "different_check", "check_key does not match"),
            ("check_type", "unknown", "unsupported model-check type"),
            (
                "data_relation",
                "published_analytic_model_applied_to_ngc3079_disk_parameters",
                "does not match",
            ),
            ("range_basis", "unknown", "unsupported model-check range basis"),
            ("endpoint_derivation", "unknown", "endpoint derivation"),
            ("sensitivity_direction", "two_sided_range", "does not match"),
            ("accepted_configuration_count", "1", "must be literal NA"),
            ("lower_mass_msun", "1400000", "must be literal NA"),
            ("lower_primary_ratio_factor", "0", "positive and ordered"),
            ("lower_primary_ratio_factor", "1.01", "does not match"),
            ("upper_primary_ratio_factor", "1", "positive and ordered"),
            ("upper_primary_ratio_factor", "1.34", "does not match"),
            ("assumed_inner_radius_pc", "0", "must be positive"),
            ("assumed_inner_radius_pc", "0.40", "does not match"),
            ("assumed_turbulent_velocity_kms", "0", "must be positive"),
            ("assumed_turbulent_velocity_kms", "61", "does not match"),
            (
                "assumed_surface_density_power_law_index",
                "NA",
                "power-law index is required",
            ),
            (
                "assumed_surface_density_power_law_index",
                "-0.9",
                "does not match",
            ),
            ("assumed_rotation_profile", "NA", "rotation_profile is required"),
            (
                "assumed_rotation_profile",
                "unknown_rotation",
                "unsupported model-check rotation profile",
            ),
            ("assumed_rotation_profile", "flat_rotation", "does not match"),
            ("mass_bias_eta", "NA", "mass_bias_eta is required"),
            ("mass_bias_eta", "0", "mass_bias_eta must be positive"),
            ("mass_bias_eta", "2", "does not match"),
            (
                "assumption_profile",
                "constant_turbulence_and_wind_pressure",
                "assumption profile",
            ),
            ("probability_status", "probabilistic", "probability status"),
            ("adoption_status", "combine_in_quadrature", "adoption status"),
            ("reference_key", "vanDerBoom2026", "does not match"),
            ("reference_doi", "10.0000/wrong", "does not match its DOI"),
            ("provenance_note", "different", "does not match"),
        )
        for field, value, message in ngc4945_mutations:
            with self.subTest(source="NGC4945", field=field), tempfile.TemporaryDirectory() as tmpdir:
                mutated = dict(ngc4945)
                mutated[field] = value
                path = self._write_model_checks(tmpdir, [mutated])
                with self.assertRaisesRegex(ValueError, message):
                    bhs.load_model_checks(path, known_sources=self.source_keys)

        ngc4388 = by_source["NGC4388"]
        ngc4388_mutations = (
            ("check_key", "different_check", "check_key does not match"),
            (
                "data_relation",
                "published_analytic_model_applied_to_ngc4945_disk_parameters",
                "does not match",
            ),
            ("sensitivity_direction", "two_sided_range", "does not match"),
            ("accepted_configuration_count", "1", "must be literal NA"),
            ("lower_mass_msun", "8400000", "must be literal NA"),
            ("lower_primary_ratio_factor", "0.999", "does not match"),
            ("upper_primary_ratio_factor", "1.002", "does not match"),
            ("assumed_inner_radius_pc", "0.25", "does not match"),
            ("assumed_turbulent_velocity_kms", "11", "does not match"),
            ("assumed_surface_density_power_law_index", "-0.9", "does not match"),
            ("assumed_rotation_profile", "flat_rotation", "does not match"),
            ("mass_bias_eta", "2", "does not match"),
            ("assumption_profile", "unknown", "assumption profile"),
            ("probability_status", "probabilistic", "probability status"),
            ("adoption_status", "combine_in_quadrature", "adoption status"),
            ("reference_key", "vanDerBoom2026", "does not match"),
            ("reference_doi", "10.0000/wrong", "does not match its DOI"),
            ("provenance_note", "different", "does not match"),
        )
        for field, value, message in ngc4388_mutations:
            with self.subTest(source="NGC4388", field=field), tempfile.TemporaryDirectory() as tmpdir:
                mutated = dict(ngc4388)
                mutated[field] = value
                path = self._write_model_checks(tmpdir, [mutated])
                with self.assertRaisesRegex(ValueError, message):
                    bhs.load_model_checks(path, known_sources=self.source_keys)

        circinus = by_source["Circinus"]
        circinus_mutations = (
            ("check_key", "different_check", "check_key does not match"),
            (
                "data_relation",
                "published_analytic_model_applied_to_ngc4945_disk_parameters",
                "does not match",
            ),
            ("sensitivity_direction", "two_sided_range", "does not match"),
            ("accepted_configuration_count", "1", "must be literal NA"),
            ("lower_mass_msun", "1700000", "must be literal NA"),
            ("lower_primary_ratio_factor", "1.01", "does not match"),
            ("upper_primary_ratio_factor", "1.04", "does not match"),
            ("assumed_inner_radius_pc", "0.12", "does not match"),
            ("assumed_turbulent_velocity_kms", "31", "does not match"),
            ("assumed_surface_density_power_law_index", "-0.9", "does not match"),
            ("assumed_rotation_profile", "flat_rotation", "does not match"),
            ("mass_bias_eta", "2", "does not match"),
            ("assumption_profile", "unknown", "assumption profile"),
            ("probability_status", "probabilistic", "probability status"),
            ("adoption_status", "combine_in_quadrature", "adoption status"),
            ("reference_doi", "10.0000/wrong", "does not match its DOI"),
            ("provenance_note", "different", "does not match"),
        )
        for field, value, message in circinus_mutations:
            with self.subTest(source="Circinus", field=field), tempfile.TemporaryDirectory() as tmpdir:
                mutated = dict(circinus)
                mutated[field] = value
                path = self._write_model_checks(tmpdir, [mutated])
                with self.assertRaisesRegex(ValueError, message):
                    bhs.load_model_checks(path, known_sources=self.source_keys)

        with tempfile.TemporaryDirectory() as tmpdir:
            unknown = dict(original)
            unknown["source_key"] = "UNKNOWN"
            path = self._write_model_checks(tmpdir, [unknown])
            with self.assertRaisesRegex(ValueError, "unknown model-sensitivity source"):
                bhs.load_model_checks(path, known_sources=self.source_keys)

        with tempfile.TemporaryDirectory() as tmpdir:
            duplicate = dict(original)
            path = self._write_model_checks(tmpdir, [original, duplicate])
            with self.assertRaisesRegex(ValueError, "multiple model-sensitivity checks"):
                bhs.load_model_checks(path, known_sources=self.source_keys)

        with tempfile.TemporaryDirectory() as tmpdir:
            wrong_reviewed_source = dict(original)
            wrong_reviewed_source["source_key"] = "NGC1194"
            path = self._write_model_checks(tmpdir, [wrong_reviewed_source])
            with self.assertRaisesRegex(ValueError, "unsupported reviewed"):
                bhs.load_model_checks(path, known_sources=self.source_keys)

        with tempfile.TemporaryDirectory() as tmpdir:
            fields = bhs.MODEL_CHECK_FIELDS[:-1]
            path = self._write_model_checks(
                tmpdir,
                [{field: original[field] for field in fields}],
                fields=fields,
            )
            with self.assertRaisesRegex(ValueError, "unexpected columns"):
                bhs.load_model_checks(path, known_sources=self.source_keys)

    def test_model_check_endpoint_changes_cannot_change_the_adopted_result(self) -> None:
        for source_key in ("Circinus", "NGC3079", "NGC4388", "NGC4945"):
            with self.subTest(source=source_key):
                source = next(
                    row for row in self.inputs if row["source_key"] == source_key
                )
                baseline = bhs.build_shared_row(source)
                mutated = dict(self.model_checks[source_key])
                mutated["upper_primary_ratio_factor"] = 9.9
                candidate = bhs.attach_model_check(dict(baseline), mutated)

                self.assertNotEqual(
                    candidate["size_model_check_upper_uas"],
                    self.by_key[source_key]["size_model_check_upper_uas"],
                )
                for field in bhs.OUTPUT_FIELDS:
                    if (
                        field.startswith("size_model_check_")
                        or field == "uncertainty_status"
                    ):
                        continue
                    self.assertEqual(candidate[field], baseline[field], field)

    def test_model_check_attachment_fails_closed(self) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "NGC3079")
        baseline = bhs.build_shared_row(source)
        check = dict(self.model_checks["NGC3079"])

        mismatched = dict(check)
        mismatched["source_key"] = "NGC4258"
        with self.assertRaisesRegex(ValueError, "does not match output source"):
            bhs.attach_model_check(dict(baseline), mismatched)

        missing_primary = dict(baseline)
        missing_primary["size_uas"] = None
        with self.assertRaisesRegex(
            ValueError, "requires the retained primary mass and size"
        ):
            bhs.attach_model_check(missing_primary, check)

        wrong_policy = dict(baseline)
        wrong_policy["size_total_status"] = "formal_component_only"
        with self.assertRaisesRegex(ValueError, "primary uncertainty policy mismatch"):
            bhs.attach_model_check(wrong_policy, check)

        already_attached = bhs.attach_model_check(dict(baseline), check)
        with self.assertRaisesRegex(ValueError, "already attached"):
            bhs.attach_model_check(already_attached, check)

        ngc4945_source = next(
            row for row in self.inputs if row["source_key"] == "NGC4945"
        )
        ngc4945_baseline = bhs.build_shared_row(ngc4945_source)
        with self.assertRaisesRegex(ValueError, "does not match output source"):
            bhs.attach_model_check(ngc4945_baseline, check)

    def test_nonadopted_ratio_fits_are_separate_and_reproducible(self) -> None:
        self.assertEqual(len(self.nonadopted_fits), 13)
        self.assertEqual(
            set(self.nonadopted_fits),
            {
                "NGC1194",
                "NGC2273",
                "NGC4388",
                "CGCG074-064",
                "NGC6264",
                "UGC3789",
                "ESO558-G009",
                "NGC2960",
                "NGC5765b",
                "NGC6323",
                "J0437",
                "Mrk1029",
                "NGC1320",
            },
        )
        expected = {
            "NGC1194": {
                "values": (1_310_000.0, 20_900.0, 21_100.0),
                "treatment": "shared_distance",
                "primary_result_status": "retained",
                "data_relation": "reanalysis_of_published_kuo2011_vlbi_spot_data",
                "reference_key": "GonzalezJuarez2024",
                "reference_doi": "10.1051/0004-6361/202450098",
                "basis": (
                    "same_data_reanalysis_conditional_on_restricted_disk_model_"
                    "and_discrepant_recession_velocity"
                ),
            },
            "NGC2273": {
                "values": (328_300.0, 5_100.0, 4_500.0),
                "treatment": "shared_distance",
                "primary_result_status": "retained",
                "data_relation": "reanalysis_of_published_kuo2011_vlbi_spot_data",
                "reference_key": "Villaraos2022",
                "reference_doi": "10.1093/mnras/stac2973",
                "basis": (
                    "same_data_reanalysis_conditional_on_fixed_vertical_center_"
                    "and_selected_scatter_angle"
                ),
            },
            "NGC4388": {
                "values": (423_300.0, 26_800.0, 18_500.0),
                "treatment": "shared_distance",
                "primary_result_status": "retained",
                "data_relation": "reanalysis_of_published_kuo2011_vlbi_spot_data",
                "reference_key": "Villaraos2022",
                "reference_doi": "10.1093/mnras/stac2973",
                "basis": (
                    "same_data_reanalysis_conditional_on_fixed_vertical_center_"
                    "and_selected_scatter_angle"
                ),
            },
            "CGCG074-064": {
                "values": (274_900.0, 1_500.0, 1_700.0),
                "treatment": "paired_posterior",
                "primary_result_status": "retained",
                "data_relation": (
                    "reanalysis_of_published_pesce2020a_vlbi_spot_data"
                ),
                "reference_key": "Villaraos2022",
                "reference_doi": "10.1093/mnras/stac2973",
                "basis": (
                    "same_data_reanalysis_conditional_on_fixed_pesce2020a_"
                    "inclination_and_selected_scatter_angle"
                ),
                "provenance": (
                    "published_table_2_ratio_and_formal_interval_table_2_scatter_"
                    "11deg_chi2red_1p289_conflicts_with_section_4p1p4_scatter_"
                    "10deg_chi2red_1p600_no_scatter_or_chi2_value_adopted"
                ),
            },
            "NGC6264": {
                "values": (212_600.0, 700.0, 700.0),
                "treatment": "paired_posterior",
                "primary_result_status": "retained",
                "data_relation": "reanalysis_of_published_kuo2011_vlbi_spot_data",
                "reference_key": "Villaraos2022",
                "reference_doi": "10.1093/mnras/stac2973",
                "basis": (
                    "older_partial_data_reanalysis_conditional_on_fixed_vertical_center_"
                    "and_selected_scatter_angle"
                ),
                "provenance": (
                    "published_table_2_ratio_and_formal_interval_table_2_scatter_"
                    "10deg_chi2red_1p304_conflicts_with_section_4p1p7_chi2red_"
                    "1p301"
                ),
            },
            "UGC3789": {
                "values": (227_700.0, 900.0, 900.0),
                "treatment": "paired_posterior",
                "primary_result_status": "retained",
                "data_relation": "reanalysis_of_published_reid2013_vlbi_spot_data",
                "reference_key": "Villaraos2022",
                "reference_doi": "10.1093/mnras/stac2973",
                "basis": (
                    "same_observation_vlbi_spot_subset_reanalysis_conditional_on_"
                    "fixed_reid2013_inclination_and_selected_scatter_angle"
                ),
                "provenance": (
                    "published_table_2_ratio_and_formal_interval_table_2_and_"
                    "section_4p1p3_scatter_10deg_chi2red_1p430"
                ),
            },
            "ESO558-G009": {
                "values": (159_400.0, 2_600.0, 2_500.0),
                "treatment": "shared_distance",
                "primary_result_status": "retained",
                "data_relation": (
                    "reanalysis_of_published_gao2017_vlbi_spot_data"
                ),
                "reference_key": "Villaraos2022",
                "reference_doi": "10.1093/mnras/stac2973",
                "basis": (
                    "same_data_reanalysis_conditional_on_fixed_gao2017_"
                    "inclination_and_selected_scatter_angle"
                ),
                "provenance": (
                    "published_table_2_ratio_and_formal_interval_table_2_and_"
                    "section_4p1p5_fixed_inclination_98deg_scatter_7deg_"
                    "chi2red_1p188"
                ),
            },
            "NGC2960": {
                "values": (173_800.0, 1_600.0, 1_500.0),
                "treatment": "shared_distance",
                "primary_result_status": "retained",
                "data_relation": "reanalysis_of_published_kuo2011_vlbi_spot_data",
                "reference_key": "Villaraos2022",
                "reference_doi": "10.1093/mnras/stac2973",
                "basis": (
                    "older_published_data_reanalysis_conditional_on_fixed_kuo2011_"
                    "inclination_vertical_center_and_selected_scatter_angle"
                ),
                "provenance": (
                    "published_table_2_ratio_and_formal_interval_fixed_i89_y0_"
                    "minus0p020_scatter15deg_table_chi2red1p488_conflicts_with_"
                    "section_4p1p6_y0_origin_x0_0p004_chi2red1p511_section_4p1p6_"
                    "y0_variation_up_to_4_percent_not_in_formal_interval_exact_spot_"
                    "overlap_with_kuo2018_unverifiable"
                ),
            },
            "NGC5765b": {
                "values": (372_700.0, 1_300.0, 1_300.0),
                "treatment": "paired_posterior",
                "primary_result_status": "retained",
                "data_relation": "reanalysis_of_published_gao2016_vlbi_spot_data",
                "reference_key": "Villaraos2022",
                "reference_doi": "10.1093/mnras/stac2973",
                "basis": (
                    "same_data_reanalysis_conditional_on_fixed_gao2016_inclination_"
                    "and_selected_scatter_angle"
                ),
            },
            "NGC6323": {
                "values": (91_600.0, 500.0, 500.0),
                "treatment": "paired_posterior",
                "primary_result_status": "retained",
                "data_relation": "reanalysis_of_published_kuo2011_vlbi_spot_data",
                "reference_key": "Villaraos2022",
                "reference_doi": "10.1093/mnras/stac2973",
                "basis": (
                    "older_partial_data_reanalysis_conditional_on_fixed_kuo2015_"
                    "inclination_and_selected_scatter_angle"
                ),
                "provenance": (
                    "published_table_2_ratio_and_formal_interval_four_parameters_"
                    "y0_fitted_table_2_scatter_11deg_chi2red_1p276_conflicts_with_"
                    "section_4p1p2_chi2red_1p251"
                ),
            },
            "J0437": {
                "values": (42_500.0, 1_200.0, 1_000.0),
                "treatment": "shared_distance",
                "primary_result_status": "retained",
                "data_relation": "reanalysis_of_published_gao2017_vlbi_spot_data",
                "reference_key": "Villaraos2022",
                "reference_doi": "10.1093/mnras/stac2973",
                "basis": (
                    "same_data_reanalysis_conditional_on_fixed_gao2017_inclination_"
                    "vertical_center_and_zero_scatter_angle"
                ),
                "provenance": (
                    "published_table_2_ratio_and_formal_interval_three_parameters_"
                    "fixed_i81_y0_minus0p040_zero_scatter_chi2red_1p038_section_"
                    "4p1p8_y0_variation_up_to_2_percent_within_formal_interval"
                ),
            },
            "Mrk1029": {
                "values": (14_400.0, 1_100.0, 1_100.0),
                "treatment": "no_robust_mass",
                "primary_result_status": "absent",
                "data_relation": "reanalysis_of_published_gao2017_vlbi_spot_data",
                "reference_key": "GonzalezJuarez2024",
                "reference_doi": "10.1051/0004-6361/202450098",
                "basis": (
                    "restricted_same_data_fit_conditional_on_fixed_gao2017_"
                    "inclination_vertical_center_and_selected_scatter_angle_"
                    "does_not_resolve_primary_complex_geometry_and_dynamics"
                ),
                "provenance": (
                    "published_table_2_ratio_and_formal_interval_three_parameters_"
                    "m_over_d_zrec_x0_fitted_y0_fixed_i79_fixed_scatter14deg_"
                    "chi2red1p385_pa_plusminus10deg_sensitivity_1p19percent_"
                    "neglected_within_formal_interval"
                ),
            },
            "NGC1320": {
                "values": (149_700.0, 6_900.0, 7_100.0),
                "treatment": "no_robust_mass",
                "primary_result_status": "absent",
                "data_relation": "reanalysis_of_published_gao2017_vlbi_spot_data",
                "reference_key": "GonzalezJuarez2024",
                "reference_doi": "10.1051/0004-6361/202450098",
                "basis": (
                    "restricted_same_data_fit_conditional_on_fixed_edge_on_"
                    "inclination_and_selected_scatter_angle_does_not_resolve_"
                    "primary_complex_geometry_and_dynamics"
                ),
                "provenance": (
                    "published_table_2_ratio_and_formal_interval_four_parameters_"
                    "m_over_d_zrec_x0_y0_fitted_reference_midpoint_between_"
                    "redshifted_and_blueshifted_group_geometric_centers_i90_fixed_"
                    "scatter17deg_chi2red1p127_pa_plusminus10deg_sensitivity_"
                    "2p64percent_within_formal_interval_bh_position_near_"
                    "blueshifted_group_new_vlbi_requested"
                ),
            },
        }
        for key, expectation in expected.items():
            ratio, minus, plus = expectation["values"]
            fit = self.nonadopted_fits[key]
            row = self.by_key[key]
            self.assertEqual(fit["ratio_msun_per_mpc"], ratio)
            self.assertEqual(fit["ratio_unc_minus_msun_per_mpc"], minus)
            self.assertEqual(fit["ratio_unc_plus_msun_per_mpc"], plus)
            self.assertEqual(fit["size_uas"], bhs.shadow_size_uas(ratio, 1.0))
            self.assertEqual(
                fit["size_unc_minus_uas"], bhs.shadow_size_uas(minus, 1.0)
            )
            self.assertEqual(
                fit["size_unc_plus_uas"], bhs.shadow_size_uas(plus, 1.0)
            )
            self.assertEqual(
                fit["interval_status"], "bayesian_formal_conditional_on_model"
            )
            self.assertEqual(
                fit["adoption_status"], "reported_separately_not_adopted"
            )
            self.assertEqual(fit["reference_key"], expectation["reference_key"])
            self.assertEqual(fit["reference_doi"], expectation["reference_doi"])
            self.assertEqual(
                fit["provenance_note"],
                expectation.get(
                    "provenance", "published_table_2_ratio_and_formal_interval"
                ),
            )
            self.assertEqual(fit["primary_treatment"], expectation["treatment"])
            self.assertEqual(
                fit["primary_result_status"], expectation["primary_result_status"]
            )
            self.assertEqual(fit["data_relation"], expectation["data_relation"])

            self.assertEqual(row["nonadopted_fit_ratio_msun_per_mpc"], ratio)
            self.assertEqual(
                row["nonadopted_fit_ratio_unc_minus_msun_per_mpc"], minus
            )
            self.assertEqual(
                row["nonadopted_fit_ratio_unc_plus_msun_per_mpc"], plus
            )
            self.assertEqual(row["nonadopted_fit_size_uas"], fit["size_uas"])
            self.assertEqual(
                row["nonadopted_fit_adoption_status"],
                "reported_separately_not_adopted",
            )
            self.assertEqual(
                row["nonadopted_fit_basis"], expectation["basis"]
            )

            if expectation["primary_result_status"] == "absent":
                # These conditional fits must not create an adopted mass or size.
                self.assertIsNone(row["mass_msun"])
                self.assertIsNone(row["size_uas"])
                self.assertIsNone(row["size_stat_unc_minus_uas"])
                self.assertIsNone(row["size_systematic_unc_minus_uas"])
                self.assertIsNone(row["size_total_unc_minus_uas"])
                self.assertEqual(row["adopted_uncertainty_model"], "not_available")
                self.assertEqual(row["ratio_provenance"], "no_robust_primary_mass")
                self.assertIn(
                    "reported separately and not adopted", row["uncertainty_status"]
                )

        for key, row in self.by_key.items():
            if key in expected:
                continue
            self.assertIsNone(row["nonadopted_fit_key"], key)
            self.assertIsNone(row["nonadopted_fit_size_uas"], key)

    def test_nonadopted_ratio_fit_loader_fails_closed(self) -> None:
        with bhs.DEFAULT_NONADOPTED_FITS.open(newline="", encoding="utf-8") as handle:
            originals = list(csv.DictReader(handle))

        mutations = (
            ("fit_key", "different_fit", "fit_key does not match"),
            ("data_relation", "unknown", "data relation"),
            ("model_basis", "unknown", "model basis"),
            ("interval_status", "unconditional", "interval status"),
            ("adoption_status", "adopted", "adoption status"),
            ("nonadoption_basis", "none", "non-adoption basis"),
            ("ratio_msun_per_mpc", "0", "must be positive"),
            (
                "ratio_unc_minus_msun_per_mpc",
                "2000000",
                "must be smaller than the ratio",
            ),
            ("ratio_unc_plus_msun_per_mpc", "1000", "does not match"),
            ("reference_key", "unknown", "reference_key"),
            ("reference_doi", "10.0000/wrong", "does not match its DOI"),
            ("provenance_note", "published_table_1", "does not match"),
            ("provenance_note", "", "provenance_note is required"),
        )
        for field, value, message in mutations:
            with self.subTest(field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(row) for row in originals]
                rows[0][field] = value
                path = self._write_nonadopted_fits(tmpdir, rows)
                with self.assertRaisesRegex(ValueError, message):
                    bhs.load_nonadopted_ratio_fits(
                        path, source_treatments=self.source_treatments
                    )

        with tempfile.TemporaryDirectory() as tmpdir:
            path = self._write_nonadopted_fits(tmpdir, originals[:1])
            with self.assertRaisesRegex(ValueError, "exactly match"):
                bhs.load_nonadopted_ratio_fits(
                    path, source_treatments=self.source_treatments
                )

        with tempfile.TemporaryDirectory() as tmpdir:
            rows = [dict(row) for row in originals]
            rows[0]["source_key"] = "UNKNOWN"
            path = self._write_nonadopted_fits(tmpdir, rows)
            with self.assertRaisesRegex(ValueError, "unknown.*source"):
                bhs.load_nonadopted_ratio_fits(
                    path, source_treatments=self.source_treatments
                )

        with tempfile.TemporaryDirectory() as tmpdir:
            rows = [dict(row) for row in originals]
            by_source = {row["source_key"]: row for row in rows}
            by_source["NGC1320"]["source_key"] = "Mrk1029"
            path = self._write_nonadopted_fits(tmpdir, rows)
            with self.assertRaisesRegex(ValueError, "duplicate.*source"):
                bhs.load_nonadopted_ratio_fits(
                    path, source_treatments=self.source_treatments
                )

        with tempfile.TemporaryDirectory() as tmpdir:
            rows = [dict(row) for row in originals]
            rows[1]["fit_key"] = rows[0]["fit_key"]
            path = self._write_nonadopted_fits(tmpdir, rows)
            with self.assertRaisesRegex(ValueError, "duplicate.*fit_key"):
                bhs.load_nonadopted_ratio_fits(
                    path, source_treatments=self.source_treatments
                )

        with tempfile.TemporaryDirectory() as tmpdir:
            fields = bhs.NONADOPTED_FIT_FIELDS[:-1]
            rows = [{field: row[field] for field in fields} for row in originals]
            path = self._write_nonadopted_fits(tmpdir, rows, fields=fields)
            with self.assertRaisesRegex(ValueError, "unexpected columns"):
                bhs.load_nonadopted_ratio_fits(
                    path, source_treatments=self.source_treatments
                )

        wrong_treatments = dict(self.source_treatments)
        wrong_treatments["Mrk1029"] = "shared_distance"
        with self.assertRaisesRegex(ValueError, "requires no_robust_mass"):
            bhs.load_nonadopted_ratio_fits(source_treatments=wrong_treatments)

        wrong_treatments = dict(self.source_treatments)
        wrong_treatments["NGC1194"] = "no_robust_mass"
        with self.assertRaisesRegex(ValueError, "requires shared_distance"):
            bhs.load_nonadopted_ratio_fits(source_treatments=wrong_treatments)

        wrong_treatments = dict(self.source_treatments)
        wrong_treatments["NGC5765b"] = "shared_distance"
        with self.assertRaisesRegex(ValueError, "requires paired_posterior"):
            bhs.load_nonadopted_ratio_fits(source_treatments=wrong_treatments)

        wrong_treatments = dict(self.source_treatments)
        wrong_treatments["NGC2273"] = "no_robust_mass"
        with self.assertRaisesRegex(ValueError, "requires shared_distance"):
            bhs.load_nonadopted_ratio_fits(source_treatments=wrong_treatments)

        wrong_treatments = dict(self.source_treatments)
        wrong_treatments["NGC4388"] = "no_robust_mass"
        with self.assertRaisesRegex(ValueError, "requires shared_distance"):
            bhs.load_nonadopted_ratio_fits(source_treatments=wrong_treatments)

        cross_source_mutations = (
            ("Mrk1029", "NGC1194", "data_relation"),
            ("Mrk1029", "NGC1194", "nonadoption_basis"),
            ("NGC2273", "NGC1194", "nonadoption_basis"),
        )
        for source_key, donor_key, field in cross_source_mutations:
            with self.subTest(
                source_key=source_key, cross_source_field=field
            ), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(row) for row in originals]
                by_source = {row["source_key"]: row for row in rows}
                by_source[source_key][field] = by_source[donor_key][field]
                path = self._write_nonadopted_fits(tmpdir, rows)
                with self.assertRaisesRegex(ValueError, "does not match the reviewed"):
                    bhs.load_nonadopted_ratio_fits(
                        path, source_treatments=self.source_treatments
                    )

        with tempfile.TemporaryDirectory() as tmpdir:
            rows = [dict(row) for row in originals]
            by_source = {row["source_key"]: row for row in rows}
            by_source["NGC2273"]["reference_key"] = "GonzalezJuarez2024"
            by_source["NGC2273"]["reference_doi"] = (
                "10.1051/0004-6361/202450098"
            )
            path = self._write_nonadopted_fits(tmpdir, rows)
            with self.assertRaisesRegex(ValueError, "does not match the reviewed"):
                bhs.load_nonadopted_ratio_fits(
                    path, source_treatments=self.source_treatments
                )

        ngc4388_mutations = (
            ("fit_key", "villaraos_2022_ngc4388_other_fit"),
            ("data_relation", "reanalysis_of_published_gao2017_vlbi_spot_data"),
            ("model_basis", "unknown_model"),
            ("ratio_msun_per_mpc", "423301"),
            ("ratio_unc_minus_msun_per_mpc", "26799"),
            ("ratio_unc_plus_msun_per_mpc", "18501"),
            ("interval_status", "unconditional"),
            ("adoption_status", "adopted"),
            (
                "nonadoption_basis",
                "does_not_resolve_primary_complex_geometry_and_dynamics",
            ),
            ("reference_key", "GonzalezJuarez2024"),
            ("reference_doi", "10.1051/0004-6361/202450098"),
            ("provenance_note", "published_table_1_ratio"),
        )
        for field, value in ngc4388_mutations:
            with self.subTest(source="NGC4388", field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(row) for row in originals]
                next(row for row in rows if row["source_key"] == "NGC4388")[field] = value
                path = self._write_nonadopted_fits(tmpdir, rows)
                with self.assertRaisesRegex(ValueError, "NGC4388"):
                    bhs.load_nonadopted_ratio_fits(
                        path, source_treatments=self.source_treatments
                    )

    def test_nonadopted_fit_attachment_cannot_change_an_adopted_field(self) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "Mrk1029")
        baseline = bhs.build_no_mass_row(source)
        fit = dict(self.nonadopted_fits["Mrk1029"])
        fit["ratio_msun_per_mpc"] = 999_999.0
        fit["size_uas"] = bhs.shadow_size_uas(999_999.0, 1.0)
        candidate = bhs.attach_nonadopted_ratio_fit(dict(baseline), fit)

        self.assertNotEqual(
            candidate["nonadopted_fit_size_uas"],
            self.by_key["Mrk1029"]["nonadopted_fit_size_uas"],
        )
        for field in bhs.OUTPUT_FIELDS:
            if field.startswith("nonadopted_fit_") or field == "uncertainty_status":
                continue
            self.assertEqual(candidate[field], baseline[field], field)

    def test_ngc1194_nonadopted_fit_preserves_the_adopted_result(self) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "NGC1194")
        baseline = bhs.build_shared_row(source)
        fit = dict(self.nonadopted_fits["NGC1194"])
        fit["ratio_msun_per_mpc"] = 999_999.0
        fit["size_uas"] = bhs.shadow_size_uas(999_999.0, 1.0)
        candidate = bhs.attach_nonadopted_ratio_fit(dict(baseline), fit)

        self.assertNotEqual(
            candidate["nonadopted_fit_size_uas"],
            self.by_key["NGC1194"]["nonadopted_fit_size_uas"],
        )
        for field in bhs.OUTPUT_FIELDS:
            if field.startswith("nonadopted_fit_"):
                continue
            self.assertEqual(candidate[field], baseline[field], field)

        self.assertEqual(self.by_key["NGC1194"]["size_uas"], 0.12533097499265997)
        self.assertEqual(
            self.by_key["NGC1194"]["size_unc_adopted_uas"],
            0.0077126753841636909,
        )
        self.assertEqual(
            self.by_key["NGC1194"]["mass_msun"], 64208389.715832204
        )
        self.assertEqual(
            self.by_key["NGC1194"]["ratio_provenance"],
            "published_mass_normalization_over_its_reference_distance",
        )
        self.assertEqual(
            self.by_key["NGC1194"]["uncertainty_status"],
            "shared distance canceled; reported undecomposed mass uncertainty propagated",
        )

    def test_ngc4388_nonadopted_fit_preserves_every_adopted_field(self) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "NGC4388")
        baseline = bhs.build_shared_row(source)
        fit = dict(self.nonadopted_fits["NGC4388"])
        fit["ratio_msun_per_mpc"] = 999_999.0
        fit["size_uas"] = bhs.shadow_size_uas(999_999.0, 1.0)
        candidate = bhs.attach_nonadopted_ratio_fit(dict(baseline), fit)

        self.assertNotEqual(
            candidate["nonadopted_fit_size_uas"],
            self.by_key["NGC4388"]["nonadopted_fit_size_uas"],
        )
        for field in bhs.OUTPUT_FIELDS:
            if field.startswith("nonadopted_fit_"):
                continue
            self.assertEqual(candidate[field], baseline[field], field)

    def test_ngc6264_nonadopted_fit_preserves_every_adopted_field(self) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "NGC6264")
        baseline = bhs.build_paired_row(source, self.quantiles["NGC6264"])
        fit = dict(self.nonadopted_fits["NGC6264"])
        fit["ratio_msun_per_mpc"] = 999_999.0
        fit["size_uas"] = bhs.shadow_size_uas(999_999.0, 1.0)
        candidate = bhs.attach_nonadopted_ratio_fit(dict(baseline), fit)

        self.assertNotEqual(
            candidate["nonadopted_fit_size_uas"],
            self.by_key["NGC6264"]["nonadopted_fit_size_uas"],
        )
        for field in bhs.OUTPUT_FIELDS:
            if field.startswith("nonadopted_fit_"):
                continue
            self.assertEqual(candidate[field], baseline[field], field)

        row = self.by_key["NGC6264"]
        self.assertEqual(row["size_uas"], 0.021455952583840227)
        self.assertEqual(row["size_stat_unc_minus_uas"], 5.196152422706768e-05)
        self.assertEqual(row["size_stat_unc_plus_uas"], 5.715767664977306e-05)
        self.assertEqual(
            row["size_systematic_unc_plus_uas"], 0.00021455952583840227
        )
        self.assertEqual(row["mass_msun"], 27_600_000.0)
        self.assertEqual(row["distance_mpc"], 132.1)
        self.assertEqual(
            row["nonadopted_fit_size_uas"], 0.021808206662753567
        )
        self.assertEqual(
            row["nonadopted_fit_size_unc_minus_uas"],
            0.00007180500782656395,
        )
        self.assertEqual(
            row["nonadopted_fit_size_unc_plus_uas"],
            0.00007180500782656395,
        )
        self.assertEqual(
            row["nonadopted_fit_basis"],
            "older_partial_data_reanalysis_conditional_on_fixed_vertical_center_"
            "and_selected_scatter_angle",
        )
        self.assertAlmostEqual(
            100.0 * (row["nonadopted_fit_size_uas"] / row["size_uas"] - 1.0),
            1.6417545552307233,
        )

    def test_ugc3789_nonadopted_fit_preserves_every_adopted_field(self) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "UGC3789")
        baseline = bhs.build_paired_row(source, self.quantiles["UGC3789"])
        fit = dict(self.nonadopted_fits["UGC3789"])
        fit["ratio_msun_per_mpc"] = 999_999.0
        fit["size_uas"] = bhs.shadow_size_uas(999_999.0, 1.0)
        candidate = bhs.attach_nonadopted_ratio_fit(dict(baseline), fit)

        self.assertNotEqual(
            candidate["nonadopted_fit_size_uas"],
            self.by_key["UGC3789"]["nonadopted_fit_size_uas"],
        )
        for field in bhs.OUTPUT_FIELDS:
            if field.startswith("nonadopted_fit_"):
                continue
            self.assertEqual(candidate[field], baseline[field], field)

        row = self.by_key["UGC3789"]
        self.assertEqual(row["size_uas"], 0.023705886582872197)
        self.assertEqual(row["size_stat_unc_minus_uas"], 7.794228634059805e-05)
        self.assertEqual(row["size_stat_unc_plus_uas"], 8.417766924784936e-05)
        self.assertEqual(
            row["size_systematic_unc_plus_uas"], 0.00023705886582872198
        )
        self.assertEqual(row["mass_msun"], 11_900_000.0)
        self.assertEqual(row["distance_mpc"], 51.5)
        self.assertEqual(row["nonadopted_fit_size_uas"], 0.02335714326015516)
        self.assertEqual(
            row["nonadopted_fit_size_unc_minus_uas"],
            0.00009232072434843937,
        )
        self.assertEqual(
            row["nonadopted_fit_size_unc_plus_uas"],
            0.00009232072434843937,
        )
        self.assertEqual(
            row["nonadopted_fit_basis"],
            "same_observation_vlbi_spot_subset_reanalysis_conditional_on_fixed_"
            "reid2013_inclination_and_selected_scatter_angle",
        )
        self.assertAlmostEqual(
            100.0 * (row["nonadopted_fit_size_uas"] / row["size_uas"] - 1.0),
            -1.47112541645672,
        )

    def test_eso558_nonadopted_fit_preserves_every_adopted_field(self) -> None:
        source = next(
            row for row in self.inputs if row["source_key"] == "ESO558-G009"
        )
        baseline = bhs.build_shared_row(source)
        fit = dict(self.nonadopted_fits["ESO558-G009"])
        fit["ratio_msun_per_mpc"] = 999_999.0
        fit["size_uas"] = bhs.shadow_size_uas(999_999.0, 1.0)
        candidate = bhs.attach_nonadopted_ratio_fit(dict(baseline), fit)

        self.assertNotEqual(
            candidate["nonadopted_fit_size_uas"],
            self.by_key["ESO558-G009"]["nonadopted_fit_size_uas"],
        )
        for field in bhs.OUTPUT_FIELDS:
            if field.startswith("nonadopted_fit_"):
                continue
            self.assertEqual(candidate[field], baseline[field], field)

        row = self.by_key["ESO558-G009"]
        self.assertEqual(row["size_uas"], 0.01620665338623987)
        self.assertEqual(row["size_total_unc_minus_uas"], 0.001041900020436918)
        self.assertEqual(row["size_total_unc_plus_uas"], 0.001041900020436918)
        self.assertEqual(row["mass_msun"], 16_792_963.464140728)
        self.assertEqual(row["distance_mpc"], 106.28958051420837)
        self.assertEqual(row["nonadopted_fit_size_uas"], 0.016351026067934704)
        self.assertEqual(
            row["nonadopted_fit_size_unc_minus_uas"],
            0.00026670431478438035,
        )
        self.assertEqual(
            row["nonadopted_fit_size_unc_plus_uas"],
            0.00025644645652344266,
        )
        self.assertEqual(
            row["nonadopted_fit_basis"],
            "same_data_reanalysis_conditional_on_fixed_gao2017_inclination_"
            "and_selected_scatter_angle",
        )
        self.assertAlmostEqual(
            100.0 * (row["nonadopted_fit_size_uas"] / row["size_uas"] - 1.0),
            0.8908235294117661,
        )

    def test_ngc2273_nonadopted_fit_preserves_the_adopted_result(self) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "NGC2273")
        baseline = bhs.build_shared_row(source)
        fit = dict(self.nonadopted_fits["NGC2273"])
        fit["ratio_msun_per_mpc"] = 999_999.0
        fit["size_uas"] = bhs.shadow_size_uas(999_999.0, 1.0)
        candidate = bhs.attach_nonadopted_ratio_fit(dict(baseline), fit)

        self.assertNotEqual(
            candidate["nonadopted_fit_size_uas"],
            self.by_key["NGC2273"]["nonadopted_fit_size_uas"],
        )
        for field in bhs.OUTPUT_FIELDS:
            if field.startswith("nonadopted_fit_"):
                continue
            self.assertEqual(candidate[field], baseline[field], field)

        row = self.by_key["NGC2273"]
        self.assertEqual(row["size_uas"], 0.029935384029973856)
        self.assertEqual(row["size_unc_adopted_uas"], 0.0019956922686649235)
        self.assertEqual(row["mass_msun"], 7408660.3518267926)
        self.assertEqual(
            row["ratio_provenance"],
            "published_mass_normalization_over_its_reference_distance",
        )
        self.assertEqual(
            row["uncertainty_status"],
            "shared distance canceled; reported undecomposed mass uncertainty propagated",
        )
        self.assertEqual(
            row["parameter_reference_doi"], "10.1088/0004-637X/727/1/20"
        )
        self.assertEqual(
            row["mass_uncertainty_reference_doi"],
            "10.3847/1538-4357/834/1/52",
        )

    def test_nonadopted_fit_attachment_fails_closed(self) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "NGC1194")
        baseline = bhs.build_shared_row(source)
        fit = dict(self.nonadopted_fits["NGC1194"])

        mismatched = dict(fit)
        mismatched["source_key"] = "Mrk1029"
        with self.assertRaisesRegex(ValueError, "does not match output source"):
            bhs.attach_nonadopted_ratio_fit(dict(baseline), mismatched)

        missing_primary = dict(baseline)
        missing_primary["mass_msun"] = None
        with self.assertRaisesRegex(ValueError, "requires the retained primary"):
            bhs.attach_nonadopted_ratio_fit(missing_primary, fit)

        already_attached = bhs.attach_nonadopted_ratio_fit(dict(baseline), fit)
        with self.assertRaisesRegex(ValueError, "already attached"):
            bhs.attach_nonadopted_ratio_fit(already_attached, fit)

    def test_shared_distance_cancels_exactly_from_shadow_ratio(self) -> None:
        for source in self.inputs:
            if source["treatment"] != "shared_distance":
                continue
            original = bhs.build_shared_row(source)
            changed = dict(source)
            changed["adopted_distance_mpc"] = str(
                float(source["adopted_distance_mpc"]) * 1.23456789
            )
            rebuilt = bhs.build_shared_row(changed)
            self.assertTrue(
                math.isclose(
                    original["size_uas"], rebuilt["size_uas"], rel_tol=2e-15, abs_tol=0.0
                ),
                source["source_key"],
            )
            self.assertIsNone(original["size_unc_independent_uas"])
            self.assertIsNone(original["size_unc_correlated_uas"])

    def test_ic2560_uses_coefficient_error_not_distance_error(self) -> None:
        row = self.by_key["IC2560"]
        expected_size = bhs.shadow_size_uas(4.14e6, 44.5)
        expected_stat = expected_size * 0.17 / 4.14
        self.assertEqual(row["size_uas"], expected_size)
        self.assertTrue(
            math.isclose(
                row["size_stat_unc_minus_uas"], expected_stat, rel_tol=2e-16, abs_tol=0.0
            )
        )
        self.assertTrue(
            math.isclose(
                row["size_stat_unc_plus_uas"], expected_stat, rel_tol=2e-16, abs_tol=0.0
            )
        )
        self.assertIsNone(row["size_systematic_unc_minus_uas"])
        self.assertIsNone(row["size_total_unc_minus_uas"])
        self.assertEqual(row["size_systematic_status"], "unquantified_model")
        self.assertEqual(
            row["size_total_status"], "incomplete_unquantified_systematic"
        )
        self.assertIn("model systematic unquantified", row["uncertainty_status"])
        self.assertAlmostEqual(row["distance_unc_mpc"], 6.1, places=15)
        self.assertAlmostEqual(row["size_unc_adopted_uas"], expected_stat, places=18)

    def test_ic2560_reviewed_contract_is_exact_and_fails_closed(self) -> None:
        source = next(row for row in self.inputs if row["source_key"] == "IC2560")
        row = self.by_key["IC2560"]
        expected_size = bhs.shadow_size_uas(4_140_000.0, 44.5)
        expected_stat = expected_size * 170_000.0 / 4_140_000.0
        expected_mass_distance = 4_140_000.0 * 6.1 / 44.5
        expected_mass_available = math.hypot(170_000.0, expected_mass_distance)

        self.assertEqual(float(source["reference_mass_msun"]), 4_140_000.0)
        self.assertEqual(float(source["reference_mass_unc_minus_msun"]), 170_000.0)
        self.assertEqual(float(source["reference_mass_unc_plus_msun"]), 170_000.0)
        self.assertEqual(float(source["reference_distance_mpc"]), 44.5)
        self.assertEqual(float(source["adopted_distance_mpc"]), 44.5)
        self.assertEqual(float(source["distance_measurement_unc_minus_mpc"]), 6.1)
        self.assertEqual(float(source["distance_measurement_unc_plus_mpc"]), 6.1)
        self.assertEqual(source["systematic_status"], "unquantified_model")
        self.assertEqual(
            source["uncertainty_note"],
            "formal_fit_error_known_no_separately_quantified_model_"
            "systematic_available",
        )
        self.assertEqual(
            source["mass_reference_doi"], "10.1017/S1743921323002971"
        )
        self.assertEqual(row["size_uas"], expected_size)
        self.assertEqual(row["size_stat_unc_minus_uas"], expected_stat)
        self.assertEqual(row["size_stat_unc_plus_uas"], expected_stat)
        self.assertEqual(row["mass_fit_unc_minus_msun"], 170_000.0)
        self.assertEqual(row["mass_fit_unc_plus_msun"], 170_000.0)
        self.assertEqual(
            row["mass_shared_distance_unc_minus_msun"], expected_mass_distance
        )
        self.assertEqual(
            row["mass_shared_distance_unc_plus_msun"], expected_mass_distance
        )
        self.assertEqual(row["mass_total_unc_minus_msun"], expected_mass_available)
        self.assertEqual(row["mass_total_unc_plus_msun"], expected_mass_available)
        self.assertEqual(row["mass_display_unc_minus_msun"], expected_mass_available)
        self.assertEqual(row["mass_display_unc_plus_msun"], expected_mass_available)
        self.assertIsNone(row["mass_systematic_unc_minus_msun"])
        self.assertIsNone(row["mass_systematic_unc_plus_msun"])
        self.assertIsNone(row["size_systematic_unc_minus_uas"])
        self.assertIsNone(row["size_total_unc_minus_uas"])
        self.assertNotIn("IC2560", self.model_checks)
        self.assertNotIn("IC2560", self.nonadopted_fits)

        source_mutations = (
            ("reference_mass_msun", "4140001"),
            ("reference_mass_unc_minus_msun", "170001"),
            ("reference_mass_unc_plus_msun", "170001"),
            ("reference_distance_mpc", "44.6"),
            ("adopted_distance_mpc", "44.6"),
            ("distance_measurement_unc_minus_mpc", "6.2"),
            ("distance_measurement_unc_plus_mpc", "6.2"),
            ("distance_flow_unc_mpc", "1"),
            ("distance_shared_systematic_unc_mpc", "1"),
            ("size_stat_fraction_minus_override", "0.04"),
            ("size_systematic_fraction_minus", "0.05"),
            ("reference_mass_systematic_unc_minus_msun", "100000"),
            ("systematic_status", "none_adopted"),
            ("mass_reference_doi", "10.0000/wrong"),
            ("mass_uncertainty_reference_doi", "10.0000/wrong"),
            ("distance_reference_doi", "10.0000/wrong"),
            ("uncertainty_note", "different_note"),
            ("mass_display_interval_status", "source_reported_total"),
            ("distance_stat_component_status", "not_available"),
            ("mass_systematic_status", "source_reported_absolute_position_systematic"),
        )
        for field, value in source_mutations:
            with self.subTest(source_field=field), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(item) for item in self.inputs]
                next(item for item in rows if item["source_key"] == "IC2560")[
                    field
                ] = value
                path = self._write_source_inputs(tmpdir, rows)
                with self.assertRaises(ValueError):
                    bhs.load_source_inputs(path)

        output_mutations = (
            ("size_uas", expected_size + 1e-12),
            ("size_stat_unc_minus_uas", expected_stat + 1e-12),
            ("mass_total_unc_minus_msun", expected_mass_available + 1.0),
            ("mass_display_unc_plus_msun", expected_mass_available + 1.0),
            ("size_systematic_status", "none_adopted"),
            ("size_total_status", "formal_component_only"),
            ("mass_systematic_unc_minus_msun", 1.0),
            ("size_total_unc_minus_uas", 1e-4),
            ("size_model_check_key", "unreviewed_check"),
            ("nonadopted_fit_key", "withdrawn_preprint"),
        )
        for field, value in output_mutations:
            with self.subTest(output_field=field):
                rows = [dict(item) for item in self.rows]
                next(item for item in rows if item["source_key"] == "IC2560")[
                    field
                ] = value
                with self.assertRaisesRegex(ValueError, "IC2560"):
                    bhs.render_csv(rows)

    def test_component_semantics_are_not_overstated(self) -> None:
        self.assertEqual(
            self.by_key["NGC1194"]["size_stat_component_status"],
            "reported_undecomposed_uncertainty",
        )
        self.assertIn(
            "not_probability_interval", self.by_key["NGC3393"]["size_total_status"]
        )
        self.assertEqual(
            self.by_key["ESO558-G009"]["size_stat_component_status"],
            "fit_plus_absolute_position_mixed_component",
        )
        self.assertEqual(
            self.by_key["NGC4388"]["size_stat_component_status"],
            "fit_plus_absolute_position_mixed_component",
        )
        self.assertEqual(
            self.by_key["ESO558-G009"]["size_total_status"],
            "quadrature_mixed_fit_position_plus_systematic",
        )
        self.assertIn(
            "mixed fit-plus-absolute-position",
            self.by_key["ESO558-G009"]["uncertainty_status"],
        )
        for key in ("NGC1068", "NGC4388", "Circinus", "IC2560", "J0437"):
            self.assertIsNone(self.by_key[key]["size_total_unc_minus_uas"], key)
            self.assertEqual(
                self.by_key[key]["size_total_status"],
                "incomplete_unquantified_systematic",
            )

    def test_circinus_display_status_inputs_and_stale_component_tokens_fail_closed(self) -> None:
        mutations = (
            (
                "mass_display_interval_status",
                "unknown",
                "invalid mass_display_interval_status",
            ),
            (
                "mass_display_interval_status",
                "decomposed_components",
                "source-reported mass display totals",
            ),
            (
                "distance_stat_component_status",
                "unknown",
                "invalid distance_stat_component_status",
            ),
            (
                "distance_stat_component_status",
                "statistical",
                "reported legacy multi-indicator distance status",
            ),
            (
                "reference_mass_unc_minus_msun",
                "300001",
                "does not match the reviewed input",
            ),
            (
                "size_stat_fraction_plus_override",
                "0.02",
                "does not match the reviewed input",
            ),
        )
        for field, value, message in mutations:
            with self.subTest(field=field, value=value), tempfile.TemporaryDirectory() as tmpdir:
                rows = [dict(row) for row in self.inputs]
                next(row for row in rows if row["source_key"] == "Circinus")[field] = value
                path = Path(tmpdir) / "inputs.csv"
                with path.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=bhs.INPUT_FIELDS)
                    writer.writeheader()
                    writer.writerows(rows)
                with self.assertRaisesRegex(ValueError, message):
                    bhs.load_source_inputs(path)

        stale = [dict(row) for row in self.rows]
        next(row for row in stale if row["source_key"] == "Circinus")[
            "size_stat_component_status"
        ] = "angular_radius_formal_term"
        with self.assertRaisesRegex(ValueError, "unsupported size statistical"):
            bhs.render_csv(stale)

    def test_mass_uncertainty_provenance_is_explicit_and_not_misattributed(self) -> None:
        input_by_key = {row["source_key"]: row for row in self.inputs}
        gao_doi = "10.3847/1538-4357/834/1/52"
        kuo_doi = "10.1088/0004-637X/727/1/20"

        for key in ("NGC1194", "NGC2273"):
            self.assertEqual(input_by_key[key]["mass_reference_doi"], kuo_doi)
            self.assertEqual(
                input_by_key[key]["mass_uncertainty_reference_doi"], gao_doi
            )
            self.assertEqual(self.by_key[key]["parameter_reference_doi"], kuo_doi)
            self.assertEqual(
                self.by_key[key]["mass_uncertainty_reference_doi"], gao_doi
            )
            self.assertEqual(
                self.by_key[key]["size_stat_component_status"],
                "reported_undecomposed_uncertainty",
            )
            self.assertEqual(
                input_by_key[key]["systematic_status"],
                "reported_undecomposed_uncertainty",
            )
            self.assertEqual(
                self.by_key[key]["size_total_status"],
                "reported_undecomposed_uncertainty_component_only",
            )
            self.assertEqual(
                self.by_key[key]["legacy_uncertainty_basis"],
                "reported_undecomposed_uncertainty_scalar",
            )
            self.assertIn("reported undecomposed", self.by_key[key]["uncertainty_status"])
            self.assertNotIn("statistical", self.by_key[key]["uncertainty_status"])
            self.assertNotIn("residual", self.by_key[key]["uncertainty_status"])

    def test_obsolete_residual_status_fails_closed(self) -> None:
        rows = [dict(row) for row in self.inputs]
        rows[1]["systematic_status"] = "conservative_undecomposed_residual"
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "inputs.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=bhs.INPUT_FIELDS)
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(ValueError, "invalid systematic_status"):
                bhs.load_source_inputs(path)

        self.assertNotIn("reference_mass_stat_unc_minus_msun", bhs.INPUT_FIELDS)
        for row in self.rows:
            provenance = row["mass_uncertainty_reference_doi"]
            if provenance is not None:
                self.assertTrue(str(provenance).startswith("10."), row["source_key"])
        self.assertEqual(
            {
                row["source_key"]
                for row in self.rows
                if row["mass_uncertainty_reference_doi"] is None
            },
            {"NGC4945", "NGC3079", "Mrk1029", "NGC1320"},
        )

    def test_two_sources_have_no_primary_mass_or_shadow(self) -> None:
        for key in ("Mrk1029", "NGC1320"):
            row = self.by_key[key]
            self.assertIsNone(row["mass_msun"])
            self.assertIsNone(row["size_uas"])
            self.assertEqual(row["adopted_uncertainty_model"], "not_available")
        for key in ("NGC4945", "NGC3079"):
            self.assertEqual(
                self.by_key[key]["size_stat_component_status"],
                "no_formal_uncertainty",
            )

    def test_six_paired_one_percent_systematics_only(self) -> None:
        one_percent = {
            row["source_key"]
            for row in self.rows
            if row["size_systematic_status"]
            == "numeric_one_percent_non_keplerian_allowance"
        }
        self.assertEqual(
            one_percent,
            {"NGC4258", "NGC5765b", "CGCG074-064", "NGC6264", "UGC3789", "NGC6323"},
        )

    def test_output_schema_and_deterministic_render(self) -> None:
        rendered_once = bhs.render_csv(self.rows)
        rendered_twice = bhs.render_csv(bhs.build_rows())
        self.assertEqual(rendered_once, rendered_twice)
        self.assertTrue(rendered_once.endswith("\n"))
        reader = csv.DictReader(io.StringIO(rendered_once))
        self.assertEqual(reader.fieldnames, bhs.OUTPUT_FIELDS)
        self.assertEqual(len(list(reader)), 21)

    def test_checked_in_output_is_exact_build(self) -> None:
        expected = bhs.render_csv(self.rows)
        self.assertEqual(bhs.DEFAULT_OUTPUT.read_text(encoding="utf-8"), expected)

    def test_missing_paired_source_is_rejected(self) -> None:
        with bhs.DEFAULT_QUANTILES.open(newline="", encoding="utf-8") as handle:
            quantile_rows = list(csv.DictReader(handle))
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "quantiles.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=bhs.QUANTILE_FIELDS)
                writer.writeheader()
                writer.writerows(quantile_rows[:-1])
            with self.assertRaisesRegex(ValueError, "exactly match"):
                bhs.build_rows(quantile_path=path)

    def test_formal_only_build_and_final_check_succeed(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            candidate = Path(tmpdir) / "candidate.csv"
            self.assertEqual(
                bhs.main(["--output", str(candidate), "--require-final"]),
                0,
            )
            self.assertTrue(candidate.is_file())
            self.assertEqual(
                bhs.main(
                    [
                        "--output",
                        str(candidate),
                        "--check",
                        "--require-final",
                    ]
                ),
                0,
            )

    def test_final_release_gate_accepts_explicit_unquantified_model_status(self) -> None:
        self.assertIsNone(bhs.require_final_release(self.rows))


if __name__ == "__main__":
    unittest.main()
