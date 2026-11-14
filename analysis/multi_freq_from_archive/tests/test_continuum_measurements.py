from __future__ import annotations

import ast
import csv
import hashlib
import json
import math
import subprocess
import sys
import unittest
from pathlib import Path

import numpy as np


PACKAGE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_DIR))

from continuum_measurements import (  # noqa: E402
    FINAL_NOISE_MIN_INDEPENDENT_BEAMS,
    SCALED_MAD_NORMAL_FACTOR,
    adopted_noise_snr,
    angular_value_to_arcsec,
    background_rms_snr,
    calendar_date_sort_key,
    classify_ordinary_fit,
    effective_independent_noise_samples,
    final_noise_review_flags,
    forced_fit_diagnostic,
    format_two_decimal,
    gaussian_beam_area_arcsec2,
    independent_beam_areas,
    normalize_right_ascension_degrees,
    portable_release_product,
    require_fixed_restoring_beam_shape,
    require_minimum_independent_beams,
    square_aperture_coverage,
    square_bmaj_region_width_arcsec,
)
from run_continuum_four_pass import (  # noqa: E402
    compare_reviewed_current_to_generated,
)


class ContinuumMeasurementTests(unittest.TestCase):
    def test_calendar_date_sort_key_ignores_suffix_and_rejects_bad_dates(self) -> None:
        self.assertLess(
            calendar_date_sort_key("1996/12/22_1"),
            calendar_date_sort_key("1997/01/05_1"),
        )
        self.assertEqual(
            calendar_date_sort_key("2019/01/18_1"),
            calendar_date_sort_key("2019/01/18_2"),
        )
        self.assertEqual(calendar_date_sort_key("NA"), (0, 0, 0))
        for invalid in ("1997/02/29_1", "1997/01/05_bad", "1997-01-05_1"):
            with self.assertRaises(ValueError):
                calendar_date_sort_key(invalid)

    def test_active_fitter_does_not_parse_suffixed_day_as_a_float(self) -> None:
        source = (PACKAGE_DIR / "get_continuum_params.py").read_text(encoding="utf-8")
        self.assertIn("calendar_date_sort_key(x[4])", source)
        self.assertNotIn("float(x[4].split('/')[2])", source)

    def test_formatter_preserves_zero_digits_and_rounds_normally(self) -> None:
        self.assertEqual(format_two_decimal(0.02063), "0.02")
        self.assertEqual(format_two_decimal(0.50607), "0.51")
        self.assertEqual(format_two_decimal(0.04072), "0.04")
        self.assertEqual(format_two_decimal(0.3053881258), "0.31")
        self.assertEqual(format_two_decimal(12), "12.00")

    def test_native_right_ascension_wrap_is_normalized_not_rejected(self) -> None:
        self.assertAlmostEqual(
            normalize_right_ascension_degrees(-149.2314080204706),
            210.7685919795294,
        )
        self.assertEqual(normalize_right_ascension_degrees(360.0), 0.0)
        self.assertEqual(normalize_right_ascension_degrees(720.25), 0.25)
        for invalid in (float("nan"), float("inf"), float("-inf")):
            with self.assertRaises(ValueError):
                normalize_right_ascension_degrees(invalid)

    def test_discovery_rms_and_final_adopted_noise_snr_are_separate(self) -> None:
        self.assertTrue(math.isclose(background_rms_snr(12.41, 1.84), 12.41 / 1.84))
        self.assertTrue(math.isclose(adopted_noise_snr(12.41, 1.70), 12.41 / 1.70))
        for invalid_rms in (0, -1, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                background_rms_snr(1, invalid_rms)

    def test_forced_fit_diagnostic_never_fabricates_a_nonconverged_flux(self) -> None:
        nonconverged = {
            "converged": [False],
            "results": {"nelements": 0},
        }
        self.assertEqual(
            forced_fit_diagnostic(nonconverged),
            (False, None, None),
        )

        converged = {
            "converged": [True],
            "results": {
                "nelements": 1,
                "component0": {
                    "flux": {"value": [0.00125], "unit": "Jy"}
                },
            },
        }
        self.assertEqual(
            forced_fit_diagnostic(converged),
            (True, 0.00125, "Jy"),
        )
        with_deconvolved = dict(converged)
        with_deconvolved["deconvolved"] = {
            "component0": {
                "flux": {"value": [0.00125], "unit": "Jy"}
            }
        }
        self.assertEqual(
            forced_fit_diagnostic(with_deconvolved),
            (True, 0.00125, "Jy"),
        )
        numpy_flag = dict(converged)
        numpy_flag["converged"] = np.array([True], dtype=np.bool_)
        self.assertEqual(
            forced_fit_diagnostic(numpy_flag),
            (True, 0.00125, "Jy"),
        )

        for malformed in (
            {},
            {"converged": []},
            {"converged": ["False"], "results": {"nelements": 0}},
            {"converged": [False], "results": {"nelements": 0.0}},
            {"converged": [True], "results": {"nelements": 0}},
            {
                "converged": [True],
                "results": {
                    "nelements": 1,
                    "component0": {
                        "flux": {"value": [0.00125], "unit": None}
                    },
                },
            },
            {
                "converged": [True],
                "results": {
                    "nelements": 1,
                    "component0": {
                        "flux": {"value": [0.00125], "unit": "Jy"}
                    },
                },
                "deconvolved": None,
            },
            {
                "converged": [False],
                "results": {"nelements": 1, "component0": {}},
            },
            {
                "converged": [True],
                "results": {"nelements": 2, "component0": {}},
            },
            {
                "converged": [True],
                "results": {
                    "component0": {
                        "flux": {"value": [float("nan")], "unit": "Jy"}
                    }
                },
            },
            {
                "converged": [True],
                "results": {
                    "nelements": 1,
                    "component0": {
                        "flux": {"value": [0.00125], "unit": "Jy"}
                    },
                },
                "deconvolved": {
                    "component0": {
                        "flux": {"value": [0.0025], "unit": "Jy"}
                    }
                },
            },
        ):
            with self.subTest(malformed=malformed), self.assertRaises(ValueError):
                forced_fit_diagnostic(malformed)

    def test_ordinary_fit_classification_reuses_one_adopted_noise_value(self) -> None:
        below = classify_ordinary_fit(2.9, 1.0)
        self.assertFalse(below.is_detection)
        self.assertEqual(below.snr, 2.9)
        self.assertEqual(below.upper_limit_mjy, 3.0)

        threshold = classify_ordinary_fit(3.0, 1.0)
        self.assertTrue(threshold.is_detection)
        self.assertEqual(threshold.snr, 3.0)
        self.assertEqual(threshold.upper_limit_mjy, 3.0)

        above = classify_ordinary_fit(4.0, 1.0)
        self.assertTrue(above.is_detection)
        self.assertEqual(above.upper_limit_mjy, 4.0)
        self.assertGreaterEqual(above.upper_limit_mjy, 4.0)

        for invalid_rms in (0, -1, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                classify_ordinary_fit(1, invalid_rms)

    def test_final_noise_area_conversion_and_minimum_guard_are_exact(self) -> None:
        self.assertAlmostEqual(angular_value_to_arcsec(math.pi / 648000, "rad"), 1.0)
        self.assertAlmostEqual(angular_value_to_arcsec(1 / 3600, "deg"), 1.0)
        beam_area = gaussian_beam_area_arcsec2(
            0.72726505994808,
            0.66569519042964,
        )
        self.assertAlmostEqual(beam_area, 0.548570643450105)
        cgcg_beam_areas = independent_beam_areas(
            2577,
            0.14000000009212799,
            0.14000000009212799,
            0.72726505994808,
            0.66569519042964,
        )
        self.assertAlmostEqual(cgcg_beam_areas, 92.07419439875643)
        cgcg = classify_ordinary_fit(
            0.22514342163586162,
            0.034753733841365486,
        )
        self.assertTrue(cgcg.is_detection)
        self.assertEqual(cgcg.snr, 6.478251305702457)
        self.assertEqual(f"{cgcg.snr:.3f}", "6.478")
        self.assertEqual(format_two_decimal(0.22514342163586162), "0.23")
        self.assertEqual(
            effective_independent_noise_samples(100, 2, 2, 1, 1),
            100,
        )
        self.assertEqual(
            require_minimum_independent_beams(FINAL_NOISE_MIN_INDEPENDENT_BEAMS),
            FINAL_NOISE_MIN_INDEPENDENT_BEAMS,
        )
        with self.assertRaisesRegex(ValueError, "at least 50"):
            require_minimum_independent_beams(49.999)

    def test_final_fixed_beam_guard_catches_shape_mutations(self) -> None:
        check = require_fixed_restoring_beam_shape(
            0.7272650601931019,
            0.6656951847265469,
            242.904620211229,
            0.72726505994808,
            0.66569519042964,
            62.90461730957,
        )
        self.assertLess(check.major_delta_arcsec, 1e-6)
        self.assertLess(check.minor_delta_arcsec, 1e-6)
        self.assertLess(check.pa_delta_deg, 1e-4)
        serialized = require_fixed_restoring_beam_shape(
            0.7272650601931019,
            0.6656951847265469,
            62.90504015957,
            0.72726505994808,
            0.66569519042964,
            62.90461730957,
        )
        self.assertAlmostEqual(serialized.pa_delta_deg, 0.00042285, places=8)
        self.assertTrue(serialized.pa_guard_applicable)
        mutations = (
            (0.73, 0.66569519042964, 62.90461730957),
            (0.72726505994808, 0.66, 62.90461730957),
            (0.72726505994808, 0.66569519042964, 63.0),
        )
        for major, minor, pa in mutations:
            with self.subTest(major=major, minor=minor, pa=pa):
                with self.assertRaisesRegex(ValueError, "not fixed"):
                    require_fixed_restoring_beam_shape(
                        major,
                        minor,
                        pa,
                        0.72726505994808,
                        0.66569519042964,
                        62.90461730957,
                    )

        circular = require_fixed_restoring_beam_shape(
            1.0,
            1.0,
            89.0,
            1.0,
            1.0,
            0.0,
        )
        self.assertFalse(circular.pa_guard_applicable)
        self.assertEqual(circular.pa_delta_deg, 89.0)

    def test_final_square_aperture_records_finite_image_truncation(self) -> None:
        coverage = square_aperture_coverage(
            49.5,
            49.5,
            10.0,
            1.0,
            1.0,
            100,
            100,
        )
        self.assertEqual(coverage.minimum_edge_clearance_pixels, 45.0)
        # A finite image may truncate the requested square.  The exact
        # intersection is retained and flagged; N_eff remains the hard gate.
        clipped = square_aperture_coverage(
            3.0,
            49.5,
            10.0,
            1.0,
            1.0,
            100,
            100,
        )
        self.assertEqual(clipped.minimum_edge_clearance_pixels, -1.5)
        self.assertIn(
            "aperture_truncated_by_image",
            final_noise_review_flags(
                0,
                1,
                1,
                1,
                aperture_minimum_edge_clearance_pixels=(
                    clipped.minimum_edge_clearance_pixels
                ),
            ),
        )

    def test_final_noise_review_flags_retain_classic_statistics_as_diagnostics(self) -> None:
        self.assertEqual(final_noise_review_flags(0, 1, 1, 1), ())
        self.assertEqual(
            final_noise_review_flags(0.21, 1, 1, 1),
            ("nonzero_background_mean",),
        )
        self.assertEqual(
            final_noise_review_flags(0, 1, 1.21, 1),
            ("rms_scaled_mad_disagreement",),
        )
        with self.assertRaisesRegex(ValueError, "scaled MAD must be positive"):
            final_noise_review_flags(0, 1, 1, 0)
        self.assertAlmostEqual(SCALED_MAD_NORMAL_FACTOR, 1.482602218505602)

    def test_release_product_paths_are_portable_and_unknown_absolutes_fail(self) -> None:
        absolute_sma = (
            "/home/author/repo/analysis/sma_calibration/230316_128/images/clean/"
            "nitercal_10/flag2_cal1/fields/science/NGC1068_test.image"
        )
        self.assertEqual(
            portable_release_product(absolute_sma),
            "analysis/sma_calibration/230316_128/images/clean/"
            "nitercal_10/flag2_cal1/fields/science/NGC1068_test.image",
        )
        relative_alma = "ALMA/NGC4945/example.fits"
        self.assertEqual(portable_release_product(relative_alma), relative_alma)
        with self.assertRaises(ValueError):
            portable_release_product("/tmp/untraceable/example.fits")

    def test_search_region_is_an_explicit_square_derived_only_from_bmaj(self) -> None:
        self.assertAlmostEqual(square_bmaj_region_width_arcsec(0.2, 10), 2.0)
        self.assertAlmostEqual(
            square_bmaj_region_width_arcsec(0.02, 10, minimum_arcsec=1), 1.0
        )
        self.assertAlmostEqual(
            square_bmaj_region_width_arcsec(0.4, 20, maximum_arcsec=3), 3.0
        )
        self.assertAlmostEqual(
            square_bmaj_region_width_arcsec(
                0.4, 20, maximum_arcsec=3, override_arcsec=2.5
            ),
            2.5,
        )
        source = (PACKAGE_DIR / "get_continuum_params.py").read_text()
        self.assertNotIn("offset2tot", source)
        self.assertIn("square_bmaj_region_width_arcsec", source)

    def test_active_fitter_has_no_mixed_frame_or_quadrature_noise_path(self) -> None:
        source = (PACKAGE_DIR / "get_continuum_params.py").read_text()
        self.assertIn(
            "[aselectcoord1, aselectcoord2]",
            source,
        )
        self.assertIn(
            "[selectcoord1, selectcoord2]",
            source,
        )
        self.assertNotIn(
            "[aselectcoord1,aselectcoord2],[selectcoord1,aselectcoord2],bestrezfile",
            source,
        )
        self.assertIn("referencecode(type='direction')", source)
        self.assertIn("active_coords = coordinates_for_image(file, acoords, coords)", source)
        self.assertNotIn('direction = f"ICRS {ra} {dec}"', source)
        self.assertNotIn("np.sqrt(fluxerr**2 + noiseval**2)", source)
        self.assertIn("classification=classify_ordinary_fit(fluxval, noiseval)", source)
        self.assertEqual(source.count("background_rms_snr(fluxval, noiseval)"), 1)

    def test_ordinary_noise_uses_decoupled_final_aperture_and_reuses_its_limit(self) -> None:
        source = (PACKAGE_DIR / "get_continuum_params.py").read_text(encoding="utf-8")
        tree = ast.parse(source)
        calls = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "printrueconv"
        ]
        self.assertEqual(len(calls), 1)
        for call in calls:
            keywords = {keyword.arg: keyword.value for keyword in call.keywords}
            self.assertNotIn("local_fit_region", keywords)
            self.assertNotIn("replace_existing_final_noise_audit", keywords)

        fitter = next(
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name == "printrueconv"
        )
        fitter_source = ast.get_source_segment(source, fitter)
        ordinary_branch = next(
            node
            for node in ast.walk(fitter)
            if isinstance(node, ast.If)
            and isinstance(node.test, ast.UnaryOp)
            and isinstance(node.test.op, ast.Not)
            and isinstance(node.test.operand, ast.Name)
            and node.test.operand.id == "doublegaussfit"
            and any(
                isinstance(child, ast.Try)
                for statement in node.body
                for child in ast.walk(statement)
            )
        )
        ordinary_source = ast.get_source_segment(source, ordinary_branch)
        self.assertIn(
            "ordinary_noise_scope == 'publication_final_10bmaj_source_excluded'",
            ordinary_source,
        )
        self.assertIn("measure_final_publication_noise", ordinary_source)
        self.assertIn("require_minimum_independent_beams", ordinary_source)
        self.assertIn("totstats = imstat(imagename=file)", ordinary_source)
        self.assertIn("'intermediate_whole_image'", fitter_source)
        self.assertIn("classification.upper_limit_mjy", source)
        self.assertGreaterEqual(source.count("upbound=ordinary_upper_limit_mjy"), 2)
        self.assertIn("Cannot measure the source-excluded local noise", source)
        self.assertIn("Duplicate final-noise audit product", source)
        self.assertNotIn(
            'if audit_product in final_noise_audit_records',
            fitter_source,
        )

        measurement = next(
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef)
            and node.name == "measure_final_publication_noise"
        )
        measurement_source = ast.get_source_segment(source, measurement)
        self.assertIn("FINAL_NOISE_APERTURE_BMAJ_MULTIPLIER", measurement_source)
        self.assertIn("centerbox[[", measurement_source)
        self.assertIn("crtf_ellipse_axes_are_semiaxes", measurement_source)
        self.assertIn("square_aperture_coverage", measurement_source)
        self.assertIn("truncated_available_valid_intersection", measurement_source)
        self.assertIn(
            '"single_initial_fixed_beam_fit_no_widened_retry"',
            measurement_source,
        )
        self.assertNotIn("local_fit_region", measurement_source)

    def test_final_audit_links_pixel_and_native_celestial_fit_centers(self) -> None:
        source = (PACKAGE_DIR / "get_continuum_params.py").read_text(encoding="utf-8")
        helper = source[
            source.index("def _native_fitted_direction_degrees") :
            source.index("def measure_final_publication_noise")
        ]
        self.assertIn('["shape"]["direction"]', helper)
        self.assertIn("normalize_celestial_reference_code", helper)
        self.assertIn("image_celestial_frame(imagename)", helper)
        self.assertIn("component_frame != image_frame", helper)
        self.assertIn("normalize_right_ascension_degrees(ra_deg)", helper)
        self.assertIn("0.0 <= ra_deg < 360.0", helper)
        self.assertIn("-90.0 <= dec_deg <= 90.0", helper)
        for field in (
            '"fit_center_reference_frame": fit_center_frame',
            '"fit_center_ra_deg": repr(fit_center_ra_deg)',
            '"fit_center_dec_deg": repr(fit_center_dec_deg)',
        ):
            self.assertIn(field, source)

    def test_ordinary_final_fit_has_no_unreachable_widened_retry_branch(self) -> None:
        source = (PACKAGE_DIR / "get_continuum_params.py").read_text(encoding="utf-8")
        tree = ast.parse(source)
        fitter_calls = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "printrueconv"
        ]
        self.assertEqual(len(fitter_calls), 1)
        for forbidden in (
            "retry_fitter",
            "retry_pending",
            "replace_final_noise_audit",
            "replace_existing_final_noise_audit",
            "Authoritative widened retry",
        ):
            self.assertNotIn(forbidden, source)

    def test_forced_dirty_sma_bypasses_ordinary_classifier_and_area_guard(self) -> None:
        source = (PACKAGE_DIR / "get_continuum_params.py").read_text(encoding="utf-8")
        tree = ast.parse(source)
        fitter = next(
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name == "printrueconv"
        )
        fitter_source = ast.get_source_segment(source, fitter)
        force_offset = fitter_source.index("if force_upper_limit:")
        classifier_offset = fitter_source.index("classification=classify_ordinary_fit")
        self.assertLess(force_offset, classifier_offset)
        force_node = next(
            node
            for node in ast.walk(fitter)
            if isinstance(node, ast.If)
            and isinstance(node.test, ast.Name)
            and node.test.id == "force_upper_limit"
        )
        force_body_calls = {
            child.func.id
            for statement in force_node.body
            for child in ast.walk(statement)
            if isinstance(child, ast.Call) and isinstance(child.func, ast.Name)
        }
        self.assertNotIn("classify_ordinary_fit", force_body_calls)
        self.assertNotIn("require_minimum_independent_beams", force_body_calls)

    def test_detection_is_not_starred_by_its_plot_color(self) -> None:
        source = (PACKAGE_DIR / "get_continuum_params.py").read_text(encoding="utf-8")
        tree = ast.parse(source)
        fitter = next(
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name == "printrueconv"
        )
        ordinary_branch = next(
            node
            for node in ast.walk(fitter)
            if isinstance(node, ast.If)
            and isinstance(node.test, ast.UnaryOp)
            and isinstance(node.test.op, ast.Not)
            and isinstance(node.test.operand, ast.Name)
            and node.test.operand.id == "doublegaussfit"
            and any(
                isinstance(child, ast.Try)
                for statement in node.body
                for child in ast.walk(statement)
            )
        )
        ordinary_source = ast.get_source_segment(source, ordinary_branch)
        self.assertNotIn("sfluxval=fluxval+'*'", ordinary_source)
        self.assertEqual(source.count("low,date,ordinary_upper_limit_mjy=printrueconv"), 1)

    def test_forced_dirty_sma_limit_remains_the_one_beam_region_path(self) -> None:
        source = (PACKAGE_DIR / "get_continuum_params.py").read_text(encoding="utf-8")
        self.assertIn(
            "if force_upper_limit:\n"
            "                            # Preserve the reviewed dirty-SMA policy: its\n"
            "                            # limit is three times the one-beam-region RMS,\n"
            "                            # separate from the ordinary-fit classifier.\n"
            "                            forced_statistics=image_statistics_mjy_per_beam(\n"
            "                                file, active_regionsmall\n"
            "                            )\n"
            "                            noise=forced_statistics[\"rms_mjy_per_beam\"]\n"
            "                            upbound=noise*3",
            source,
        )
        self.assertIn("record_forced_sma_noise_audit", source)
        self.assertIn("forced_fit=forced_fit_diagnostic(imfitter)", source)
        self.assertIn("if forced_fit.converged:", source)
        self.assertIn("fit_converged=forced_fit.converged", source)
        self.assertNotIn("imfitter['deconvolved']['component0']", source)
        self.assertIn(
            "Nonconverged forced dirty-SMA fit must not carry a diagnostic flux",
            source,
        )
        self.assertIn(
            "if fit_converged is not True and fit_converged is not False:",
            source,
        )
        self.assertNotIn("type(fit_converged)", source)
        self.assertIn("math.isfinite(fitted_flux_value)", source)

    def test_released_correction_manifest_records_exact_scope(self) -> None:
        manifest = PACKAGE_DIR / "validated_continuum_corrections.csv"
        with manifest.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 84)
        self.assertEqual(
            sum(row["old_flux_field"] != row["new_flux_field"] for row in rows),
            68,
        )
        self.assertEqual(sum(row["old_snr"] != row["new_snr"] for row in rows), 72)
        reasons = "+".join(row["correction_reasons"] for row in rows)
        for reason in (
            "ordinary_two_decimal_formatting",
            "coherent_j2000_anchor_pair",
            "single_background_rms",
            "single_local_scaled_mad",
            "fixed_nuclear_component_source18",
            "source18_10bmaj_local_noise_aperture",
            "header_aware_target_frame",
            "header_aware_anchor_frame",
        ):
            self.assertIn(reason, reasons)
        self.assertFalse(any("/home/" in row["product"] for row in rows))
        ngc4945_reasons = {
            row["correction_reasons"]
            for row in rows
            if row["source"] == "NGC4945"
        }
        self.assertEqual(
            ngc4945_reasons,
            {
                "fixed_nuclear_component_source18+single_local_scaled_mad+"
                "source18_10bmaj_local_noise_aperture"
            },
        )
        self.assertEqual({row["schema_version"] for row in rows}, {"2"})

        historical_path = (
            PACKAGE_DIR
            / "fitsumfiles/historical/"
            "fitsummary_fixtochosen_sorted_409_reference.txt"
        )
        historical_rows = [
            ast.literal_eval(line)
            for line in historical_path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        self.assertEqual(len(historical_rows), 409)

        sorted_path = (
            PACKAGE_DIR / "fitsumfiles/sorted/fitsummary_fixtochosen_sorted.txt"
        )
        sorted_rows = [
            ast.literal_eval(line)
            for line in sorted_path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]

        historical_by_product = {
            portable_release_product(row[7]): row for row in historical_rows
        }
        live_by_product = {
            portable_release_product(row[7]): row for row in sorted_rows
        }
        legacy_rows = [row for row in rows if row["source"] != "NGC4945"]
        source18_rows = [row for row in rows if row["source"] == "NGC4945"]
        self.assertEqual(len(legacy_rows), 39)
        self.assertEqual(len(source18_rows), 45)

        # The 39 non-NGC4945 rows are immutable evidence from the correction
        # phase that produced the frozen 409-row table.  A later, globally
        # applied noise policy must not retroactively rewrite their values.
        for correction in legacy_rows:
            historical = historical_by_product[correction["product"]]
            self.assertEqual(historical[3], correction["new_flux_field"])
            self.assertEqual(historical[6], correction["new_snr"])

        # The NGC 4945 records are intentionally different: their schema-2 rows
        # are the replayable source-18 replacement and therefore describe the
        # live table produced by that replay.
        for correction in source18_rows:
            live = live_by_product[correction["product"]]
            self.assertEqual(live[3], correction["new_flux_field"])
            self.assertEqual(live[6], correction["new_snr"])

        # All final live changes relative to the frozen baseline are a separate
        # phase.  Require the promoted, product-complete reviewed ledger rather
        # than treating the cumulative correction manifest as the complete
        # inventory-and-admission ledger.
        adjudication = compare_reviewed_current_to_generated(
            sorted_path,
            historical_path,
            PACKAGE_DIR / "continuum_product_exclusion_policy.csv",
            PACKAGE_DIR / "reviewed_current_to_generated_deltas.csv",
            PACKAGE_DIR / "continuum_product_addition_policy.csv",
        )
        self.assertEqual(adjudication["status"], "exact_reviewed_match")
        self.assertEqual(adjudication["canonical_rows"], 409)
        self.assertEqual(adjudication["reviewed_exclusion_rows"], 3)
        self.assertEqual(adjudication["validated_product_identities"], 406)
        self.assertEqual(adjudication["validated_addition_identities"], 89)

    def test_ngc2960_icrs_upper_limit_audit_and_regression(self) -> None:
        audit = json.loads(
            (PACKAGE_DIR / "ngc2960_frame_correction_audit.json").read_text(
                encoding="utf-8"
            )
        )
        product = (
            "ALMA/NGC2960/"
            "member.uid___A001_X5a3_X167.ari_l.Mrk_1419_sci.spw0_1_2_3_"
            "85885MHz.12m.cont.I.pbcor.fits"
        )
        self.assertEqual(
            audit["status"],
            "historical_frame_rms_evidence_superseded_by_adjudicated_scaled_mad",
        )
        self.assertEqual(
            audit["measurement"]["scope"],
            "historical_pre_scaled_mad_frame_rms_decision",
        )
        downstream = audit["downstream_regeneration"]
        self.assertEqual(
            downstream["status"],
            "final_scaled_mad_release_and_manuscript_synchronized",
        )
        self.assertEqual(downstream["historical_flux_field"], "0.47*")

        sorted_path = (
            PACKAGE_DIR / "fitsumfiles/sorted/fitsummary_fixtochosen_sorted.txt"
        )
        sorted_rows = [
            ast.literal_eval(line)
            for line in sorted_path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        live_matches = [
            row
            for row in sorted_rows
            if portable_release_product(row[7]) == product
        ]
        self.assertEqual(len(live_matches), 1)
        live = live_matches[0]
        current_flux_field = str(live[3])
        current_snr = str(live[6])
        self.assertEqual(
            downstream["current_canonical_flux_field"], current_flux_field
        )
        self.assertFalse(current_flux_field.endswith("*"))
        self.assertNotEqual(current_snr, "NA")

        superseding_ledger = downstream["superseding_ledger"]
        self.assertEqual(set(superseding_ledger), {"path", "sha256"})
        self.assertEqual(
            superseding_ledger["path"],
            "reviewed_current_to_generated_deltas.csv",
        )
        ledger_path = PACKAGE_DIR / superseding_ledger["path"]
        self.assertTrue(ledger_path.is_file())
        self.assertEqual(
            hashlib.sha256(ledger_path.read_bytes()).hexdigest(),
            superseding_ledger["sha256"],
        )
        with ledger_path.open(encoding="utf-8", newline="") as handle:
            ledger_rows = list(csv.DictReader(handle))
        transition = [row for row in ledger_rows if row["product"] == product]
        self.assertEqual(len(transition), 1)
        self.assertEqual(transition[0]["current_flux_field"], "0.47*")
        self.assertEqual(transition[0]["generated_flux_field"], current_flux_field)
        self.assertEqual(transition[0]["current_snr"], "NA")
        self.assertEqual(transition[0]["generated_snr"], current_snr)
        self.assertIn("flux_field", transition[0]["changed_fields"].split("+"))

        required_release_products = [
            "fitsumfiles/sorted/fitsummary_fixtochosen_sorted.txt",
            "fitsumfiles/final/fitsummary_final_withextras.txt",
            "machinetables/fitsummary_machine.txt",
            "machinetables/fitsummary_aas_mrt.txt",
            "machinetables/AAS72796R2_datafile2_corrected.txt",
            "machinetables/correlations_machine.txt",
        ]
        self.assertEqual(
            downstream["required_release_products"], required_release_products
        )
        self.assertEqual(
            set(downstream["release_product_sha256"]),
            set(required_release_products),
        )
        for relative in required_release_products:
            expected_hash = downstream["release_product_sha256"][relative]
            actual_hash = hashlib.sha256(
                (PACKAGE_DIR / relative).read_bytes()
            ).hexdigest()
            self.assertEqual(actual_hash, expected_hash, relative)

        required_manuscript_products = ["tab2.txt", "datafile2.txt"]
        self.assertEqual(
            downstream["required_manuscript_products"],
            required_manuscript_products,
        )
        self.assertEqual(
            set(downstream["manuscript_product_sha256"]),
            set(required_manuscript_products),
        )
        for expected_hash in downstream["manuscript_product_sha256"].values():
            self.assertRegex(expected_hash, r"^[0-9a-f]{64}$")
        manuscript_dir = PACKAGE_DIR.parents[2] / "revision_package/working_clean"
        if manuscript_dir.is_dir():
            for relative in required_manuscript_products:
                expected_hash = downstream["manuscript_product_sha256"][relative]
                actual_hash = hashlib.sha256(
                    (manuscript_dir / relative).read_bytes()
                ).hexdigest()
                self.assertEqual(actual_hash, expected_hash, relative)
        self.assertIsNone(downstream["remaining_action"])
        self.assertEqual(audit["input"]["product"], product)
        self.assertEqual(audit["input"]["radesys"], "ICRS")
        self.assertEqual(audit["coordinate_frame_evidence"]["catalog_frame"], "J2000")

        old = audit["region_evidence"]["old_wrong_frame_region"]
        accepted = audit["region_evidence"]["accepted_native_frame_region"]
        self.assertEqual(old["pixel_count"], 182)
        self.assertEqual(accepted["pixel_count"], 169)
        self.assertTrue(
            math.isclose(
                3 * old["rms_mjy_per_beam"],
                old["three_sigma_mjy"],
                rel_tol=0,
                abs_tol=1e-15,
            )
        )
        self.assertTrue(
            math.isclose(
                3 * accepted["rms_mjy_per_beam"],
                accepted["three_sigma_mjy"],
                rel_tol=0,
                abs_tol=1e-15,
            )
        )
        self.assertEqual(format_two_decimal(old["three_sigma_mjy"]), "0.45")
        self.assertEqual(format_two_decimal(accepted["three_sigma_mjy"]), "0.47")
        self.assertEqual(audit["decision"]["accepted_flux_field"], "0.47*")
        self.assertFalse(audit["decision"]["classification_changed"])
        self.assertFalse(audit["decision"]["selected_s230_observation_changed"])

        with (PACKAGE_DIR / "validated_continuum_corrections.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle))
        correction = [row for row in rows if row["product"] == product]
        self.assertEqual(len(correction), 1)
        self.assertEqual(correction[0]["old_flux_field"], "0.45*")
        self.assertEqual(correction[0]["new_flux_field"], "0.47*")
        self.assertEqual(correction[0]["correction_reasons"], "header_aware_target_frame")
        self.assertTrue(
            math.isclose(
                float(correction[0]["validated_unrounded_flux_mjy"]),
                accepted["three_sigma_mjy"],
                rel_tol=0,
                abs_tol=1e-15,
            )
        )

        source = (PACKAGE_DIR / "get_continuum_params.py").read_text(encoding="utf-8")
        self.assertIn("measure_final_publication_noise", source)
        self.assertNotIn(
            "noise=image_rms_mjy_per_beam(file, regiontot)",
            source,
        )

    def test_ngc4945_fixed_nucleus_sidecars_are_complete_and_consistent(self) -> None:
        photometry_path = PACKAGE_DIR / "ngc4945_nuclear_photometry.csv"
        with photometry_path.open(encoding="utf-8", newline="") as handle:
            photometry = list(csv.DictReader(handle))
        self.assertEqual(len(photometry), 45)
        self.assertEqual({row["status"] for row in photometry}, {"detection"})
        self.assertEqual({row["fixed_ra"] for row in photometry}, {"13:05:27.4798"})
        self.assertEqual({row["fixed_dec"] for row in photometry}, {"-49:28:05.404"})
        self.assertFalse(any("/home/" in row["product"] for row in photometry))
        self.assertEqual(
            {row["component_identity"] for row in photometry},
            {"Emig2020_source18_AGN_maser_nucleus"},
        )
        self.assertEqual(
            {row["row_disposition"] for row in photometry},
            {"retained_source18_replacement"},
        )
        self.assertEqual({row["noise_guard_status"] for row in photometry}, {"passed"})
        self.assertEqual(
            {row["adopted_noise_estimator"] for row in photometry}, {"scaled_mad"}
        )
        self.assertTrue(
            all(
                math.isclose(
                    float(row["adopted_noise_mjy_per_beam"]),
                    float(row["noise_scaled_mad_mjy_per_beam"]),
                    rel_tol=0.0,
                    abs_tol=1e-15,
                )
                for row in photometry
            )
        )
        self.assertGreaterEqual(
            min(float(row["effective_independent_samples"]) for row in photometry),
            50.0,
        )
        self.assertEqual(
            sum(
                "aperture_truncated_by_image" in row["noise_review_flags"]
                for row in photometry
            ),
            6,
        )
        self.assertEqual(
            sum(
                row["noise_review_disposition"]
                == "diagnostic_flags_recorded_not_classification_overrides"
                for row in photometry
            ),
            43,
        )

        anchor = [
            row for row in photometry if math.isclose(float(row["frequency_ghz"]), 349.4)
        ]
        self.assertEqual(len(anchor), 1)
        self.assertEqual(anchor[0]["old_flux_field"], "15.84")
        self.assertEqual(anchor[0]["new_flux_field"], "13.08")
        self.assertEqual(anchor[0]["new_snr"], "7.051")
        self.assertTrue(
            math.isclose(
                float(anchor[0]["adopted_noise_mjy_per_beam"]),
                1.8551760085006326,
                rel_tol=0.0,
                abs_tol=1e-15,
            )
        )
        self.assertTrue(
            math.isclose(
                float(anchor[0]["noise_rms_mjy_per_beam"]),
                3.110572446005386,
                rel_tol=0.0,
                abs_tol=1e-15,
            )
        )

        with (PACKAGE_DIR / "ngc4945_component_audit.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            components = {row["record_id"]: row for row in csv.DictReader(handle)}
        rejected = components["rejected-source17-archived-anchor"]
        accepted = components["accepted-source18-fixed-canonical"]
        self.assertEqual(rejected["decision"], "rejected")
        self.assertEqual(rejected["included_in_canonical_table"], "no")
        self.assertEqual(accepted["decision"], "accepted")
        self.assertEqual(accepted["included_in_canonical_table"], "yes")
        self.assertEqual(accepted["component_identity"], "Emig2020_source18_AGN_maser_nucleus")

        result = subprocess.run(
            [sys.executable, str(PACKAGE_DIR / "rerun_ngc4945_nuclear.py"), "--check"],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
