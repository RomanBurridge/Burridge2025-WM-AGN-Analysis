from __future__ import annotations

import ast
import csv
import hashlib
import math
import tempfile
import unittest
from pathlib import Path
from unittest import mock


PACKAGE_DIR = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(PACKAGE_DIR))

from continuum_noise_inventory import (  # noqa: E402
    EXPECTED_FORCED_FIT_POLICY_BY_PRODUCT,
    FORCED_FIT_POLICY_CONVERGED,
    FORCED_FIT_POLICY_NONCONVERGED,
    GENERIC_FINAL_NOISE_AUDIT_FIELDS,
    validate_final_noise_inventory,
)
from rerun_ngc4945_nuclear import PHOTOMETRY_FIELDS  # noqa: E402
from continuum_product_additions import (  # noqa: E402
    ADDITION_POLICY_FIELDS,
    ADDITION_POLICY_VERSION,
    ADMISSION_DECISION,
)
from continuum_boundary_rejections import (  # noqa: E402
    BOUNDARY_REJECTION_EVIDENCE_FIELDS,
    BOUNDARY_REJECTION_EVIDENCE_NAME,
    BOUNDARY_REJECTION_EVIDENCE_VERSION,
    BOUNDARY_REJECTION_POLICY_FIELDS,
    BOUNDARY_REJECTION_POLICY_VERSION,
    EXPECTED_CLASSIFICATION,
    EXPECTED_ROW_DISPOSITION,
    REJECTION_DECISION,
    REJECTION_GATE,
)


def write_csv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class FinalNoiseInventoryTests(unittest.TestCase):
    BOUNDARY_SPECS = (
        ("BoundaryA", "VLA/BoundaryA/a.fits", "1.50", "1985/01/01", "91.0"),
        ("BoundaryB", "VLA/BoundaryB/b.fits", "5.00", "1986/02/02", "95.0"),
        ("BoundaryC", "VLA/BoundaryC/c.fits", "15.00", "1987/03/03", "93.0"),
    )

    def test_release_forced_fit_policy_map_pins_exact_ten_two_split(self) -> None:
        self.assertEqual(len(EXPECTED_FORCED_FIT_POLICY_BY_PRODUCT), 12)
        nonconverged = {
            product
            for product, policy in EXPECTED_FORCED_FIT_POLICY_BY_PRODUCT.items()
            if policy == FORCED_FIT_POLICY_NONCONVERGED
        }
        self.assertEqual(
            nonconverged,
            {
                "analysis/sma_calibration/230316_128/images/nitercal_0/"
                "flag2_cal1/fields/science/NGC6264_test.image",
                "analysis/sma_calibration/230316_128/images/nitercal_0/"
                "flag2_cal1/fields/science/UGC6093_test.image",
            },
        )
        self.assertEqual(
            sum(
                policy == FORCED_FIT_POLICY_CONVERGED
                for policy in EXPECTED_FORCED_FIT_POLICY_BY_PRODUCT.values()
            ),
            10,
        )

    def make_inventory(
        self,
        directory: Path,
    ) -> tuple[Path, Path, Path, list[dict[str, str]], list[dict[str, str]]]:
        table_path = directory / "publication.txt"
        generic_path = directory / "generic.csv"
        source18_path = directory / "source18.csv"
        table_rows = []
        generic = []

        def ordinary(source: str, product: str) -> None:
            table_rows.append([source, "ALMA", "230.00", "4.00", "2020/01/01_1", "1", "4.000", product])
            record = {field: "" for field in GENERIC_FINAL_NOISE_AUDIT_FIELDS}
            beam_area = math.pi / (4.0 * math.log(2.0))
            beam_areas = 95.0 / beam_area
            mask_semiaxis = 3.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))
            record.update(
                {
                    "schema_version": "1",
                    "run_context": "test",
                    "source": source,
                    "telescope": "ALMA",
                    "frequency_ghz": "230",
                    "product": product,
                    "forced_policy": "no",
                    "measurement_scope": "publication_final_10bmaj_source_excluded",
                    "fit_selection_policy": (
                        "single_initial_fixed_beam_fit_no_widened_retry"
                    ),
                    "aperture_center_policy": "fitted_component_pixelcoords",
                    "aperture_center_x_pixel": "49.5",
                    "aperture_center_y_pixel": "49.5",
                    "fit_center_reference_frame": "ICRS",
                    "fit_center_ra_deg": "195.0",
                    "fit_center_dec_deg": "-49.0",
                    "aperture_width_arcsec": "10",
                    # Pixel-center rasterization can select slightly more area
                    # than the requested continuous square; this is a ratio,
                    # not a probability constrained to one.
                    "aperture_pre_exclusion_npts": "101",
                    "aperture_selected_pixel_area_over_requested_area": "1.01",
                    "mask_semantics": "crtf_ellipse_axes_are_semiaxes",
                    "mask_major_semiaxis_arcsec": repr(mask_semiaxis),
                    "mask_minor_semiaxis_arcsec": repr(mask_semiaxis),
                    "mask_pa_deg": "0",
                    "beam_major_arcsec": "1",
                    "beam_minor_arcsec": "1",
                    "beam_pa_deg": "0",
                    "fit_beam_major_delta_arcsec": "0",
                    "fit_beam_minor_delta_arcsec": "0",
                    "fit_beam_pa_delta_deg": "0",
                    "fit_beam_pa_guard_applicable": "no",
                    "image_width_pixels": "100",
                    "image_height_pixels": "100",
                    "aperture_half_width_x_pixels": "5",
                    "aperture_half_width_y_pixels": "5",
                    "aperture_left_clearance_pixels": "45",
                    "aperture_right_clearance_pixels": "45",
                    "aperture_bottom_clearance_pixels": "45",
                    "aperture_top_clearance_pixels": "45",
                    "aperture_minimum_edge_clearance_pixels": "45",
                    "aperture_edge_status": "fully_contained",
                    "npts": "95",
                    "pixel_width_arcsec": "1",
                    "pixel_height_arcsec": "1",
                    "pixel_area_arcsec2": "1",
                    "valid_unmasked_area_arcsec2": "95",
                    "gaussian_beam_area_arcsec2": repr(beam_area),
                    "independent_beam_areas": repr(beam_areas),
                    "effective_independent_samples": repr(beam_areas),
                    "minimum_independent_beams": "50",
                    "brightness_unit": "Jy/beam",
                    "rms_mjy_per_beam": "1",
                    "sigma_mjy_per_beam": "1",
                    "mean_mjy_per_beam": "0",
                    "scaled_mad_mjy_per_beam": "1",
                    "adopted_noise_estimator": "scaled_mad",
                    "adopted_noise_mjy_per_beam": "1",
                    "mean_over_sigma": "0",
                    "rms_over_scaled_mad": "1",
                    "review_flags": "none",
                    "review_disposition": "no_diagnostic_flags",
                    "flux_mjy": "4",
                    "snr": "4",
                    "classification": "detection",
                    "ordinary_upper_limit_mjy": "4",
                    "reported_limit_mjy": "",
                    "reported_flux_field": "4.00",
                    "row_disposition": (
                        "rejected_source17_requires_source18_join"
                        if source == "NGC4945"
                        else "generic_table_retained"
                    ),
                    "guard_status": "passed",
                }
            )
            generic.append(record)

        for index in range(45):
            ordinary("NGC4945", f"ALMA/NGC4945/product-{index:03d}.fits")
        for index in range(349):
            ordinary("TEST", f"ALMA/TEST/product-{index:03d}.fits")
        for product, fit_policy in EXPECTED_FORCED_FIT_POLICY_BY_PRODUCT.items():
            source = Path(product).name.removesuffix("_test.image")
            if source == "NGC5765b":
                source = "NGC5765B"
            table_rows.append(
                [source, "SMA", "230.00", "3.00*", "2020/01/01_1", "1", "NA", product]
            )
            record = {field: "" for field in GENERIC_FINAL_NOISE_AUDIT_FIELDS}
            forced_beam_area = math.pi / (4.0 * math.log(2.0))
            forced_beam_areas = 1.0 / forced_beam_area
            record.update(
                {
                    "schema_version": "1",
                    "run_context": "test",
                    "source": source,
                    "telescope": "SMA",
                    "frequency_ghz": "230",
                    "product": product,
                    "forced_policy": "yes_dirty_sma_nitercal_0",
                    "measurement_scope": "reviewed_forced_dirty_sma_one_beam_region",
                    "fit_selection_policy": fit_policy,
                    "aperture_center_policy": "catalog_or_anchor_coordinate",
                    "aperture_width_arcsec": "one_restoring_beam_rotbox",
                    "mask_semantics": (
                        "none_source_not_excluded_reviewed_forced_policy"
                    ),
                    "beam_major_arcsec": "1",
                    "beam_minor_arcsec": "1",
                    "beam_pa_deg": "0",
                    "npts": "1",
                    "pixel_width_arcsec": "1",
                    "pixel_height_arcsec": "1",
                    "pixel_area_arcsec2": "1",
                    "valid_unmasked_area_arcsec2": "1",
                    "gaussian_beam_area_arcsec2": repr(forced_beam_area),
                    "independent_beam_areas": repr(forced_beam_areas),
                    "effective_independent_samples": repr(forced_beam_areas),
                    "minimum_independent_beams": "not_applicable",
                    "brightness_unit": "Jy/beam",
                    "rms_mjy_per_beam": "1",
                    "sigma_mjy_per_beam": "1",
                    "mean_mjy_per_beam": "0",
                    "scaled_mad_mjy_per_beam": "1",
                    "flux_mjy": (
                        "0.5" if fit_policy == FORCED_FIT_POLICY_CONVERGED else ""
                    ),
                    "classification": "forced_upper_limit",
                    "adopted_noise_estimator": "classic_rms_reviewed_forced_policy",
                    "adopted_noise_mjy_per_beam": "1",
                    "reported_limit_mjy": "3",
                    "reported_flux_field": "3.00*",
                    "review_flags": "reviewed_forced_policy",
                    "review_disposition": "preapproved_forced_policy",
                    "row_disposition": "generic_table_retained",
                    "guard_status": "not_applicable_reviewed_forced_policy",
                }
            )
            generic.append(record)

        source18 = []
        for record in generic:
            if record["source"] != "NGC4945":
                continue
            row = {field: "" for field in PHOTOMETRY_FIELDS}
            row.update(
                {
                    "source": "NGC4945",
                    "telescope": "ALMA",
                    "frequency_ghz": "230.00",
                    "date_label": "2020/01/01_1",
                    "product": record["product"],
                    "new_flux_field": "4.00",
                    "new_snr": "4.000",
                    "component_identity": "Emig2020_source18_AGN_maser_nucleus",
                    "rejected_seed_identity": "Emig2020_source17_SSC_SNR_associated",
                    "row_disposition": "retained_source18_replacement",
                    "noise_guard_status": "passed",
                }
            )
            source18.append(row)

        table_path.write_text(
            "".join(f"{row!r}\n" for row in table_rows),
            encoding="utf-8",
        )
        generic.sort(key=lambda row: row["product"])
        write_csv(generic_path, GENERIC_FINAL_NOISE_AUDIT_FIELDS, generic)
        write_csv(source18_path, PHOTOMETRY_FIELDS, source18)
        return table_path, generic_path, source18_path, generic, source18

    def validate(
        self,
        table_path: Path,
        generic_path: Path,
        source18_path: Path,
        boundary_policy_path: Path | None = None,
    ) -> dict[str, int]:
        with mock.patch(
            "rerun_ngc4945_nuclear.validate_photometry_records",
            autospec=True,
        ):
            return validate_final_noise_inventory(
                table_path,
                generic_path,
                source18_path,
                boundary_rejection_policy_path=boundary_policy_path,
            )

    def add_boundary_contract(
        self,
        directory: Path,
        generic_path: Path,
        generic: list[dict[str, str]],
    ) -> Path:
        evidence_path = directory / BOUNDARY_REJECTION_EVIDENCE_NAME
        evidence_rows: list[dict[str, str]] = []
        for index, (_, product, _, _, percent) in enumerate(
            self.BOUNDARY_SPECS, start=1
        ):
            model_hash = f"{index}" * 64
            evidence_rows.append(
                {
                    "schema_version": "1",
                    "evidence_version": BOUNDARY_REJECTION_EVIDENCE_VERSION,
                    "rejection_id": f"CBR-{index:03d}",
                    "product": product,
                    "independent_source_file": "synthetic/excluded_results.csv",
                    "independent_source_sha256": "a" * 64,
                    "independent_source_record_index": str(index),
                    "fit_converged": "true",
                    "boundary_result": "true",
                    "rerun_outcome": "on_boundary",
                    "historical_boundary_diagnostic": "true",
                    "boundary_highest_percent": percent,
                    "archived_final_model_exists": "true",
                    "archived_final_model_flux_mjy": str(index),
                    "rerun_model_flux_mjy": str(index),
                    "model_flux_delta_mjy": "0.0",
                    "retry3_model_payload_sha256": model_hash,
                    "authoritative7_model_payload_sha256": model_hash,
                    "model_payload_size_bytes": str(100 + index),
                    "input_metadata_unchanged": "true",
                    "status": "ok",
                }
            )
        write_csv(
            evidence_path,
            BOUNDARY_REJECTION_EVIDENCE_FIELDS,
            evidence_rows,
        )
        evidence_sha = hashlib.sha256(evidence_path.read_bytes()).hexdigest()

        policy_rows: list[dict[str, str]] = []
        ordinary_template = next(
            row for row in generic if row["source"] == "TEST"
        )
        for index, (source, product, frequency, date, percent) in enumerate(
            self.BOUNDARY_SPECS, start=1
        ):
            model_hash = f"{index}" * 64
            policy_rows.append(
                {
                    "schema_version": "1",
                    "policy_version": BOUNDARY_REJECTION_POLICY_VERSION,
                    "rejection_id": f"CBR-{index:03d}",
                    "source": source,
                    "telescope": "VLA",
                    "expected_frequency_ghz": frequency,
                    "expected_observation_date": date,
                    "product": product,
                    "product_sha256": f"{index + 3}" * 64,
                    "size_bytes": str(200 + index),
                    "independent_evidence_record": (
                        f"{BOUNDARY_REJECTION_EVIDENCE_NAME}#CBR-{index:03d}"
                    ),
                    "independent_evidence_sha256": evidence_sha,
                    "historical_boundary_highest_percent": percent,
                    "historical_model_payload_sha256": model_hash,
                    "expected_row_disposition": EXPECTED_ROW_DISPOSITION,
                    "expected_classification": EXPECTED_CLASSIFICATION,
                    "rejection_gate": REJECTION_GATE,
                    "rejection_decision": REJECTION_DECISION,
                    "reviewer": "noise-inventory-test",
                    "review_date": "2026-08-20",
                    "historical_evidence": "synthetic exact boundary reproduction",
                    "rationale": "Synthetic reviewed boundary rejection.",
                }
            )
            record = dict(ordinary_template)
            record.update(
                {
                    "source": source,
                    "telescope": "VLA",
                    "frequency_ghz": frequency,
                    "product": product,
                    "row_disposition": EXPECTED_ROW_DISPOSITION,
                    "classification": EXPECTED_CLASSIFICATION,
                }
            )
            generic.append(record)
        policy_path = directory / "boundary.csv"
        write_csv(policy_path, BOUNDARY_REJECTION_POLICY_FIELDS, policy_rows)
        generic.sort(key=lambda row: row["product"])
        write_csv(generic_path, GENERIC_FINAL_NOISE_AUDIT_FIELDS, generic)
        return policy_path

    def test_exact_generic_and_source18_union_passes(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            paths = self.make_inventory(Path(raw))
            summary = self.validate(*paths[:3])
        self.assertEqual(
            summary,
            {
                "publication_rows": 406,
                "admitted_input_additions": 0,
                "publication_additions": 0,
                "input_only_additions": 0,
                "generic_audit_rows": 406,
                "generic_publication_product_rows": 406,
                "ordinary_generic_rows": 394,
                "ordinary_publication_product_rows": 394,
                "forced_dirty_sma_rows": 12,
                "boundary_rejection_policy_rows": 0,
                "boundary_rejected_noise_audit_rows": 0,
                "generic_ngc4945_rejected_rows": 45,
                "source18_retained_rows": 45,
                "retained_union_rows": 406,
            },
        )

    def test_exact_boundary_union_keeps_three_diagnostics_out_of_publication(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            table, generic_path, source18_path, generic, _ = self.make_inventory(
                directory
            )
            boundary_policy = self.add_boundary_contract(
                directory, generic_path, generic
            )
            summary = self.validate(
                table, generic_path, source18_path, boundary_policy
            )
        self.assertEqual(summary["publication_rows"], 406)
        self.assertEqual(summary["generic_audit_rows"], 409)
        self.assertEqual(summary["generic_publication_product_rows"], 406)
        self.assertEqual(summary["ordinary_generic_rows"], 397)
        self.assertEqual(summary["ordinary_publication_product_rows"], 394)
        self.assertEqual(summary["boundary_rejection_policy_rows"], 3)
        self.assertEqual(summary["boundary_rejected_noise_audit_rows"], 3)

    def test_boundary_inventory_and_disposition_mutations_fail_closed(self) -> None:
        cases = (
            ("missing audit row", "product inventories differ"),
            ("arbitrary extra row", "product inventories differ"),
            ("wrong disposition", "identity/disposition is stale"),
            ("published boundary row", "must be absent from the publication table"),
        )
        for mutation, message in cases:
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as raw:
                directory = Path(raw)
                table, generic_path, source18_path, generic, _ = self.make_inventory(
                    directory
                )
                boundary_policy = self.add_boundary_contract(
                    directory, generic_path, generic
                )
                boundary_products = {spec[1] for spec in self.BOUNDARY_SPECS}
                if mutation == "missing audit row":
                    generic[:] = [
                        row
                        for row in generic
                        if row["product"] != self.BOUNDARY_SPECS[0][1]
                    ]
                elif mutation == "arbitrary extra row":
                    extra = dict(
                        next(row for row in generic if row["source"] == "TEST")
                    )
                    extra["product"] = "VLA/Unreviewed/extra.fits"
                    generic.append(extra)
                elif mutation == "wrong disposition":
                    boundary_row = next(
                        row for row in generic if row["product"] in boundary_products
                    )
                    boundary_row["row_disposition"] = "generic_table_retained"
                else:
                    rows = [
                        ast.literal_eval(line)
                        for line in table.read_text(encoding="utf-8").splitlines()
                        if line.strip()
                    ]
                    rows[45][7] = self.BOUNDARY_SPECS[0][1]
                    table.write_text(
                        "".join(f"{row!r}\n" for row in rows), encoding="utf-8"
                    )
                generic.sort(key=lambda row: row["product"])
                write_csv(
                    generic_path, GENERIC_FINAL_NOISE_AUDIT_FIELDS, generic
                )
                with self.assertRaisesRegex(RuntimeError, message):
                    self.validate(
                        table, generic_path, source18_path, boundary_policy
                    )

    def test_nonconverged_forced_fit_has_no_fabricated_flux(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            paths = self.make_inventory(Path(raw))
            generic = paths[3]
            forced = next(
                row
                for row in generic
                if row["fit_selection_policy"] == FORCED_FIT_POLICY_NONCONVERGED
            )
            write_csv(paths[1], GENERIC_FINAL_NOISE_AUDIT_FIELDS, generic)
            summary = self.validate(*paths[:3])
        self.assertEqual(summary["forced_dirty_sma_rows"], 12)

    def test_nonconverged_forced_fit_rejects_a_fabricated_flux(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            paths = self.make_inventory(Path(raw))
            generic = paths[3]
            forced = next(
                row
                for row in generic
                if row["fit_selection_policy"] == FORCED_FIT_POLICY_NONCONVERGED
            )
            forced["flux_mjy"] = "0.5"
            write_csv(paths[1], GENERIC_FINAL_NOISE_AUDIT_FIELDS, generic)
            with self.assertRaisesRegex(RuntimeError, "carries a fitted flux"):
                self.validate(*paths[:3])

    def test_forced_fit_policy_and_flux_pairing_is_bidirectional(self) -> None:
        cases = (
            (FORCED_FIT_POLICY_CONVERGED, "", "invalid flux_mjy"),
            ("unknown_forced_fit_policy", "", "fit-selection policy is invalid"),
        )
        for policy, flux, message in cases:
            with self.subTest(policy=policy), tempfile.TemporaryDirectory() as raw:
                paths = self.make_inventory(Path(raw))
                generic = paths[3]
                forced = next(
                    row for row in generic if row["forced_policy"] != "no"
                )
                forced["fit_selection_policy"] = policy
                forced["flux_mjy"] = flux
                write_csv(paths[1], GENERIC_FINAL_NOISE_AUDIT_FIELDS, generic)
                with self.assertRaisesRegex(RuntimeError, message):
                    self.validate(*paths[:3])

    def test_forced_convergence_inventory_is_the_exact_reviewed_ten_two_split(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            paths = self.make_inventory(Path(raw))
            generic = paths[3]
            nonconverged = next(
                row
                for row in generic
                if row["fit_selection_policy"] == FORCED_FIT_POLICY_NONCONVERGED
            )
            nonconverged["fit_selection_policy"] = FORCED_FIT_POLICY_CONVERGED
            nonconverged["flux_mjy"] = "0.5"
            write_csv(paths[1], GENERIC_FINAL_NOISE_AUDIT_FIELDS, generic)
            with self.assertRaisesRegex(RuntimeError, "10-converged/2-nonconverged"):
                self.validate(*paths[:3])

    def test_publication_and_input_only_addition_counts_are_dynamic(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            table, generic_path, source18_path, generic, _ = self.make_inventory(directory)
            product = "ALMA/New/new.fits"
            table_rows = [
                ast.literal_eval(line)
                for line in table.read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
            table_rows.append(
                ["New", "ALMA", "500.00", "4.00", "2021/02/03_1", "1", "4.000", product]
            )
            table.write_text("".join(f"{row!r}\n" for row in table_rows), encoding="utf-8")
            record = dict(next(row for row in generic if row["forced_policy"] == "no"))
            record.update(
                {
                    "source": "New",
                    "frequency_ghz": "500.00",
                    "product": product,
                    "reported_flux_field": "4.00",
                    "row_disposition": "generic_table_retained",
                }
            )
            generic.append(record)
            generic.sort(key=lambda row: row["product"])
            write_csv(generic_path, GENERIC_FINAL_NOISE_AUDIT_FIELDS, generic)
            policy = directory / "additions.csv"
            rows = []
            for index, (source, policy_product, disposition) in enumerate(
                (
                    ("New", product, "publication_row_required"),
                    ("NoRow", "ALMA/NoRow/no.fits", "input_only_expected_no_publication_row"),
                )
            ):
                rows.append(
                    {
                        "schema_version": "1",
                        "policy_version": ADDITION_POLICY_VERSION,
                        "admission_id": f"ADD-{index + 1:03d}",
                        "source": source,
                        "telescope": "ALMA",
                        "expected_frequency_ghz": f"{500 + index * 10:.2f}",
                        "expected_date_label": f"202{index + 1}/02/03_1",
                        "product": policy_product,
                        "product_sha256": f"{index + 1}" * 64,
                        "expected_run_disposition": disposition,
                        "admission_decision": ADMISSION_DECISION,
                        "reviewer": "test",
                        "review_date": "2026-08-20",
                        "rationale": "Synthetic addition.",
                    }
                )
            rows.sort(key=lambda row: row["product"])
            write_csv(policy, ADDITION_POLICY_FIELDS, rows)
            with mock.patch(
                "rerun_ngc4945_nuclear.validate_photometry_records",
                autospec=True,
            ):
                summary = validate_final_noise_inventory(
                    table, generic_path, source18_path, policy
                )
        self.assertEqual(summary["publication_rows"], 407)
        self.assertEqual(summary["ordinary_generic_rows"], 395)
        self.assertEqual(summary["admitted_input_additions"], 2)
        self.assertEqual(summary["publication_additions"], 1)
        self.assertEqual(summary["input_only_additions"], 1)

    def test_missing_extra_duplicate_and_nonterminal_generic_rows_fail(self) -> None:
        cases = (
            ("missing", "product inventories differ"),
            ("duplicate", "duplicate products"),
            ("pending", "nonterminal or incorrect"),
            ("failed_guard", "guard did not pass"),
            ("zero_mad", "invalid scaled_mad"),
            ("neff_math", "effective independent samples"),
            ("edge_status", "aperture edge status"),
            ("fixed_delta", "fit-beam major delta"),
            ("upper_limit", "upper-limit invariant"),
            ("review_flags", "review flags"),
            ("estimator", "adopted-noise estimator"),
            ("fit_policy", "aperture/mask policy"),
            ("fit_frame", "fit-center reference frame"),
            ("fit_ra", "fit-center RA is out of range"),
            ("fit_dec", "fit-center Dec is out of range"),
        )
        for mutation, message in cases:
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as raw:
                directory = Path(raw)
                table, generic_path, source18_path, generic, _ = self.make_inventory(directory)
                if mutation == "missing":
                    generic.pop()
                elif mutation == "duplicate":
                    generic[-1] = dict(generic[-2])
                elif mutation == "pending":
                    generic[0]["row_disposition"] = "candidate_pending"
                elif mutation == "failed_guard":
                    ordinary = next(row for row in generic if row["forced_policy"] == "no")
                    ordinary["guard_status"] = "failed_insufficient_independent_samples"
                elif mutation == "zero_mad":
                    ordinary = next(row for row in generic if row["forced_policy"] == "no")
                    ordinary["scaled_mad_mjy_per_beam"] = "0"
                else:
                    ordinary = next(row for row in generic if row["forced_policy"] == "no")
                    field, value = {
                        "neff_math": ("effective_independent_samples", "50"),
                        "edge_status": (
                            "aperture_edge_status",
                            "truncated_available_valid_intersection",
                        ),
                        "fixed_delta": ("fit_beam_major_delta_arcsec", "0.1"),
                        "upper_limit": ("ordinary_upper_limit_mjy", "5"),
                        "review_flags": ("review_flags", "nonzero_background_mean"),
                        "estimator": ("adopted_noise_estimator", "classic_rms"),
                        "fit_policy": ("fit_selection_policy", "widened_retry"),
                        "fit_frame": ("fit_center_reference_frame", "GALACTIC"),
                        "fit_ra": ("fit_center_ra_deg", "360"),
                        "fit_dec": ("fit_center_dec_deg", "-91"),
                    }[mutation]
                    ordinary[field] = value
                generic.sort(key=lambda row: row["product"])
                write_csv(generic_path, GENERIC_FINAL_NOISE_AUDIT_FIELDS, generic)
                with self.assertRaisesRegex(RuntimeError, message):
                    self.validate(table, generic_path, source18_path)

    def test_source18_missing_or_component_identity_collision_fails(self) -> None:
        cases = (
            ("missing", "does not contain 45 unique"),
            ("collision", "component identity/disposition is invalid"),
        )
        for mutation, message in cases:
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as raw:
                directory = Path(raw)
                table, generic, source18_path, _, source18 = self.make_inventory(directory)
                if mutation == "missing":
                    source18.pop()
                else:
                    source18[0]["component_identity"] = source18[0]["rejected_seed_identity"]
                write_csv(source18_path, PHOTOMETRY_FIELDS, source18)
                with self.assertRaisesRegex(RuntimeError, message):
                    self.validate(table, generic, source18_path)

    def test_forced_dirty_sma_audit_math_and_not_applicable_fields_are_strict(self) -> None:
        cases = (
            ("effective_independent_samples", "1", "forced effective independent samples"),
            ("fit_beam_major_delta_arcsec", "0", "not-applicable fields"),
            ("scaled_mad_mjy_per_beam", "-1", "scaled MAD is negative"),
        )
        for field, value, message in cases:
            with self.subTest(field=field), tempfile.TemporaryDirectory() as raw:
                directory = Path(raw)
                table, generic_path, source18_path, generic, _ = self.make_inventory(
                    directory
                )
                forced = next(row for row in generic if row["forced_policy"] != "no")
                forced[field] = value
                generic.sort(key=lambda row: row["product"])
                write_csv(generic_path, GENERIC_FINAL_NOISE_AUDIT_FIELDS, generic)
                with self.assertRaisesRegex(RuntimeError, message):
                    self.validate(table, generic_path, source18_path)

    def test_rejected_generic_may_differ_but_source18_must_match_table(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            directory = Path(raw)
            table, generic_path, source18_path, generic, source18 = self.make_inventory(
                directory
            )
            rejected = next(row for row in generic if row["source"] == "NGC4945")
            rejected.update(
                {
                    "flux_mjy": "5",
                    "snr": "5",
                    "ordinary_upper_limit_mjy": "5",
                    "reported_flux_field": "5.00",
                }
            )
            generic.sort(key=lambda row: row["product"])
            write_csv(generic_path, GENERIC_FINAL_NOISE_AUDIT_FIELDS, generic)
            self.validate(table, generic_path, source18_path)

            source18[0]["new_flux_field"] = "5.00"
            write_csv(source18_path, PHOTOMETRY_FIELDS, source18)
            with self.assertRaisesRegex(RuntimeError, "disagrees with publication row"):
                self.validate(table, generic_path, source18_path)


if __name__ == "__main__":
    unittest.main()
