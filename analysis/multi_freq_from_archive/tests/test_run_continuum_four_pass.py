import argparse
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


PACKAGE_DIR = Path(__file__).resolve().parents[1]
if str(PACKAGE_DIR) not in sys.path:
    sys.path.insert(0, str(PACKAGE_DIR))

from run_continuum_four_pass import (  # noqa: E402
    CORRECTION_FIELDS,
    HISTORICAL_CANONICAL_RELATIVE,
    LEGACY_CORRECTION_FIELDS,
    MARKER_NAME,
    PHOTOMETRY_FIELDS,
    REVIEW_DELTA_FIELDS,
    REVIEW_ADDITION_FIELDS,
    REVIEW_FREEZE_ARTIFACT_NAMES,
    REVIEW_SCOPE,
    REVIEW_STATUS,
    STAGES,
    WORKFLOW_INPUT_NAMES,
    _artifact_snapshot,
    adjudicate_workflow,
    analysis_utils_bundle_inventory,
    assert_clean_outputs,
    assert_isolated_layout,
    compare_current_corrected_table,
    compare_reviewed_current_to_generated,
    freeze_analysis_utils_bundle,
    initialize_workspace,
    load_review_freeze,
    load_reviewed_exclusion_products,
    preflight_python,
    prepare_workspace,
    sha256_file,
    stage_commands,
    validate_inputs,
    validate_correction_manifest_transition,
    validate_prepared_archive_inputs,
    workflow_input_hashes,
)
from continuum_classification_identity import POLICY_FIELDS  # noqa: E402
from continuum_boundary_rejections import (  # noqa: E402
    BOUNDARY_REJECTION_EVIDENCE_FIELDS,
    BOUNDARY_REJECTION_EVIDENCE_VERSION,
    BOUNDARY_REJECTION_POLICY_FIELDS,
    BOUNDARY_REJECTION_POLICY_VERSION,
)
from continuum_product_additions import ADDITION_POLICY_FIELDS  # noqa: E402
try:
    from .frozen_input_fixture import fixture_load, fixture_select, write_fixture_manifest
except ImportError:  # direct unittest discovery adds tests/ rather than package root
    from frozen_input_fixture import fixture_load, fixture_select, write_fixture_manifest


class FourPassOrchestratorTests(unittest.TestCase):
    def make_workspace(self, root: Path) -> tuple[Path, Path, Path]:
        analysis = root / "analysis/multi_freq_from_archive"
        analysis.mkdir(parents=True)
        (analysis / "get_continuum_params.py").write_text("print('dummy')\n")
        (analysis / "rerun_ngc4945_nuclear.py").write_text("print('dummy')\n")
        (analysis / "ALMA").mkdir()
        (analysis / "VLA").mkdir()
        sma = root / "analysis/sma_calibration/run/images/science"
        sma.mkdir(parents=True)
        selected = root / "analysis/sma_calibration/selected.csv"
        selected_image = root / "analysis/sma_calibration/clean/Test_test.image"
        selected_residual = root / "analysis/sma_calibration/clean/Test_test.residual"
        selected_image.mkdir(parents=True)
        selected_residual.mkdir()
        selected.write_text(
            "date,flagnum,source,image_path,residual_path\n"
            "230316,2,Test,analysis/sma_calibration/clean/Test_test.image,"
            "analysis/sma_calibration/clean/Test_test.residual\n"
        )
        return analysis, sma, selected

    @staticmethod
    def write_sorted_rows(path: Path, rows: list[list[object]]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("".join(f"{row!r}\n" for row in rows), encoding="utf-8")

    @staticmethod
    def write_csv_rows(
        path: Path,
        fields: tuple[str, ...] | list[str],
        rows: list[dict[str, object]],
    ) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    @staticmethod
    def _write_minimal_fits(path: Path, frequency_ghz: float, date: str) -> None:
        def card(keyword: str, value: str) -> bytes:
            return f"{keyword:<8}= {value:<70}".encode("ascii")[:80].ljust(80)

        cards = [
            card("SIMPLE", "T"),
            card("DATE-OBS", f"'{date}T00:00:00.000'"),
            card("CTYPE3", "'FREQ'"),
            card("CRVAL3", repr(frequency_ghz * 1.0e9)),
            card("CUNIT3", "'Hz'"),
            b"END".ljust(80),
        ]
        payload = b"".join(cards).ljust(2880, b" ")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)

    def write_boundary_fixture(self, analysis: Path) -> list[dict[str, str]]:
        products = [
            ("CBR-001", "1.49", "1986-06-19", "boundary_001.fits"),
            ("CBR-002", "4.86", "1985-03-30", "boundary_002.fits"),
            ("CBR-003", "14.94", "1985-05-08", "boundary_003.fits"),
        ]
        evidence_rows: list[dict[str, str]] = []
        policy_inputs: list[tuple[str, str, str, str, str, int]] = []
        for index, (rejection_id, frequency, date, basename) in enumerate(
            products, start=1
        ):
            product = f"VLA/Boundary/{basename}"
            product_path = analysis / product
            self._write_minimal_fits(product_path, float(frequency), date)
            model_hash = f"{index}" * 64
            policy_inputs.append(
                (
                    rejection_id,
                    frequency,
                    date.replace("-", "/"),
                    product,
                    sha256_file(product_path),
                    product_path.stat().st_size,
                )
            )
            evidence_rows.append(
                {
                    "schema_version": "1",
                    "evidence_version": BOUNDARY_REJECTION_EVIDENCE_VERSION,
                    "rejection_id": rejection_id,
                    "product": product,
                    "independent_source_file": "synthetic/excluded_results.csv",
                    "independent_source_sha256": "a" * 64,
                    "independent_source_record_index": str(index),
                    "fit_converged": "true",
                    "boundary_result": "true",
                    "rerun_outcome": "on_boundary",
                    "historical_boundary_diagnostic": "true",
                    "boundary_highest_percent": str(90.0 + index),
                    "archived_final_model_exists": "true",
                    "archived_final_model_flux_mjy": str(float(index)),
                    "rerun_model_flux_mjy": str(float(index)),
                    "model_flux_delta_mjy": "0.0",
                    "retry3_model_payload_sha256": model_hash,
                    "authoritative7_model_payload_sha256": model_hash,
                    "model_payload_size_bytes": str(100 + index),
                    "input_metadata_unchanged": "true",
                    "status": "ok",
                }
            )
        evidence_path = analysis / "continuum_boundary_rejection_evidence.csv"
        self.write_csv_rows(
            evidence_path,
            BOUNDARY_REJECTION_EVIDENCE_FIELDS,
            evidence_rows,
        )
        evidence_hash = sha256_file(evidence_path)
        policy_rows: list[dict[str, str]] = []
        for index, values in enumerate(policy_inputs, start=1):
            rejection_id, frequency, date, product, product_hash, size = values
            policy_rows.append(
                {
                    "schema_version": "1",
                    "policy_version": BOUNDARY_REJECTION_POLICY_VERSION,
                    "rejection_id": rejection_id,
                    "source": "Boundary",
                    "telescope": "VLA",
                    "expected_frequency_ghz": frequency,
                    "expected_observation_date": date,
                    "product": product,
                    "product_sha256": product_hash,
                    "size_bytes": str(size),
                    "independent_evidence_record": (
                        f"continuum_boundary_rejection_evidence.csv#{rejection_id}"
                    ),
                    "independent_evidence_sha256": evidence_hash,
                    "historical_boundary_highest_percent": str(90.0 + index),
                    "historical_model_payload_sha256": f"{index}" * 64,
                    "expected_row_disposition": "rejected_on_boundary",
                    "expected_classification": "detection",
                    "rejection_gate": (
                        "half_maximum_model_contour_intersects_outer_ten_percent"
                    ),
                    "rejection_decision": (
                        "retain_noise_diagnostic_exclude_publication_row"
                    ),
                    "reviewer": "synthetic_test",
                    "review_date": "2026-08-20",
                    "historical_evidence": "synthetic_independent_reproduction",
                    "rationale": "Synthetic exact boundary contract fixture",
                }
            )
        self.write_csv_rows(
            analysis / "continuum_boundary_rejection_policy.csv",
            BOUNDARY_REJECTION_POLICY_FIELDS,
            policy_rows,
        )
        return policy_rows

    def make_reviewed_delta_fixture(
        self, root: Path
    ) -> tuple[Path, Path, Path, Path, list[list[object]], list[list[object]]]:
        canonical = root / "canonical.txt"
        generated = root / "generated.txt"
        policy = root / "exclusions.csv"
        ledger = root / "reviewed_current_to_generated_deltas.csv"
        canonical_rows = [
            [
                "Test",
                "ALMA",
                f"{index:.2f}",
                "1.00",
                "2020/01/01_1",
                "0.10000",
                "10.000",
                f"ALMA/Test/product_{index:03d}.fits\n",
            ]
            for index in range(409)
        ]
        generated_rows = [list(row) for row in canonical_rows[:406]]
        generated_rows[0][3] = "1.23"
        generated_rows[0][6] = "12.345"
        generated_rows[10][6] = "11.111"
        self.write_sorted_rows(canonical, canonical_rows)
        self.write_sorted_rows(generated, generated_rows)
        self.write_csv_rows(
            policy,
            ("product",),
            [
                {"product": f"ALMA/Test/product_{index:03d}.fits"}
                for index in range(406, 409)
            ],
        )
        self.write_csv_rows(
            ledger,
            REVIEW_DELTA_FIELDS,
            [
                {
                    "corrected_row_number": "1",
                    "source": "Test",
                    "telescope": "ALMA",
                    "frequency_ghz": "0.00",
                    "date_label": "2020/01/01_1",
                    "product": "ALMA/Test/product_000.fits",
                    "current_flux_field": "1.00",
                    "generated_flux_field": "1.23",
                    "current_snr": "10.000",
                    "generated_snr": "12.345",
                    "changed_fields": "flux_field+snr",
                    "review_reason": REVIEW_SCOPE,
                },
                {
                    "corrected_row_number": "11",
                    "source": "Test",
                    "telescope": "ALMA",
                    "frequency_ghz": "10.00",
                    "date_label": "2020/01/01_1",
                    "product": "ALMA/Test/product_010.fits",
                    "current_flux_field": "1.00",
                    "generated_flux_field": "1.00",
                    "current_snr": "10.000",
                    "generated_snr": "11.111",
                    "changed_fields": "snr",
                    "review_reason": REVIEW_SCOPE,
                },
            ],
        )
        return canonical, generated, policy, ledger, canonical_rows, generated_rows

    def make_adjudication_fixture(
        self, temporary: Path
    ) -> tuple[argparse.Namespace, dict[str, object], Path, Path, Path]:
        workspace = temporary / "release"
        root = workspace / "analysis/multi_freq_from_archive"
        evidence = root / "workflow_runs/local-rms-test"
        evidence.mkdir(parents=True)
        for name in WORKFLOW_INPUT_NAMES:
            path = root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"fixture for {name}\n", encoding="utf-8")
        boundary_rows = self.write_boundary_fixture(root)
        boundary_policy_path = root / "continuum_boundary_rejection_policy.csv"
        addition_policy_path = root / "continuum_product_addition_policy.csv"
        self.write_csv_rows(addition_policy_path, ADDITION_POLICY_FIELDS, [])

        fixture = temporary / "table-fixture"
        fixture.mkdir()
        canonical, generated, policy, ledger, _, _ = self.make_reviewed_delta_fixture(
            fixture
        )
        canonical_path = root / HISTORICAL_CANONICAL_RELATIVE
        generated_path = (
            root
            / "fitsumfiles/sorted/230316_2_fitsummary_fixedtochosen_sorted.txt"
        )
        policy_path = root / "continuum_product_exclusion_policy.csv"
        ledger_path = evidence / "reviewed_current_to_generated_deltas.csv"
        addition_ledger_path = evidence / "reviewed_new_product_additions.csv"
        self.write_csv_rows(addition_ledger_path, REVIEW_ADDITION_FIELDS, [])
        for source, destination in (
            (canonical, canonical_path),
            (generated, generated_path),
            (policy, policy_path),
            (ledger, ledger_path),
        ):
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(source.read_bytes())

        noise_path = (
            root
            / "fitsumfiles/noise_audits/230316_2_final_publication_noise_audit.csv"
        )
        self.write_csv_rows(
            noise_path,
            ("product", "row_disposition"),
            [
                {
                    "product": f"noise-{index:03d}",
                    "row_disposition": "generic_table_retained",
                }
                for index in range(406)
            ]
            + [
                {
                    "product": row["product"],
                    "row_disposition": "rejected_on_boundary",
                }
                for row in boundary_rows
            ],
        )
        photometry_path = root / "ngc4945_nuclear_photometry.csv"
        self.write_csv_rows(
            photometry_path,
            ("product",),
            [{"product": f"p{index}"} for index in range(45)],
        )
        component_path = root / "ngc4945_component_audit.csv"
        self.write_csv_rows(
            component_path,
            ("record_id",),
            [{"record_id": f"c{index}"} for index in range(7)],
        )
        corrections_path = root / "validated_continuum_corrections.csv"
        self.write_csv_rows(
            corrections_path,
            ("product",),
            [{"product": "ALMA/Test/product_000.fits"}],
        )
        preapplication_corrections_path = (
            evidence / "preapplication_validated_continuum_corrections.csv"
        )
        preapplication_corrections_path.write_bytes(corrections_path.read_bytes())
        inventory_validation_path = evidence / "final_noise_inventory_validation.json"
        inventory_validation = {
            "schema_version": 1,
            "status": "exact_validated_union",
            "publication_rows": 406,
            "admitted_input_additions": 0,
            "publication_additions": 0,
            "input_only_additions": 0,
            "generic_audit_rows": 409,
            "generic_publication_product_rows": 406,
            "ordinary_generic_rows": 397,
            "ordinary_publication_product_rows": 394,
            "forced_dirty_sma_rows": 12,
            "boundary_rejection_policy_rows": 3,
            "boundary_rejected_noise_audit_rows": 3,
            "generic_ngc4945_rejected_rows": 45,
            "source18_retained_rows": 45,
            "retained_union_rows": 406,
            "ngc4945_photometry_semantic_status": "exact_validated",
            "ngc4945_photometry_semantic_rows": 45,
            "ngc4945_component_semantic_status": "exact_validated",
            "ngc4945_component_semantic_rows": 7,
            "ngc4945_component_record_ids": [f"c{index}" for index in range(7)],
        }
        inventory_validation_path.write_text(
            json.dumps(inventory_validation, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        correction_transition_path = (
            evidence / "correction_manifest_transition_validation.json"
        )
        correction_transition = {
            "schema_version": 1,
            "status": "exact_source18_only_transition",
            "preapplication_rows": 66,
            "postrun_rows": 84,
            "unchanged_non_ngc4945_rows": 39,
            "added_ngc4945_rows": 18,
            "reviewed_ngc4945_rows": 45,
            "reviewed_ngc4945_products": [f"p{index}" for index in range(45)],
        }
        correction_transition_path.write_text(
            json.dumps(correction_transition, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        identity_policy_path = root / "continuum_classification_identity_policy.csv"
        identity_policy_row = {field: "fixture" for field in POLICY_FIELDS}
        identity_policy_row["schema_version"] = "1"
        self.write_csv_rows(
            identity_policy_path,
            POLICY_FIELDS,
            [identity_policy_row],
        )
        identity_validation_path = (
            evidence / "classification_identity_validation.json"
        )
        identity_validation = {
            "schema_version": 1,
            "status": "exact_required_identity_reviews_resolved",
            "required_review_rows": 1,
            "new_detection_rows": 1,
            "rms_mad_threshold_crossing_rows": 0,
            "new_product_admission_rows": 0,
            "admitted_product_upper_limit_rows": 0,
            "validated_reviews": [{"review_id": "fixture"}],
        }
        identity_validation_path.write_text(
            json.dumps(identity_validation, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        analysis_utils_execution_dir = evidence / "analysis_scripts"
        analysis_utils_execution_dir.mkdir()
        (analysis_utils_execution_dir / "analysisUtils.py").write_text(
            "# frozen analysisUtils fixture\n"
        )
        (analysis_utils_execution_dir / "mpfit.py").write_text(
            "# frozen mpfit fixture\n"
        )
        analysis_utils_inventory = analysis_utils_bundle_inventory(
            analysis_utils_execution_dir
        )
        analysis_utils_snapshot = evidence / "analysis_utils_snapshot.json"
        analysis_utils_snapshot.write_text(
            json.dumps(
                {
                    **analysis_utils_inventory,
                    "source_directory": str(
                        (Path(temporary) / "analysis-utils-source").resolve()
                    ),
                    "execution_directory": str(
                        analysis_utils_execution_dir.resolve()
                    ),
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )

        initialize_workspace(workspace)
        preparation_path = workspace / ".continuum-fit-inputs.json"
        preparation_path.write_text("{}\n", encoding="utf-8")
        selected_manifest = (
            workspace / "analysis/sma_calibration/selected_clean_iterations.csv"
        )
        selected_manifest.parent.mkdir(parents=True, exist_ok=True)
        selected_manifest.write_text(
            "date_code,flagnum,source,image_path,residual_path\n",
            encoding="utf-8",
        )
        inventory = {
            "archive_input_mode": "copy",
            "counts": {"ALMA": 133, "VLA": 307},
            "sma_discovery_inputs": 15,
            "total_discovery_inputs": 455,
            "addition_policy_rows": 0,
            "publication_addition_products": [],
            "boundary_rejection_policy_rows": 3,
            "boundary_rejection_products": [
                row["product"] for row in boundary_rows
            ],
            "file_metadata_sha256": "1" * 64,
            "file_content_sha256": "2" * 64,
            "file_metadata": [],
        }
        inventory_summary = {
            key: value for key, value in inventory.items() if key != "file_metadata"
        }
        integrity_path = evidence / "prepared_input_integrity.json"
        integrity_path.write_text(
            json.dumps(
                {
                    "status": "unchanged",
                    "before_file_metadata_sha256": "1" * 64,
                    "after_file_metadata_sha256": "1" * 64,
                    "archive_input_mode": "copy",
                    "total_discovery_inputs": 455,
                    "addition_policy_rows": 0,
                    "publication_addition_rows": 0,
                    "boundary_rejection_policy_rows": 3,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        plan_path = evidence / "plan.json"
        plan = {
            "workspace": str(workspace.resolve()),
            "run_id": "local-rms-test",
            "date": "230316",
            "flagnum": 2,
            "workflow_input_sha256": workflow_input_hashes(root),
            "prepared_input_inventory": inventory_summary,
            "canonical_reference": str(canonical_path.resolve()),
            "canonical_reference_sha256": sha256_file(canonical_path),
            "exclusion_policy": str(policy_path.resolve()),
            "exclusion_policy_sha256": sha256_file(policy_path),
            "addition_policy": str(addition_policy_path.resolve()),
            "addition_policy_sha256": sha256_file(addition_policy_path),
            "addition_policy_rows": 0,
            "publication_addition_rows": 0,
            "boundary_rejection_policy": str(boundary_policy_path.resolve()),
            "boundary_rejection_policy_sha256": sha256_file(
                boundary_policy_path
            ),
            "boundary_rejection_policy_rows": 3,
            "preapplication_cumulative_corrections": str(
                preapplication_corrections_path.resolve()
            ),
            "preapplication_cumulative_corrections_sha256": sha256_file(
                preapplication_corrections_path
            ),
            "analysis_utils_snapshot": str(analysis_utils_snapshot.resolve()),
            "analysis_utils_execution_dir": str(
                analysis_utils_execution_dir.resolve()
            ),
            "analysis_utils_sha256": analysis_utils_inventory[
                "file_inventory_sha256"
            ],
            "selected_clean_manifest": str(selected_manifest.resolve()),
            "selected_clean_manifest_sha256": sha256_file(selected_manifest),
        }
        plan_path.write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")
        comparison_path = evidence / "current_corrected_table_comparison.json"
        comparison_path.write_text(
            json.dumps(
                compare_current_corrected_table(
                    generated_path, canonical_path, policy_path, addition_policy_path
                ),
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        snapshot_path = evidence / "ngc4945_fixed_nucleus_check.artifacts.json"
        snapshot_path.write_text(
            json.dumps(
                _artifact_snapshot(root, "ngc4945_fixed_nucleus_check", 2),
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )

        artifacts = {
            "workflow_plan": plan_path,
            "prepared_archive_inventory": preparation_path,
            "prepared_input_integrity": integrity_path,
            "initial_comparison": comparison_path,
            "final_stage_snapshot": snapshot_path,
            "canonical_table": canonical_path,
            "generated_table": generated_path,
            "exclusion_policy": policy_path,
            "addition_policy": addition_policy_path,
            "boundary_rejection_policy": boundary_policy_path,
            "generic_final_noise_audit": noise_path,
            "ngc4945_nuclear_photometry": photometry_path,
            "ngc4945_component_audit": component_path,
            "final_noise_inventory_validation": inventory_validation_path,
            "correction_manifest_transition_validation": correction_transition_path,
            "classification_identity_policy": identity_policy_path,
            "classification_identity_validation": identity_validation_path,
            "analysis_utils_snapshot": analysis_utils_snapshot,
            "preapplication_cumulative_corrections": preapplication_corrections_path,
            "postrun_cumulative_corrections": corrections_path,
            "reviewed_delta_csv": ledger_path,
            "reviewed_addition_csv": addition_ledger_path,
        }
        self.assertEqual(set(artifacts), set(REVIEW_FREEZE_ARTIFACT_NAMES))
        freeze = {
            "schema_version": 3,
            "review_scope": REVIEW_SCOPE,
            "review_status": REVIEW_STATUS,
            "run_id": "local-rms-test",
            "expected_counts": {
                "canonical_rows": 409,
                "reviewed_exclusion_rows": 3,
                "generated_rows": 406,
                "generic_noise_audit_rows": 409,
                "generic_publication_noise_audit_rows": 406,
                "boundary_rejection_policy_rows": 3,
                "boundary_rejected_noise_audit_rows": 3,
                "ngc4945_photometry_rows": 45,
                "ngc4945_component_rows": 7,
                "preapplication_correction_rows": 1,
                "postrun_correction_rows": 1,
                "classification_identity_rows": 1,
                "reviewed_delta_rows": 2,
                "addition_policy_rows": 0,
                "publication_addition_rows": 0,
                "reviewed_addition_rows": 0,
            },
            "artifacts": {
                name: {
                    "path": path.resolve().relative_to(workspace.resolve()).as_posix(),
                    "sha256": sha256_file(path),
                }
                for name, path in artifacts.items()
            },
        }
        freeze_path = evidence / "reviewed_local_scaled_mad_freeze.json"
        freeze_path.write_text(
            json.dumps(freeze, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        args = argparse.Namespace(
            workspace=workspace,
            run_id="local-rms-test",
            review_freeze=freeze_path,
            review_freeze_sha256=sha256_file(freeze_path),
        )
        return args, inventory, canonical_path, generated_path, root / "run_continuum_four_pass.py"

    def test_stage_order_and_mandatory_source18_guard(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            analysis, sma, selected = self.make_workspace(workspace)
            commands = stage_commands(
                python_bin=Path("/casa/python"),
                fitter=analysis / "get_continuum_params.py",
                date="230316",
                flagnum=2,
                calnum=1,
                smapaths=[sma],
                selected_manifest=selected,
                use_selected_clean=True,
            )
            self.assertEqual(
                [name for name, _ in commands[:4]],
                ["prelimchecker", "fit_sizenoloc", "getchoosercoords", "fixtochosenfreq"],
            )
            self.assertEqual(tuple(name for name, _ in STAGES), tuple(name for name, _ in commands[:4]))
            self.assertEqual(commands[4][0], "ngc4945_fixed_nucleus")
            self.assertIn("--apply", commands[4][1])
            self.assertIn("--archive-root", commands[4][1])
            self.assertEqual(commands[5][0], "ngc4945_fixed_nucleus_check")
            self.assertIn("--check", commands[5][1])
            for _, command in commands[:4]:
                self.assertNotIn("--enable_notifications", command)
                moveniter = command.index("--moveniter")
                self.assertEqual(command[moveniter + 1], "0")
            self.assertIn("--use_selected_clean_manifest", commands[3][1])
            for _, command in commands[:3]:
                self.assertNotIn("--use_selected_clean_manifest", command)

    def test_workspace_requires_explicit_marker_and_refuses_git_worktree(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            self.make_workspace(workspace)
            with self.assertRaisesRegex(ValueError, "explicitly initialized"):
                assert_isolated_layout(workspace, require_marker=True)
            marker = initialize_workspace(workspace)
            self.assertEqual(marker.name, MARKER_NAME)
            payload = json.loads(marker.read_text())
            self.assertEqual(Path(payload["workspace"]), workspace.resolve())
            self.assertEqual(
                assert_isolated_layout(workspace, require_marker=True), workspace.resolve()
            )
            (workspace / "analysis/multi_freq_from_archive/get_continuum_params.py").write_text(
                "print('changed after initialization')\n"
            )
            with self.assertRaisesRegex(ValueError, "changed after explicit initialization"):
                assert_isolated_layout(workspace, require_marker=True)

        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            self.make_workspace(workspace)
            (workspace / ".git").mkdir()
            with self.assertRaisesRegex(ValueError, "Git worktree"):
                initialize_workspace(workspace)

    def test_sma_directory_symlink_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            _, sma, selected = self.make_workspace(workspace)
            linked = workspace / "analysis/sma_calibration/linked"
            linked.symlink_to(sma, target_is_directory=True)
            with self.assertRaisesRegex(ValueError, "directory symlink"):
                validate_inputs(workspace, [linked], selected)

    def test_archive_source_directory_symlink_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            analysis, sma, selected = self.make_workspace(workspace)
            outside_source = workspace / "outside_source"
            outside_source.mkdir()
            (analysis / "ALMA/UnsafeSource").symlink_to(
                outside_source, target_is_directory=True
            )
            with self.assertRaisesRegex(ValueError, "directory symlink"):
                validate_inputs(workspace, [sma], selected)

    def test_existing_generated_output_fails_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            analysis, _, _ = self.make_workspace(workspace)
            output = analysis / "images/chosenfreq"
            output.mkdir(parents=True)
            with self.assertRaisesRegex(FileExistsError, "Refusing to overwrite"):
                assert_clean_outputs(analysis, "230316_2", 2)

    def test_workflow_input_hashes_cover_code_and_review_policies(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in WORKFLOW_INPUT_NAMES:
                (root / name).parent.mkdir(parents=True, exist_ok=True)
                (root / name).write_text(f"{name}\n")
            hashes = workflow_input_hashes(root)
            self.assertEqual(set(hashes), set(WORKFLOW_INPUT_NAMES))
            self.assertTrue(all(len(digest) == 64 for digest in hashes.values()))
            self.assertIn("continuum_product_exclusion_policy.csv", hashes)
            self.assertIn("continuum_classification_identity_policy.csv", hashes)

    def test_reviewed_exclusion_policy_fails_closed_on_empty_or_duplicates(self):
        with tempfile.TemporaryDirectory() as temporary:
            policy = Path(temporary) / "policy.csv"
            policy.write_text("product\n")
            with self.assertRaisesRegex(ValueError, "policy is empty"):
                load_reviewed_exclusion_products(policy)

            policy.write_text(
                "product\n"
                "VLA/Test/one.fits\n"
                "VLA/Test/one.fits\n"
                "VLA/Test/three.fits\n"
            )
            with self.assertRaisesRegex(ValueError, "nonblank and unique"):
                load_reviewed_exclusion_products(policy)

    def test_corrected_row_count_is_derived_from_reviewed_policy(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            canonical = root / "canonical.txt"
            generated = root / "generated.txt"
            policy = root / "policy.csv"
            rows = [
                [
                    "Test",
                    "VLA",
                    f"{index:.2f}",
                    "1.00",
                    "2020/01/01_1",
                    "0.10000",
                    "10.000",
                    f"VLA/Test/product_{index:03d}.fits\n",
                ]
                for index in range(409)
            ]
            canonical.write_text("".join(f"{row!r}\n" for row in rows))
            generated.write_text("".join(f"{row!r}\n" for row in rows[:407]))
            policy.write_text(
                "product\n"
                "VLA/Test/product_407.fits\n"
                "VLA/Test/product_408.fits\n"
            )

            report = compare_current_corrected_table(generated, canonical, policy)

            self.assertEqual(report["status"], "exact_match")
            self.assertEqual(report["expected_corrected_rows"], 407)

    def test_casa_preflight_confines_import_logs_to_temporary_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            modules = root / "modules"
            modules.mkdir()
            (modules / "casatasks.py").write_text(
                "from pathlib import Path\nPath('casa-test.log').write_text('log')\n"
            )
            (modules / "casatools.py").write_text("")
            (modules / "analysisUtils.py").write_text("")
            for module in ("astropy", "matplotlib", "numpy", "reproject", "requests"):
                (modules / f"{module}.py").write_text("")
            (modules / "astroquery").mkdir()
            (modules / "astroquery/__init__.py").write_text("")
            (modules / "astroquery/alma.py").write_text("")
            (modules / "astroquery/nvas.py").write_text("")
            preflight_python(Path(sys.executable), {"PYTHONPATH": str(modules)})
            self.assertFalse(any(root.rglob("casa-*.log")))

    def test_analysis_utils_bundle_is_complete_copied_and_hash_locked(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            (source / "analysisUtils.py").write_text("import mpfit\n")
            (source / "mpfit.py").write_text("VALUE = 1\n")
            (source / "helper.py").write_text("VALUE = 2\n")
            (source / "__pycache__").mkdir()
            (source / "__pycache__/helper.pyc").write_bytes(b"ignored")

            expected = analysis_utils_bundle_inventory(source)
            destination = root / "frozen"
            manifest_path = root / "manifest.json"
            manifest = freeze_analysis_utils_bundle(
                source,
                destination,
                manifest_path,
            )

            self.assertEqual(
                {key: manifest[key] for key in expected},
                expected,
            )
            self.assertEqual(analysis_utils_bundle_inventory(destination), expected)
            self.assertFalse((destination / "__pycache__").exists())
            self.assertEqual(
                json.loads(manifest_path.read_text(encoding="utf-8")),
                manifest,
            )

            (destination / "mpfit.py").write_text("VALUE = 3\n")
            self.assertNotEqual(analysis_utils_bundle_inventory(destination), expected)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            (source / "analysisUtils.py").write_text("import mpfit\n")
            (source / "mpfit.py").write_text("VALUE = 1\n")
            (source / "linked.py").symlink_to(source / "mpfit.py")
            with self.assertRaisesRegex(ValueError, "symlink"):
                analysis_utils_bundle_inventory(source)

    # This fixture tests archive-copy isolation only. The independent
    # test_ngc4945_bootstrap exercises real authenticated preparation metadata
    # and the original source18 preservation/merge/transition APIs.
    @patch('run_continuum_four_pass.derive_source18_bootstrap', new=lambda root: ({}, {}))
    @patch('run_continuum_four_pass.initialize_source18_bootstrap', new=lambda *args: {})
    @patch('run_continuum_four_pass.validate_source18_bootstrap_snapshot', new=lambda *args: {})
    @patch('run_continuum_four_pass.load_frozen_fitter_inputs', new=fixture_load)
    @patch('run_continuum_four_pass.validate_and_select_frozen_inputs', new=fixture_select)
    def test_prepare_builds_independent_455_input_workspace(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            release = root / "release"
            analysis = release / "analysis/multi_freq_from_archive"
            analysis.mkdir(parents=True)
            (analysis / "get_continuum_params.py").write_text("print('dummy')\n")
            stale_run_output = (
                analysis
                / "fitsumfiles/noise_audits/"
                "230316_2_final_publication_noise_audit.csv"
            )
            stale_run_output.parent.mkdir(parents=True)
            stale_run_output.write_text("historical run output\n")
            self.write_csv_rows(
                analysis / "continuum_product_addition_policy.csv",
                ADDITION_POLICY_FIELDS,
                [],
            )
            self.write_sorted_rows(
                analysis
                / "fitsumfiles/historical/"
                "fitsummary_fixtochosen_sorted_409_reference.txt",
                [
                    [
                        "Baseline",
                        "ALMA",
                        f"{index}.00",
                        "1.00",
                        "2020/01/01_1",
                        "1.00000",
                        "1.000",
                        f"ALMA/Baseline/baseline_{index:03d}.fits",
                    ]
                    for index in range(409)
                ],
            )
            (release / ".git").mkdir()
            manifest = release / "analysis/sma_calibration/selected_clean_iterations.csv"
            manifest.parent.mkdir(parents=True)
            selected_relative = Path(
                "analysis/sma_calibration/230316_128/images/clean/nitercal_10/"
                "flag2_cal1/fields/science/Test_test"
            )
            manifest.write_text(
                "date_code,flagnum,source,image_path,residual_path\n"
                f"230316,2,Test,./{selected_relative}.image,./{selected_relative}.residual\n"
            )

            archive = root / "archive"
            archive_analysis = archive / "analysis/multi_freq_from_archive"
            for telescope, count in (("ALMA", 133), ("VLA", 307)):
                source = archive_analysis / telescope / "Source"
                source.mkdir(parents=True)
                for index in range(count):
                    (source / f"image_{index:03d}.fits").write_bytes(b"fits")
            for index in range(304, 307):
                (archive_analysis / "VLA/Source" / f"image_{index:03d}.fits").unlink()
            self.write_boundary_fixture(archive_analysis)
            for name in (
                "continuum_boundary_rejection_policy.csv",
                "continuum_boundary_rejection_evidence.csv",
            ):
                (analysis / name).write_bytes((archive_analysis / name).read_bytes())
            discovery_relative = Path(
                "analysis/sma_calibration/230316_128/images/nitercal_0/"
                "flag2_cal1/fields/science"
            )
            discovery = archive / discovery_relative
            discovery.mkdir(parents=True)
            for index in range(15):
                casa = discovery / f"Source{index:02d}_test.image"
                casa.mkdir()
                (casa / "table.dat").write_bytes(b"casa")
            for suffix in ("image", "residual"):
                casa = archive / Path(f"{selected_relative}.{suffix}")
                casa.mkdir(parents=True)
                (casa / "table.dat").write_bytes(suffix.encode())

            write_fixture_manifest(analysis, archive, discovery_relative,
                [Path(f'{selected_relative}.{suffix}') for suffix in ('image','residual')])
            # A later archive-query superset must not expand historical fitter
            # membership. This extra file is deliberately absent from manifest.
            (archive_analysis/'VLA/Source/new_unlisted.fits').write_bytes(b'new archive product')

            workspace = root / "isolated"
            sma_path = prepare_workspace(
                release=release,
                archive_project_root=archive,
                workspace=workspace,
                date="230316",
                flagnum=2,
                sma_discovery_relative=discovery_relative,
            )
            self.assertEqual(sma_path, workspace / discovery_relative)
            self.assertTrue((workspace / MARKER_NAME).is_file())
            self.assertFalse((workspace / ".git").exists())
            self.assertFalse(
                (
                    workspace
                    / "analysis/multi_freq_from_archive/fitsumfiles/noise_audits/"
                    "230316_2_final_publication_noise_audit.csv"
                ).exists()
            )
            archive_links = list(
                (workspace / "analysis/multi_freq_from_archive/ALMA").rglob("*.fits")
            ) + list(
                (workspace / "analysis/multi_freq_from_archive/VLA").rglob("*.fits")
            )
            self.assertEqual(len(archive_links), 440)
            self.assertFalse((workspace/'analysis/multi_freq_from_archive/VLA/Source/new_unlisted.fits').exists())
            self.assertTrue(all(path.is_file() and not path.is_symlink() for path in archive_links))
            self.assertFalse(any(path.is_symlink() for path in sma_path.iterdir()))
            self.assertTrue((workspace / Path(f"{selected_relative}.image")).is_dir())
            inventory = validate_prepared_archive_inputs(workspace)
            self.assertEqual(inventory["archive_input_mode"], "copy")
            self.assertEqual(inventory["counts"], {"ALMA": 133, "VLA": 307})
            self.assertEqual(inventory["sma_discovery_inputs"], 15)
            self.assertEqual(inventory["total_discovery_inputs"], 455)
            preparation = json.loads(
                (workspace / ".continuum-fit-inputs.json").read_text()
            )
            self.assertEqual(len(preparation["archive_input_sha256"]), 440)
            self.assertEqual(len(inventory["file_content_sha256"]), 64)
            self.assertEqual(len(preparation["sma_discovery_sha256"]), 15)
            self.assertEqual(len(preparation["selected_clean_sha256"]), 2)

            discovery_table = sma_path / "Source00_test.image/table.dat"
            discovery_original = discovery_table.read_bytes()
            discovery_table.write_bytes(b"mutated")
            with self.assertRaisesRegex(ValueError, "SMA discovery tree hashes"):
                validate_prepared_archive_inputs(workspace)
            discovery_table.write_bytes(discovery_original)

            selected_table = workspace / Path(f"{selected_relative}.image/table.dat")
            selected_original = selected_table.read_bytes()
            selected_table.write_bytes(b"mutated")
            with self.assertRaisesRegex(ValueError, "selected-clean CASA tree hashes"):
                validate_prepared_archive_inputs(workspace)
            selected_table.write_bytes(selected_original)

            manifest_original = (
                workspace / "analysis/sma_calibration/selected_clean_iterations.csv"
            ).read_bytes()
            (workspace / "analysis/sma_calibration/selected_clean_iterations.csv").write_bytes(
                manifest_original + b"\n"
            )
            with self.assertRaisesRegex(ValueError, "selected-clean manifest hash"):
                validate_prepared_archive_inputs(workspace)

    def test_corrected_table_comparison_filters_only_reviewed_exclusions(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            canonical = root / "canonical.txt"
            generated = root / "generated.txt"
            policy = root / "policy.csv"
            rows = [
                [
                    "Test",
                    "ALMA",
                    f"{index:.2f}",
                    "1.00",
                    "2020/01/01_1",
                    "0.10000",
                    "10.000",
                    f"ALMA/Test/product_{index:03d}.fits\n",
                ]
                for index in range(409)
            ]
            canonical.write_text("".join(f"{row!r}\n" for row in rows))
            generated.write_text("".join(f"{row!r}\n" for row in rows[:406]))
            policy.write_text(
                "product\n"
                "ALMA/Test/product_406.fits\n"
                "ALMA/Test/product_407.fits\n"
                "ALMA/Test/product_408.fits\n"
            )
            report = compare_current_corrected_table(generated, canonical, policy)
            self.assertEqual(report["status"], "exact_match")
            self.assertEqual(report["canonical_rows"], 409)
            self.assertEqual(report["generated_rows"], 406)
            self.assertEqual(report["expected_corrected_rows"], 406)

            changed = list(rows[:406])
            changed[10] = list(changed[10])
            changed[10][3] = "2.00"
            generated.write_text("".join(f"{row!r}\n" for row in changed))
            report = compare_current_corrected_table(generated, canonical, policy)
            self.assertEqual(report["status"], "review_required")
            self.assertEqual(report["field_mismatch_count"], 1)

    def test_corrected_table_comparison_matches_absolute_sma_to_relative_product(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            canonical = root / "canonical.txt"
            generated = root / "generated.txt"
            policy = root / "policy.csv"
            rows = [
                [
                    "Test",
                    "ALMA",
                    f"{index:.2f}",
                    "1.00",
                    "2020/01/01_1",
                    "0.10000",
                    "10.000",
                    f"ALMA/Test/product_{index:03d}.fits\n",
                ]
                for index in range(409)
            ]
            sma_product = (
                "analysis/sma_calibration/230316_128/images/nitercal_0/"
                "flag2_cal1/fields/science/Test_test.image"
            )
            rows[0][1] = "SMA"
            rows[0][7] = f"{sma_product}\n"
            generated_rows = [list(row) for row in rows[:406]]
            generated_rows[0][7] = f"/home/author/release/{sma_product}\n"
            canonical.write_text("".join(f"{row!r}\n" for row in rows))
            generated.write_text(
                "".join(f"{row!r}\n" for row in generated_rows)
            )
            policy.write_text(
                "product\n"
                "ALMA/Test/product_406.fits\n"
                "ALMA/Test/product_407.fits\n"
                "ALMA/Test/product_408.fits\n"
            )

            report = compare_current_corrected_table(generated, canonical, policy)

            self.assertEqual(report["status"], "exact_match")
            self.assertEqual(report["missing_products"], [])
            self.assertEqual(report["extra_products"], [])

    def test_corrected_table_comparison_rejects_unknown_absolute_product(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            canonical = root / "canonical.txt"
            generated = root / "generated.txt"
            policy = root / "policy.csv"
            rows = [
                [
                    "Test",
                    "ALMA",
                    f"{index:.2f}",
                    "1.00",
                    "2020/01/01_1",
                    "0.10000",
                    "10.000",
                    f"ALMA/Test/product_{index:03d}.fits\n",
                ]
                for index in range(409)
            ]
            generated_rows = [list(row) for row in rows[:406]]
            generated_rows[0][7] = "/tmp/untraceable/product.fits\n"
            canonical.write_text("".join(f"{row!r}\n" for row in rows))
            generated.write_text(
                "".join(f"{row!r}\n" for row in generated_rows)
            )
            policy.write_text(
                "product\n"
                "ALMA/Test/product_406.fits\n"
                "ALMA/Test/product_407.fits\n"
                "ALMA/Test/product_408.fits\n"
            )

            with self.assertRaisesRegex(ValueError, "unknown absolute"):
                compare_current_corrected_table(generated, canonical, policy)

    def test_initial_comparator_records_product_order_drift(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            canonical, generated, policy, _, _, generated_rows = (
                self.make_reviewed_delta_fixture(root)
            )
            generated_rows[2], generated_rows[3] = generated_rows[3], generated_rows[2]
            self.write_sorted_rows(generated, generated_rows)

            report = compare_current_corrected_table(generated, canonical, policy)

            self.assertEqual(report["status"], "review_required")
            self.assertFalse(report["product_order_matches"])
            self.assertEqual(
                report["first_order_mismatch"]["corrected_row_number"], 3
            )

    def test_exact_reviewed_current_to_generated_ledger_passes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            canonical, generated, policy, ledger, _, _ = (
                self.make_reviewed_delta_fixture(root)
            )

            initial = compare_current_corrected_table(generated, canonical, policy)
            adjudicated = compare_reviewed_current_to_generated(
                generated, canonical, policy, ledger
            )

            self.assertEqual(initial["status"], "review_required")
            self.assertTrue(initial["product_order_matches"])
            self.assertEqual(adjudicated["status"], "exact_reviewed_match")
            self.assertEqual(adjudicated["validated_product_identities"], 406)
            self.assertEqual(adjudicated["changed_product_count"], 2)
            self.assertEqual(adjudicated["flux_field_change_count"], 1)
            self.assertEqual(adjudicated["snr_change_count"], 2)

    def test_reviewed_delta_ledger_fails_closed_on_inventory_and_field_mutations(self):
        mutations = {
            "missing actual delta": lambda rows: rows[:1],
            "stale extra delta": lambda rows: rows
            + [
                {
                    **rows[-1],
                    "corrected_row_number": "12",
                    "frequency_ghz": "11.00",
                    "product": "ALMA/Test/product_011.fits",
                    "current_snr": "10.000",
                    "generated_snr": "10.000",
                }
            ],
            "wrong current field": lambda rows: [
                {**rows[0], "current_flux_field": "9.99"}, rows[1]
            ],
            "wrong generated field": lambda rows: [
                {**rows[0], "generated_snr": "99.999"}, rows[1]
            ],
            "wrong reason": lambda rows: [
                {**rows[0], "review_reason": "uncontrolled_reason"}, rows[1]
            ],
            "wrong row number": lambda rows: [
                {**rows[0], "corrected_row_number": "2"}, rows[1]
            ],
            "reordered approvals": lambda rows: list(reversed(rows)),
            "duplicate product": lambda rows: [rows[0], dict(rows[0])],
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                canonical, generated, policy, ledger, _, _ = (
                    self.make_reviewed_delta_fixture(root)
                )
                with ledger.open(encoding="utf-8", newline="") as handle:
                    rows = list(csv.DictReader(handle))
                self.write_csv_rows(ledger, REVIEW_DELTA_FIELDS, mutate(rows))
                with self.assertRaises(ValueError):
                    compare_reviewed_current_to_generated(
                        generated, canonical, policy, ledger
                    )

    def test_reviewed_adjudication_never_allows_identity_order_or_other_field_drift(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            canonical, generated, policy, ledger, _, generated_rows = (
                self.make_reviewed_delta_fixture(root)
            )
            generated_rows[1][5] = "0.20000"
            self.write_sorted_rows(generated, generated_rows)
            with self.assertRaisesRegex(ValueError, "non-flux/SNR"):
                compare_reviewed_current_to_generated(
                    generated, canonical, policy, ledger
                )

            generated_rows[1][5] = "0.10000"
            generated_rows[1], generated_rows[2] = generated_rows[2], generated_rows[1]
            self.write_sorted_rows(generated, generated_rows)
            with self.assertRaisesRegex(ValueError, "order differs"):
                compare_reviewed_current_to_generated(
                    generated, canonical, policy, ledger
                )

    def test_correction_manifest_transition_is_exactly_source18_only(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            before = root / "before.csv"
            after = root / "after.csv"
            photometry = root / "photometry.csv"
            legacy_rows = []
            v2_rows = []
            photometry_rows = []
            for index in range(66):
                is_ngc4945 = index >= 39
                source = "NGC4945" if is_ngc4945 else "Test"
                product = f"ALMA/{source}/product_{index:03d}.fits"
                legacy = {
                    "source": source,
                    "telescope": "ALMA",
                    "frequency_ghz": f"{100.0 + index:.2f}",
                    "date_label": "2020/01/01_1",
                    "product": product,
                    "old_flux_field": "1.00",
                    "new_flux_field": "1.10",
                    "old_snr": "9.000",
                    "new_snr": "10.000",
                    "correction_reasons": (
                        "fixed_nuclear_component_source18+single_background_rms"
                        if is_ngc4945
                        else "header_aware_target_frame"
                    ),
                    "validated_unrounded_flux_mjy": "1.1",
                    "validated_background_only_snr": "10.000",
                }
                legacy_rows.append(legacy)
                v2 = {
                    "schema_version": "2",
                    **{
                        field: legacy[field]
                        for field in LEGACY_CORRECTION_FIELDS[:-1]
                    },
                    "validated_adopted_noise_snr": "10.000",
                }
                if is_ngc4945:
                    v2.update(
                        {
                            "new_flux_field": "1.23",
                            "new_snr": "12.345",
                            "correction_reasons": (
                                "fixed_nuclear_component_source18+"
                                "single_local_scaled_mad+"
                                "source18_10bmaj_local_noise_aperture"
                            ),
                            "validated_unrounded_flux_mjy": "1.2345",
                            "validated_adopted_noise_snr": "12.345",
                        }
                    )
                    measurement = {field: "" for field in PHOTOMETRY_FIELDS}
                    measurement.update(
                        {
                            "product": product,
                            "old_flux_field": "1.00",
                            "old_snr": "9.000",
                            "new_flux_field": "1.23",
                            "new_snr": "12.345",
                            "unrounded_flux_mjy": "1.2345",
                        }
                    )
                    photometry_rows.append(measurement)
                v2_rows.append(v2)

            # Task 28 adds 18 NGC 4945 products after the 66-row historical
            # correction snapshot.  They enter only the post-run manifest and
            # the expanded source-18 photometry ledger.
            for index in range(66, 84):
                product = f"ALMA/NGC4945/product_{index:03d}.fits"
                v2 = {
                    "schema_version": "2",
                    "source": "NGC4945",
                    "telescope": "ALMA",
                    "frequency_ghz": f"{100.0 + index:.2f}",
                    "date_label": "2020/01/01_1",
                    "product": product,
                    "old_flux_field": "1.00",
                    "new_flux_field": "1.23",
                    "old_snr": "9.000",
                    "new_snr": "12.345",
                    "correction_reasons": (
                        "fixed_nuclear_component_source18+"
                        "single_local_scaled_mad+"
                        "source18_10bmaj_local_noise_aperture"
                    ),
                    "validated_unrounded_flux_mjy": "1.2345",
                    "validated_adopted_noise_snr": "12.345",
                }
                measurement = {field: "" for field in PHOTOMETRY_FIELDS}
                measurement.update(
                    {
                        "product": product,
                        "old_flux_field": "1.00",
                        "old_snr": "9.000",
                        "new_flux_field": "1.23",
                        "new_snr": "12.345",
                        "unrounded_flux_mjy": "1.2345",
                    }
                )
                v2_rows.append(v2)
                photometry_rows.append(measurement)

            self.write_csv_rows(before, LEGACY_CORRECTION_FIELDS, legacy_rows)
            self.write_csv_rows(after, CORRECTION_FIELDS, v2_rows)
            self.write_csv_rows(photometry, PHOTOMETRY_FIELDS, photometry_rows)
            report = validate_correction_manifest_transition(
                before, after, photometry
            )
            self.assertEqual(report["status"], "exact_source18_only_transition")
            self.assertEqual(report["unchanged_non_ngc4945_rows"], 39)
            self.assertEqual(report["added_ngc4945_rows"], 18)
            self.assertEqual(report["reviewed_ngc4945_rows"], 45)

            changed_non_ngc = [dict(row) for row in v2_rows]
            changed_non_ngc[0]["new_snr"] = "99.000"
            self.write_csv_rows(after, CORRECTION_FIELDS, changed_non_ngc)
            with self.assertRaisesRegex(ValueError, "Non-NGC4945"):
                validate_correction_manifest_transition(before, after, photometry)

            wrong_reason = [dict(row) for row in v2_rows]
            wrong_reason[39]["correction_reasons"] = "hidden_override"
            self.write_csv_rows(after, CORRECTION_FIELDS, wrong_reason)
            with self.assertRaisesRegex(ValueError, "reasons"):
                validate_correction_manifest_transition(before, after, photometry)

            wrong_photometry = [dict(row) for row in photometry_rows]
            wrong_photometry[0]["new_snr"] = "12.344"
            self.write_csv_rows(after, CORRECTION_FIELDS, v2_rows)
            self.write_csv_rows(
                photometry, PHOTOMETRY_FIELDS, wrong_photometry
            )
            with self.assertRaisesRegex(ValueError, "new_snr"):
                validate_correction_manifest_transition(before, after, photometry)

            wrong_old_lineage = [dict(row) for row in photometry_rows]
            wrong_old_lineage[0]["old_snr"] = "999.000"
            self.write_csv_rows(
                photometry, PHOTOMETRY_FIELDS, wrong_old_lineage
            )
            with self.assertRaisesRegex(ValueError, "old_snr"):
                validate_correction_manifest_transition(before, after, photometry)

    def test_review_freeze_requires_exact_schema_and_explicit_sha(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "freeze.json"
            payload = {
                "schema_version": 3,
                "review_scope": REVIEW_SCOPE,
                "review_status": REVIEW_STATUS,
                "run_id": "run-1",
                "expected_counts": {
                    "canonical_rows": 409,
                    "reviewed_exclusion_rows": 3,
                    "generated_rows": 406,
                    "generic_noise_audit_rows": 409,
                    "generic_publication_noise_audit_rows": 406,
                    "boundary_rejection_policy_rows": 3,
                    "boundary_rejected_noise_audit_rows": 3,
                    "ngc4945_photometry_rows": 45,
                    "ngc4945_component_rows": 7,
                    "preapplication_correction_rows": 66,
                    "postrun_correction_rows": 67,
                    "classification_identity_rows": 12,
                    "reviewed_delta_rows": 2,
                    "addition_policy_rows": 0,
                    "publication_addition_rows": 0,
                    "reviewed_addition_rows": 0,
                },
                "artifacts": {
                    name: {
                        "path": f"analysis/review/{name}",
                        "sha256": "0" * 64,
                    }
                    for name in REVIEW_FREEZE_ARTIFACT_NAMES
                },
            }
            path.write_text(json.dumps(payload) + "\n", encoding="utf-8")
            digest = sha256_file(path)
            loaded = load_review_freeze(path, digest, "run-1")
            self.assertEqual(loaded, payload)

            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                load_review_freeze(path, "f" * 64, "run-1")
            payload["unexpected"] = True
            path.write_text(json.dumps(payload) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "exact schema"):
                load_review_freeze(path, sha256_file(path), "run-1")

    def test_adjudicate_is_evidence_only_atomic_and_never_calls_casa(self):
        with tempfile.TemporaryDirectory() as temporary:
            args, inventory, canonical, generated, workspace_runner = (
                self.make_adjudication_fixture(Path(temporary))
            )
            canonical_before = canonical.read_bytes()
            generated_before = generated.read_bytes()
            with (
                patch(
                    "run_continuum_four_pass.validate_prepared_archive_inputs",
                    return_value=inventory,
                ),
                patch("run_continuum_four_pass.__file__", str(workspace_runner)),
                patch(
                    "run_continuum_four_pass.preflight_python",
                    side_effect=AssertionError("CASA preflight must not run"),
                ),
                patch(
                    "run_continuum_four_pass.subprocess.run",
                    side_effect=AssertionError("no subprocess may run"),
                ),
                patch(
                    "run_continuum_four_pass.validate_final_noise_inventory",
                    return_value={
                        "schema_version": 1,
                        "status": "exact_validated_union",
                        "publication_rows": 406,
                        "admitted_input_additions": 0,
                        "publication_additions": 0,
                        "input_only_additions": 0,
                        "generic_audit_rows": 409,
                        "generic_publication_product_rows": 406,
                        "ordinary_generic_rows": 397,
                        "ordinary_publication_product_rows": 394,
                        "forced_dirty_sma_rows": 12,
                        "boundary_rejection_policy_rows": 3,
                        "boundary_rejected_noise_audit_rows": 3,
                        "generic_ngc4945_rejected_rows": 45,
                        "source18_retained_rows": 45,
                        "retained_union_rows": 406,
                    },
                ) as validate_inventory,
                patch(
                    "run_continuum_four_pass.validate_ngc4945_component_artifacts",
                    return_value={
                        "ngc4945_photometry_semantic_status": "exact_validated",
                        "ngc4945_photometry_semantic_rows": 45,
                        "ngc4945_component_semantic_status": "exact_validated",
                        "ngc4945_component_semantic_rows": 7,
                        "ngc4945_component_record_ids": [
                            f"c{index}" for index in range(7)
                        ],
                    },
                ) as validate_components,
                patch(
                    "run_continuum_four_pass.validate_correction_manifest_transition",
                    return_value={
                        "schema_version": 1,
                        "status": "exact_source18_only_transition",
                        "preapplication_rows": 66,
                        "postrun_rows": 84,
                        "unchanged_non_ngc4945_rows": 39,
                        "added_ngc4945_rows": 18,
                        "reviewed_ngc4945_rows": 45,
                        "reviewed_ngc4945_products": [
                            f"p{index}" for index in range(45)
                        ],
                    },
                ) as validate_transition,
                patch(
                    "run_continuum_four_pass.validate_classification_identity_policy",
                    return_value={
                        "schema_version": 1,
                        "status": "exact_required_identity_reviews_resolved",
                        "required_review_rows": 1,
                        "new_detection_rows": 1,
                        "rms_mad_threshold_crossing_rows": 0,
                        "new_product_admission_rows": 0,
                        "admitted_product_upper_limit_rows": 0,
                        "validated_reviews": [{"review_id": "fixture"}],
                    },
                ) as validate_identity,
            ):
                self.assertEqual(adjudicate_workflow(args), 0)
                with self.assertRaises(FileExistsError):
                    adjudicate_workflow(args)

            self.assertEqual(validate_inventory.call_count, 2)
            self.assertEqual(validate_components.call_count, 2)
            self.assertEqual(validate_transition.call_count, 2)
            self.assertEqual(validate_identity.call_count, 2)

            self.assertEqual(canonical.read_bytes(), canonical_before)
            self.assertEqual(generated.read_bytes(), generated_before)
            output = (
                args.workspace
                / "analysis/multi_freq_from_archive/workflow_runs/local-rms-test/"
                "reviewed_delta_adjudication.json"
            )
            evidence = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(evidence["status"], "exact_reviewed_match")
            self.assertEqual(
                evidence["delta_comparison"]["validated_product_identities"], 406
            )

    def test_adjudicate_rejects_any_frozen_audit_input_hash_drift(self):
        with tempfile.TemporaryDirectory() as temporary:
            args, inventory, _, _, workspace_runner = self.make_adjudication_fixture(
                Path(temporary)
            )
            noise = (
                args.workspace
                / "analysis/multi_freq_from_archive/fitsumfiles/noise_audits/"
                "230316_2_final_publication_noise_audit.csv"
            )
            noise.write_text("product\nmutated\n", encoding="utf-8")
            with (
                patch(
                    "run_continuum_four_pass.validate_prepared_archive_inputs",
                    return_value=inventory,
                ),
                patch("run_continuum_four_pass.__file__", str(workspace_runner)),
            ):
                with self.assertRaisesRegex(ValueError, "SHA-256 differs"):
                    adjudicate_workflow(args)


if __name__ == "__main__":
    unittest.main()
