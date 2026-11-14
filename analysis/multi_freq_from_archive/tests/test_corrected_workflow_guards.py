from __future__ import annotations

import ast
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ANALYSIS_DIR.parents[1]
sys.path.insert(0, str(ANALYSIS_DIR))
import continuum_boundary_rejections as boundary_rejections
from continuum_measurements import portable_release_product


class CorrectedWorkflowGuardTests(unittest.TestCase):
    @staticmethod
    def _load_fitter_function(name: str):
        fitter_path = ANALYSIS_DIR / "get_continuum_params.py"
        tree = ast.parse(fitter_path.read_text(), filename=str(fitter_path))
        function = next(
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name == name
        )
        namespace = {
            "csv": csv,
            "hashlib": hashlib,
            "os": os,
            "Path": Path,
            "portable_release_product": portable_release_product,
            "shutil": shutil,
        }
        exec(
            compile(
                ast.fix_missing_locations(ast.Module(body=[function], type_ignores=[])),
                str(fitter_path),
                "exec",
            ),
            namespace,
        )
        return namespace[name]

    def test_black_hole_export_uses_the_deterministic_builder(self) -> None:
        exporter = (ANALYSIS_DIR / "export_bhs_params.py").read_text()
        notebook = (REPOSITORY_ROOT / "BH_Parameters.ipynb").read_text()

        self.assertIn("build_bhs_params.py", exporter)
        self.assertIn('"--require-final"', exporter)
        self.assertNotIn("exec(compile", exporter)
        self.assertNotIn("NOTEBOOK_PATH", exporter)
        self.assertNotIn("rho=1", notebook)
        self.assertNotIn("add_uncorrelated", notebook)

        inputs = REPOSITORY_ROOT / "analysis" / "CountingBHS_and_Spectra"
        self.assertTrue((inputs / "build_bhs_params.py").is_file())
        self.assertTrue((inputs / "bhs_source_inputs.csv").is_file())
        self.assertTrue((inputs / "pesce_theta_g_quantiles.csv").is_file())

    def test_black_hole_release_export_accepts_explicit_formal_only_treatment(self) -> None:
        result = subprocess.run(
            [sys.executable, str(ANALYSIS_DIR / "export_bhs_params.py")],
            cwd=REPOSITORY_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Refreshed", result.stdout)

    def test_post_processing_uses_corrected_builders_without_size_regeneration(self) -> None:
        path = ANALYSIS_DIR / "post_processing_finalcontinuum.ipynb"
        notebook = json.loads(path.read_text())
        source = "".join(
            line
            for cell in notebook["cells"]
            for line in cell.get("source", [])
        )
        self.assertIn("build_variability_tables.py", source)
        self.assertIn("build_contamination_tables.py", source)
        self.assertIn("rerun_ngc4945_nuclear.py", source)
        self.assertIn("--apply-recorded", source)
        self.assertNotIn("subprocess.run(['python3', 'combine_fitsummaries.py']", source)
        self.assertNotIn("subprocess.run(['python3', 'export_bhs_params.py']", source)
        self.assertNotIn('with open("Variability_Search.ipynb")', source)
        self.assertNotIn('with open("get_contams.ipynb")', source)
        self.assertNotIn('with open("flux_size_graph.ipynb")', source)
        self.assertNotIn('with open("Beam_Date_Flux_Graphs.ipynb")', source)

    def test_legacy_notebooks_cannot_overwrite_corrected_machine_products(self) -> None:
        variability = (ANALYSIS_DIR / "Variability_Search.ipynb").read_text()
        contamination = (ANALYSIS_DIR / "get_contams.ipynb").read_text()
        self.assertNotIn("corr_table='machinetables/correlations_machine.txt'", variability)
        self.assertNotIn("corr_table='machinetables/variability_machine.txt'", variability)
        self.assertNotIn(
            "output_csv_path='machinetables/S230_extrapolation_comparison.csv'",
            contamination,
        )
        self.assertIn("legacy_contam_", contamination)

    def test_alma_discovery_does_not_treat_eso_source_names_as_so_line_products(self) -> None:
        notebook = json.loads((ANALYSIS_DIR / "get_urls.ipynb").read_text())
        source = "".join(
            line
            for cell in notebook["cells"]
            for line in cell.get("source", [])
        )
        self.assertIn("product_name=os.path.basename(url)", source)
        self.assertIn("re.search(r'(^|[^A-Za-z])SO([^A-Za-z]|$)'", source)
        self.assertNotIn("if chem in url:", source)

    def test_coordinate_reader_accepts_sources_without_alternate_names(self) -> None:
        coordinate_rows = ast.literal_eval(
            (ANALYSIS_DIR / "namesandcoords.txt").read_text()
        )
        circinus = next(
            row for row in coordinate_rows if row[0][0] == "CircinusGalaxy"
        )
        self.assertEqual(circinus[1][0], [])

        fitter = (ANALYSIS_DIR / "get_continuum_params.py").read_text()
        self.assertIn("nname=inamecoords[1][0]", fitter)
        self.assertNotIn("nname=inamecoords[1][0][0]", fitter)

    def test_failed_ngc4258_nvas_products_are_skipped_not_forced_to_limits(self) -> None:
        fitter = (ANALYSIS_DIR / "get_continuum_params.py").read_text()
        self.assertIn("reviewed_archive_photometry_exclusions", fitter)
        self.assertIn("load_reviewed_exclusion_policy(", fitter)
        self.assertIn("verify_reviewed_exclusion_hash(file", fitter)
        self.assertIn("continuum_product_exclusion_policy.csv", fitter)
        self.assertNotIn("bad1=", fitter)
        self.assertNotIn("totmasterbad", fitter)
        self.assertNotIn(
            "VLA/NGC4258/22.4I1.92_AH0847_2004MAY22_1_2.38M55.7S.imfits",
            fitter,
        )
        self.assertNotIn("lowfiles=", fitter)
        self.assertNotIn("for ilow in lowfiles", fitter)

    def test_reviewed_exclusion_hash_guard_is_policy_driven_and_fail_closed(self) -> None:
        load_policy = self._load_fitter_function("load_reviewed_exclusion_policy")
        verify_hash = self._load_fitter_function("verify_reviewed_exclusion_hash")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            product = root / "artifact.fits"
            product.write_bytes(b"reviewed artifact")
            digest = hashlib.sha256(product.read_bytes()).hexdigest()
            policy = root / "policy.csv"
            policy.write_text(
                "exclusion_id,product,sha256,decision_reason\n"
                f"test,VLA/Test/artifact.fits,{digest},structured_artifact\n"
            )

            self.assertEqual(
                load_policy(policy),
                {"VLA/Test/artifact.fits": digest},
            )
            verify_hash(product, digest)
            product.write_bytes(b"changed artifact")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                verify_hash(product, digest)

    def test_boundary_rejections_are_policy_driven_at_exactly_two_sites(self) -> None:
        fitter_path = ANALYSIS_DIR / "get_continuum_params.py"
        fitter = fitter_path.read_text()
        tree = ast.parse(fitter, filename=str(fitter_path))

        imports = [
            node
            for node in tree.body
            if isinstance(node, ast.ImportFrom)
            and node.module == "continuum_boundary_rejections"
        ]
        self.assertEqual(len(imports), 1)
        self.assertTrue(
            {
                "load_boundary_rejection_policy",
                "validate_boundary_rejection_record_identity",
            }
            <= {alias.name for alias in imports[0].names}
        )

        calls = [node for node in ast.walk(tree) if isinstance(node, ast.Call)]
        loader_calls = [
            node
            for node in calls
            if isinstance(node.func, ast.Name)
            and node.func.id == "load_boundary_rejection_policy"
        ]
        self.assertEqual(len(loader_calls), 1)
        loader_call = loader_calls[0]
        self.assertEqual(
            {keyword.arg for keyword in loader_call.keywords}, {"analysis_root"}
        )
        self.assertIn(
            "continuum_boundary_rejection_policy.csv",
            ast.unparse(loader_call.args[0]),
        )
        self.assertEqual(
            ast.unparse(loader_call.keywords[0].value), "Path('.')"
        )

        guarded_calls = [
            node
            for node in calls
            if isinstance(node.func, ast.Name)
            and node.func.id == "record_reviewed_boundary_rejection"
        ]
        self.assertEqual(len(guarded_calls), 2)
        self.assertTrue(
            all(
                isinstance(node.args[0], ast.Name) and node.args[0].id == "file"
                for node in guarded_calls
            )
        )

        final_fit_guards = []
        for node in ast.walk(tree):
            if not (
                isinstance(node, ast.If)
                and isinstance(node.test, ast.Name)
                and node.test.id == "fixtochosenfreq"
            ):
                continue
            direct_guard_calls = [
                statement.value
                for statement in node.body
                if isinstance(statement, ast.Expr)
                and isinstance(statement.value, ast.Call)
                and isinstance(statement.value.func, ast.Name)
                and statement.value.func.id == "record_reviewed_boundary_rejection"
            ]
            final_fit_guards.extend(direct_guard_calls)
        self.assertEqual(
            {node.lineno for node in final_fit_guards},
            {node.lineno for node in guarded_calls},
        )

        direct_dispositions = [
            node
            for node in calls
            if isinstance(node.func, ast.Name)
            and node.func.id == "update_final_noise_audit"
            and any(
                keyword.arg == "row_disposition"
                and isinstance(keyword.value, ast.Constant)
                and keyword.value.value == "rejected_on_boundary"
                for keyword in node.keywords
            )
        ]
        self.assertEqual(len(direct_dispositions), 1)
        helper = next(
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef)
            and node.name == "record_reviewed_boundary_rejection"
        )
        self.assertIn(direct_dispositions[0], list(ast.walk(helper)))

        policy_path = ANALYSIS_DIR / "continuum_boundary_rejection_policy.csv"
        with policy_path.open(encoding="utf-8", newline="") as handle:
            products = [row["product"] for row in csv.DictReader(handle)]
        self.assertEqual(len(products), 3)
        for product in products:
            self.assertNotIn(product, fitter)

        # Policy membership is checked only after a fit reports ONBOUNDRY. It
        # must never become a discovery-loop pre-skip condition.
        for node in ast.walk(tree):
            if not isinstance(node, ast.Compare) or not any(
                isinstance(operator, (ast.In, ast.NotIn))
                for operator in node.ops
            ):
                continue
            self.assertNotIn(
                "reviewed_boundary_rejections_by_product", ast.unparse(node)
            )

    def test_boundary_disposition_guard_rejects_unreviewed_duplicate_and_stale_rows(
        self,
    ) -> None:
        guard = self._load_fitter_function("record_reviewed_boundary_rejection")
        product = "VLA/Test/allowed.imfits"
        policy = {
            "product": product,
            "source": "Test",
            "telescope": "VLA",
            "expected_frequency_ghz": "1.49",
            "expected_row_disposition": "rejected_on_boundary",
            "expected_classification": "detection",
        }
        audit_records = {
            product: {
                "product": product,
                "source": "Test",
                "telescope": "VLA",
                "frequency_ghz": "1.49",
                "row_disposition": "pending",
                "classification": "detection",
            }
        }

        def update(product_path: str, **updates: str) -> None:
            normalized = portable_release_product(product_path)
            audit_records[normalized].update(updates)

        observed: set[str] = set()
        guard.__globals__.update(
            {
                "reviewed_boundary_rejections_by_product": {product: policy},
                "observed_boundary_rejection_products": observed,
                "final_noise_audit_records": audit_records,
                "update_final_noise_audit": update,
                "validate_boundary_rejection_record_identity": (
                    boundary_rejections.validate_boundary_rejection_record_identity
                ),
            }
        )

        guard(product)
        self.assertEqual(observed, {product})
        self.assertEqual(
            audit_records[product]["row_disposition"], "rejected_on_boundary"
        )
        with self.assertRaisesRegex(RuntimeError, "Duplicate final-fit"):
            guard(product)
        with self.assertRaisesRegex(RuntimeError, "Unreviewed final-fit"):
            guard("VLA/Test/unreviewed.imfits")

        observed.clear()
        audit_records[product]["source"] = "WrongSource"
        with self.assertRaisesRegex(RuntimeError, "identity/disposition is stale"):
            guard(product)
        self.assertEqual(observed, set())

    def test_boundary_observed_inventory_guard_requires_exact_set_equality(self) -> None:
        fitter_path = ANALYSIS_DIR / "get_continuum_params.py"
        tree = ast.parse(fitter_path.read_text(), filename=str(fitter_path))
        inventory_guards = [
            node
            for node in tree.body
            if isinstance(node, ast.If)
            and isinstance(node.test, ast.Name)
            and node.test.id == "fixtochosenfreq"
            and any(
                isinstance(child, ast.Assign)
                and any(
                    isinstance(target, ast.Name)
                    and target.id == "expected_boundary_rejections"
                    for target in child.targets
                )
                for child in node.body
            )
        ]
        self.assertEqual(len(inventory_guards), 1)
        inventory_guard = inventory_guards[0]
        exact_comparisons = [
            node
            for node in ast.walk(inventory_guard)
            if isinstance(node, ast.Compare)
            and isinstance(node.left, ast.Name)
            and node.left.id == "observed_boundary_rejection_products"
            and len(node.ops) == 1
            and isinstance(node.ops[0], ast.NotEq)
            and len(node.comparators) == 1
            and isinstance(node.comparators[0], ast.Name)
            and node.comparators[0].id == "expected_boundary_rejections"
        ]
        self.assertEqual(len(exact_comparisons), 1)

        executable = compile(
            ast.fix_missing_locations(
                ast.Module(body=[inventory_guard], type_ignores=[])
            ),
            str(fitter_path),
            "exec",
        )
        expected = {"VLA/Test/one.imfits", "VLA/Test/two.imfits"}
        for label, observed, should_raise in (
            ("exact", set(expected), False),
            ("missing", {"VLA/Test/one.imfits"}, True),
            ("extra", expected | {"VLA/Test/extra.imfits"}, True),
        ):
            namespace = {
                "fixtochosenfreq": True,
                "reviewed_boundary_rejections_by_product": {
                    product: {} for product in expected
                },
                "observed_boundary_rejection_products": observed,
            }
            with self.subTest(label=label):
                if should_raise:
                    with self.assertRaisesRegex(
                        RuntimeError,
                        "outcome inventory differs from policy",
                    ):
                        exec(executable, namespace)
                else:
                    exec(executable, namespace)

    def test_boundary_policy_product_hash_drift_fails_closed(self) -> None:
        source_policy = ANALYSIS_DIR / "continuum_boundary_rejection_policy.csv"
        source_evidence = ANALYSIS_DIR / "continuum_boundary_rejection_evidence.csv"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            policy_path = root / source_policy.name
            evidence_path = root / source_evidence.name
            shutil.copyfile(source_evidence, evidence_path)
            with source_policy.open(encoding="utf-8", newline="") as handle:
                rows = [dict(row) for row in csv.DictReader(handle)]

            identities = {}
            for row in rows:
                product_path = root / row["product"]
                product_path.parent.mkdir(parents=True, exist_ok=True)
                payload = f"immutable-{row['rejection_id']}".encode("ascii")
                product_path.write_bytes(payload)
                row["size_bytes"] = str(len(payload))
                row["product_sha256"] = hashlib.sha256(payload).hexdigest()
                identities[row["product"]] = (
                    row["expected_frequency_ghz"],
                    row["expected_observation_date"],
                )
            with policy_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=boundary_rejections.BOUNDARY_REJECTION_POLICY_FIELDS,
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(rows)

            def header_identity(path: Path) -> tuple[str, str]:
                product = path.resolve().relative_to(root.resolve()).as_posix()
                return identities[product]

            with mock.patch.object(
                boundary_rejections,
                "fits_header_output_identity",
                side_effect=header_identity,
            ):
                loaded = boundary_rejections.load_boundary_rejection_policy(
                    policy_path,
                    analysis_root=root,
                )
                self.assertEqual(len(loaded), 3)
                changed_path = root / rows[0]["product"]
                changed_path.write_bytes(b"X" * changed_path.stat().st_size)
                with self.assertRaisesRegex(ValueError, "SHA-256 differs"):
                    boundary_rejections.load_boundary_rejection_policy(
                        policy_path,
                        analysis_root=root,
                    )

    def test_coordinate_image_selection_is_bounded_and_fails_closed(self) -> None:
        fitter = (ANALYSIS_DIR / "get_continuum_params.py").read_text()
        self.assertNotIn("while 0==0:", fitter)
        self.assertIn("for _attempt in range(attempt_limit):", fitter)
        self.assertIn("discard_failed_coordinate_candidate_or_raise(", fitter)

    def test_getchooser_bad_files_have_an_unconditional_stage_label(self) -> None:
        tree = ast.parse((ANALYSIS_DIR / "get_continuum_params.py").read_text())
        stage_labels = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.If):
                continue
            test = node.test
            if not (
                isinstance(test, ast.Compare)
                and isinstance(test.left, ast.Name)
                and test.left.id == "getchoosercoords"
                and len(test.ops) == 1
                and isinstance(test.ops[0], ast.Eq)
                and len(test.comparators) == 1
                and isinstance(test.comparators[0], ast.Constant)
                and test.comparators[0].value == 1
            ):
                continue
            for statement in node.body:
                if not isinstance(statement, ast.Assign):
                    continue
                if any(
                    isinstance(target, ast.Name) and target.id == "badfileversion"
                    for target in statement.targets
                ):
                    self.assertIsInstance(statement.value, ast.Constant)
                    stage_labels.append(statement.value.value)

        self.assertIn("getchooser", stage_labels)

    def test_every_bad_file_record_retains_its_source_path(self) -> None:
        tree = ast.parse((ANALYSIS_DIR / "get_continuum_params.py").read_text())
        records = []
        for node in ast.walk(tree):
            if not (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "badfiles"
                and node.func.attr == "append"
            ):
                continue
            self.assertEqual(len(node.args), 1)
            self.assertIsInstance(node.args[0], ast.List)
            records.append(node.args[0])

        self.assertTrue(records)
        for record in records:
            self.assertEqual(
                len(record.elts),
                2,
                f"bad-file record at line {record.lineno} does not include a path",
            )
            self.assertIsInstance(record.elts[1], ast.Name)
            self.assertEqual(record.elts[1].id, "file")

    def test_source_coordinate_state_is_initialized_and_fails_closed(self) -> None:
        fitter = (ANALYSIS_DIR / "get_continuum_params.py").read_text()

        self.assertNotIn("iibest=='NA'", fitter)
        self.assertIn("coordinate_record_count=0", fitter)
        self.assertIn("if coordinate_record_count!=1:", fitter)
        self.assertIn("iibest=None", fitter)
        self.assertIn("if i!=1 or iibest is None:", fitter)
        self.assertIn("Expected exactly one assembled coordinate record", fitter)

    def test_unselected_chosen_source_with_multiple_plots_is_removed_once(self) -> None:
        cleanup = self._load_fitter_function(
            "remove_unselected_chosen_frequency_source_folders"
        )
        with tempfile.TemporaryDirectory() as temporary:
            chosen = Path(temporary) / "chosenfreq"
            rejected = chosen / "Rejected"
            retained = chosen / "Retained"
            rejected.mkdir(parents=True)
            retained.mkdir()
            (rejected / "first.pdf").write_bytes(b"first")
            (rejected / "second.pdf").write_bytes(b"second")
            (rejected / "candidate.fit").mkdir()
            (retained / "candidate.fit").mkdir()
            (chosen / "summary.pdf").write_bytes(b"summary")

            removed = cleanup(str(chosen))

            self.assertEqual(removed, [str(rejected)])
            self.assertFalse(rejected.exists())
            self.assertTrue(retained.is_dir())
            self.assertTrue((chosen / "summary.pdf").is_file())

    def test_edited_notebooks_remain_valid_json(self) -> None:
        for name in (
            "../../BH_Parameters.ipynb",
            "post_processing_finalcontinuum.ipynb",
            "Variability_Search.ipynb",
            "get_contams.ipynb",
            "get_urls.ipynb",
        ):
            json.loads((ANALYSIS_DIR / name).resolve().read_text())


if __name__ == "__main__":
    unittest.main()
