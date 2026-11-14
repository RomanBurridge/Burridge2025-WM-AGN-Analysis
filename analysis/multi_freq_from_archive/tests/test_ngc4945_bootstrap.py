"""Metadata-only preparation/transition regression; no CASA or fitted replay."""
import csv
import io
import json
from pathlib import Path
from types import SimpleNamespace
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

PACKAGE = Path(__file__).resolve().parents[1]
if str(PACKAGE) not in sys.path:
    sys.path.insert(0, str(PACKAGE))

import ngc4945_bootstrap as bootstrap
import run_continuum_four_pass as runner
from rerun_ngc4945_nuclear import (
    merge_correction_manifest, preserve_pre_component_fields, write_csv,
    PHOTOMETRY_FIELDS,
)
try:
    from .frozen_input_fixture import fixture_load, fixture_select, write_fixture_manifest
except ImportError:  # direct unittest discovery adds tests/ rather than package root
    from frozen_input_fixture import fixture_load, fixture_select, write_fixture_manifest

PUBLIC_TEST_FIXTURES = Path(__file__).resolve().parent/'fixtures/source18_public_v1_0_10'


def authenticated_public_test_fixture(name, root=PUBLIC_TEST_FIXTURES):
    """Original PUBLIC test data only; never a current final-state assumption."""
    if name not in (bootstrap.PHOTOMETRY, bootstrap.CORRECTIONS):
        raise ValueError('unsupported original-public test fixture')
    payload=(root/name).read_bytes()
    if bootstrap.sha256(payload)!=bootstrap.PUBLIC_SOURCE_SHA256[name]:
        raise ValueError('original-public test fixture SHA256 mismatch: '+name)
    return payload


class Source18BootstrapTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.root = self.base / "analysis"
        self.root.mkdir()
        shutil.copytree(PACKAGE/bootstrap.INITIAL_HISTORY_DIR,
                        self.root/bootstrap.INITIAL_HISTORY_DIR)
        self.original = {}
        for name in (*bootstrap.PUBLIC_SOURCE_SHA256, bootstrap.POLICY):
            payload = (PACKAGE / name).read_bytes()
            self.original[name] = payload
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)

    def initialize(self):
        payloads, evidence = bootstrap.derive_bootstrap(self.root)
        return bootstrap.initialize_bootstrap(self.root, payloads, evidence)

    def original_public_operation(self, name='original-public-operation'):
        public=self.base/name
        shutil.copytree(self.root,public)
        # ONLY this release-build operation requires original45/84 bytes. All
        # normal preparation tests continue to use CURRENT final metadata.
        for filename in (bootstrap.PHOTOMETRY,bootstrap.CORRECTIONS):
            (public/filename).write_bytes(authenticated_public_test_fixture(filename))
        return bootstrap.derive_initial_history_from_public(public)

    def test_public_derivation_matches_authenticated_historical_bytes(self):
        payloads = self.original_public_operation()
        prepared, evidence = bootstrap.derive_bootstrap(self.root)
        self.assertEqual(prepared, payloads)
        self.assertEqual({n: bootstrap.sha256(b) for n, b in payloads.items()},
                         bootstrap.HISTORICAL_SHA256)
        self.assertEqual(evidence['counts']['historical_photometry'], 27)
        self.assertEqual(evidence['counts']['historical_corrections'], 66)
        self.assertEqual({n: (self.root/n).read_bytes() for n in self.original}, self.original)

    def test_original_public_operation_survives_current_diagnostic_promotion(self):
        _,_,records=bootstrap.csv_records(self.original[bootstrap.PHOTOMETRY])
        rows=[dict(row) for row,_ in records]
        rows[0]['noise_rms_mjy_per_beam']=str(float(rows[0]['noise_rms_mjy_per_beam'])+1e-13)
        write_csv(self.root/bootstrap.PHOTOMETRY,PHOTOMETRY_FIELDS,rows)
        promoted=(self.root/bootstrap.PHOTOMETRY).read_bytes()
        self.assertNotEqual(promoted,self.original[bootstrap.PHOTOMETRY])
        expected=self.original_public_operation()
        actual,evidence=bootstrap.derive_bootstrap(self.root)
        self.assertEqual(actual,expected)
        self.assertEqual(evidence['release_metadata_sha256'][bootstrap.PHOTOMETRY],bootstrap.sha256(promoted))
        self.assertEqual((self.root/bootstrap.PHOTOMETRY).read_bytes(),promoted)

    def test_original_public_fixture_tamper_is_not_skipped_or_relaxed(self):
        fixture=self.base/'tampered-public-fixture'; fixture.mkdir()
        for filename in (bootstrap.PHOTOMETRY,bootstrap.CORRECTIONS):
            payload=authenticated_public_test_fixture(filename)
            (fixture/filename).write_bytes(payload+b'\n')
            with self.assertRaisesRegex(ValueError,'fixture SHA256 mismatch'):
                authenticated_public_test_fixture(filename,fixture)

    def test_csv_preserves_multiline_quoting_crlf_and_missing_not_zero(self):
        data = b'a,b\r\n"one\r\ntwo","NA"\r\n"",0\r\n'
        fields, header, rows = bootstrap.csv_records(data)
        self.assertEqual(fields, ('a', 'b'))
        self.assertEqual(rows[0][0]['b'], 'NA')
        self.assertEqual(rows[1][0], {'a': '', 'b': '0'})
        self.assertEqual(header + b''.join(raw for _, raw in rows), data)

    def test_bad_csv_grammars_fail(self):
        for raw in (b'a,a\nx,y\n', b'a,b\n\nx,y\n', b'a,b\nx,y,z\n', b''):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                bootstrap.csv_records(raw)

    def test_duplicate_products_fail(self):
        with self.assertRaisesRegex(ValueError, 'duplicate product'):
            bootstrap.unique_products([({'product': 'ALMA/a.fits'}, b'')]*2, 'test')

    def test_corrupted_final_sources_fail_before_snapshot_writes(self):
        for name in bootstrap.PUBLIC_SOURCE_SHA256:
            with self.subTest(name=name):
                path = self.root/name
                path.write_bytes(self.original[name] + b'\n')
                with self.assertRaises(ValueError):
                    bootstrap.derive_bootstrap(self.root)
                self.assertFalse((self.root/bootstrap.SNAPSHOT_DIR).exists())
                path.write_bytes(self.original[name])

    def test_duplicate_or_wrong_disposition_admission_fails(self):
        fields, header, rows = bootstrap.csv_records(self.original[bootstrap.POLICY])
        selected = next((r, raw) for r, raw in rows if r['source'] == 'NGC4945')
        (self.root/bootstrap.POLICY).write_bytes(self.original[bootstrap.POLICY] + selected[1])
        with self.assertRaisesRegex(ValueError, 'duplicate product'):
            bootstrap.derive_bootstrap(self.root)
        values = [dict(r) for r, _ in rows]
        next(r for r in values if r['source'] == 'NGC4945')['expected_run_disposition'] = 'rejected'
        buf = io.StringIO(newline='')
        writer = csv.DictWriter(buf, fields, lineterminator='\n')
        writer.writeheader(); writer.writerows(values)
        (self.root/bootstrap.POLICY).write_text(buf.getvalue())
        with self.assertRaisesRegex(ValueError, 'disposition'):
            bootstrap.derive_bootstrap(self.root)

    def test_initializer_refuses_reuse_and_wrong_destination(self):
        payloads, evidence = bootstrap.derive_bootstrap(self.root)
        (self.root/bootstrap.PHOTOMETRY).write_bytes(b'wrong destination')
        with self.assertRaisesRegex(ValueError, 'unchanged release metadata'):
            bootstrap.initialize_bootstrap(self.root, payloads, evidence)
        self.assertFalse((self.root/bootstrap.SNAPSHOT_DIR).exists())
        (self.root/bootstrap.PHOTOMETRY).write_bytes(self.original[bootstrap.PHOTOMETRY])
        bootstrap.initialize_bootstrap(self.root, payloads, evidence)
        with self.assertRaises(ValueError):
            bootstrap.initialize_bootstrap(self.root, payloads, evidence)

    def test_bad_initialization_evidence_fails_before_any_write(self):
        payloads, evidence = bootstrap.derive_bootstrap(self.root)
        evidence['status'] = 'completed'
        with self.assertRaisesRegex(ValueError, 'initialization evidence mismatch'):
            bootstrap.initialize_bootstrap(self.root, payloads, evidence)
        self.assertFalse((self.root/bootstrap.SNAPSHOT_DIR).exists())
        self.assertEqual((self.root/bootstrap.PHOTOMETRY).read_bytes(), self.original[bootstrap.PHOTOMETRY])

    def test_wrong_historical_payload_fails_before_any_write(self):
        payloads, evidence = bootstrap.derive_bootstrap(self.root)
        payloads[bootstrap.PHOTOMETRY] += b'\n'
        with self.assertRaisesRegex(ValueError, 'initialization payload mismatch'):
            bootstrap.initialize_bootstrap(self.root, payloads, evidence)
        self.assertFalse((self.root/bootstrap.SNAPSHOT_DIR).exists())

    def test_preapplication_refuses_final_state_instead_of_weakening45_guard(self):
        record = self.initialize()
        bootstrap.validate_preapplication_lineage(self.root, record)
        (self.root/bootstrap.PHOTOMETRY).write_bytes(self.original[bootstrap.PHOTOMETRY])
        with self.assertRaisesRegex(ValueError, 'not authenticated27/66'):
            bootstrap.validate_preapplication_lineage(self.root, record)
        rows = list(csv.DictReader(io.StringIO(self.original[bootstrap.PHOTOMETRY].decode())))
        with self.assertRaisesRegex(RuntimeError, 'exactly 27'):
            preserve_pre_component_fields(rows, self.root/bootstrap.PHOTOMETRY,
                                          self.root/bootstrap.CORRECTIONS)

    def test_snapshot_and_preapplication_corruption_fail(self):
        record = self.initialize()
        for name in bootstrap.HISTORICAL_SHA256:
            path = self.root/bootstrap.SNAPSHOT_DIR/name
            before = path.read_bytes(); path.write_bytes(before+b'\n')
            with self.assertRaisesRegex(ValueError, 'immutable snapshot changed'):
                bootstrap.validate_bootstrap_snapshot(self.root, record)
            path.write_bytes(before)
        path = self.root/bootstrap.CORRECTIONS
        path.write_bytes(path.read_bytes()+b'\n')
        with self.assertRaisesRegex(ValueError, 'not authenticated27/66'):
            bootstrap.validate_preapplication_lineage(self.root, record)

    def test_symlinks_rejected(self):
        path = self.root/bootstrap.PHOTOMETRY
        other = self.base/'photometry'
        path.rename(other); path.symlink_to(other)
        with self.assertRaisesRegex(ValueError, 'symlink'):
            bootstrap.derive_bootstrap(self.root)

    def test_initial_history_asset_and_provenance_tampering_fail(self):
        for name in (*bootstrap.HISTORICAL_SHA256, bootstrap.PROVENANCE_NAME):
            path = self.root/bootstrap.INITIAL_HISTORY_DIR/name
            before = path.read_bytes()
            path.write_bytes(before+b'\n')
            with self.assertRaisesRegex(ValueError, 'initial-history'):
                bootstrap.derive_bootstrap(self.root)
            self.assertFalse((self.root/bootstrap.SNAPSHOT_DIR).exists())
            path.write_bytes(before)

    def exercise_transition(self, root, record):
        """Simulate build_records metadata, not a measurement or CASA test."""
        before = self.base/'preapplication.csv'
        shutil.copyfile(root/bootstrap.CORRECTIONS, before)
        rows = list(csv.DictReader(io.StringIO(self.original[bootstrap.PHOTOMETRY].decode())))
        historical = {r['product']: r for r in csv.DictReader(
            io.StringIO((root/bootstrap.PHOTOMETRY).read_text()))}
        # Distinct synthetic fresh fields prove initializer does not replay NEW
        # values. Existing27 OLD fields must survive; new18 keep generic OLD.
        for index, row in enumerate(rows):
            row.update(old_flux_field='123.45', old_snr='7.890',
                       new_flux_field=f'{index+1}.00', new_snr=f'{index+10}.000',
                       unrounded_flux_mjy=f'{index+1}.0')
        preserve_pre_component_fields(rows, root/bootstrap.PHOTOMETRY, root/bootstrap.CORRECTIONS)
        for row in rows:
            if row['product'] in historical:
                self.assertEqual(row['old_flux_field'], historical[row['product']]['old_flux_field'])
                self.assertEqual(row['old_snr'], historical[row['product']]['old_snr'])
            else:
                self.assertEqual((row['old_flux_field'],row['old_snr']), ('123.45','7.890'))
        merge_correction_manifest(root/bootstrap.CORRECTIONS, rows)
        write_csv(root/bootstrap.PHOTOMETRY, PHOTOMETRY_FIELDS, rows)
        report = runner.validate_correction_manifest_transition(
            before, root/bootstrap.CORRECTIONS, root/bootstrap.PHOTOMETRY)
        self.assertEqual((report['preapplication_rows'],report['postrun_rows']), (66,84))
        self.assertEqual(report['unchanged_non_ngc4945_rows'],39)
        self.assertEqual(report['added_ngc4945_rows'],18)
        bootstrap.validate_bootstrap_snapshot(root, record)  # immutable after stages
        with self.assertRaisesRegex(ValueError, 'not authenticated27/66'):
            bootstrap.validate_preapplication_lineage(root, record)

    def test_prepare_then_original_source18_metadata_transition(self):
        release = self.base/'release'
        analysis = release/'analysis/multi_freq_from_archive'
        shutil.copytree(self.root, analysis)
        (analysis/'get_continuum_params.py').write_text('# isolated synthetic archive fixture\n')
        manifest = release/'analysis/sma_calibration/selected_clean_iterations.csv'
        manifest.parent.mkdir(parents=True)
        manifest.write_text('date_code,flagnum,source,image_path,residual_path\n'
            '230316,2,Test,analysis/sma_calibration/selected.image,analysis/sma_calibration/selected.residual\n')
        (analysis/'continuum_boundary_rejection_policy.csv').write_text('synthetic\n')
        archive = self.base/'archive'
        for telescope in ('ALMA','VLA'):
            (archive/'analysis/multi_freq_from_archive'/telescope).mkdir(parents=True)
        discovery = Path('analysis/sma_calibration/discovery')
        (archive/discovery).mkdir(parents=True)
        for name in ('selected.image','selected.residual'):
            (archive/'analysis/sma_calibration'/name).mkdir()
            (archive/'analysis/sma_calibration'/name/'table.dat').write_bytes(b'fixture only')
        write_fixture_manifest(analysis, archive, discovery,
            [Path('analysis/sma_calibration')/name for name in ('selected.image','selected.residual')])
        workspace = self.base/'prepared'
        # Only unrelated archive-admission/copy volumes are synthetic. Public
        # source18 derivation, prepare initialization, snapshot validation and
        # original preserve/merge/transition APIs all execute without mocks.
        with patch.object(runner, 'EXPECTED_ARCHIVE_INPUTS', {'ALMA':0,'VLA':0}), \
             patch.object(runner, 'EXPECTED_SMA_DISCOVERY_INPUTS', 0), \
             patch.object(runner, 'load_reviewed_product_additions', return_value=[]), \
             patch.object(runner, 'load_boundary_rejection_policy', return_value=[]), \
             patch.object(runner, 'load_frozen_fitter_inputs', new=fixture_load), \
             patch.object(runner, 'validate_and_select_frozen_inputs', new=fixture_select):
            runner.prepare_workspace(release=release, archive_project_root=archive,
                workspace=workspace,date='230316',flagnum=2,sma_discovery_relative=discovery)
            inventory = runner.validate_prepared_archive_inputs(workspace)
            # Legitimately promoted final photometry may have different unused
            # RMS/sigma diagnostics. A NEW prepare must still use fixed history.
            _, _, promoted_records = bootstrap.csv_records(self.original[bootstrap.PHOTOMETRY])
            promoted = [dict(r) for r, _ in promoted_records]
            promoted[0]['noise_rms_mjy_per_beam'] = str(float(promoted[0]['noise_rms_mjy_per_beam'])+1e-13)
            write_csv(analysis/bootstrap.PHOTOMETRY, PHOTOMETRY_FIELDS, promoted)
            changed = (analysis/bootstrap.PHOTOMETRY).read_bytes()
            self.assertNotEqual(changed, self.original[bootstrap.PHOTOMETRY])
            second_workspace = self.base/'prepared-after-promotion'
            runner.prepare_workspace(release=release,archive_project_root=archive,
                workspace=second_workspace,date='230316',flagnum=2,sma_discovery_relative=discovery)
            second = runner.analysis_root(second_workspace)
            self.assertEqual(bootstrap.sha256((second/bootstrap.PHOTOMETRY).read_bytes()),
                             bootstrap.HISTORICAL_SHA256[bootstrap.PHOTOMETRY])
            self.assertEqual((analysis/bootstrap.PHOTOMETRY).read_bytes(), changed)
        root = runner.analysis_root(workspace)
        record = inventory['source18_bootstrap']
        self.assertEqual(record, json.loads((workspace/'.continuum-fit-inputs.json').read_text())['source18_bootstrap'])
        self.assertEqual((analysis/bootstrap.PHOTOMETRY).read_bytes(), changed)
        self.assertEqual((analysis/bootstrap.CORRECTIONS).read_bytes(), self.original[bootstrap.CORRECTIONS])
        self.exercise_transition(root, record)
        # A measured/final-state workspace must fail BEFORE any new preflight
        # or subprocess, rather than resetting/replaying the prior correction.
        with patch.object(runner, 'workflow_input_hashes', return_value={}), \
             patch.object(runner, 'validate_prepared_archive_inputs', return_value=inventory), \
             patch.object(runner, 'preflight_python') as preflight, \
             patch.object(runner.subprocess, 'run') as subprocess_run:
            with self.assertRaisesRegex(ValueError, 'not authenticated27/66'):
                runner.run_workflow(SimpleNamespace(workspace=workspace), dry_run=True)
            preflight.assert_not_called()
            subprocess_run.assert_not_called()


if __name__ == '__main__':
    unittest.main()
