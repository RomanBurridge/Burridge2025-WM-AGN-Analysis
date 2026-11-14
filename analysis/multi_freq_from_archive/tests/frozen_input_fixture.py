"""Synthetic archive-volume adapter; production frozen-manifest tests are separate."""
import csv
from pathlib import Path
from frozen_fitter_inputs import FrozenInputSelection, sha256_file, sha256_workflow_tree


def write_fixture_manifest(analysis, archive, discovery, selected):
    rows = []
    archive_analysis = archive/'analysis/multi_freq_from_archive'
    for telescope in ('ALMA', 'VLA'):
        for path in sorted((archive_analysis/telescope).rglob('*.fits')):
            rows.append(dict(input_role='archive_fitter_input',
                relative_path=path.relative_to(archive_analysis).as_posix(),
                expected_sha256=sha256_file(path)))
    for path in sorted((archive/discovery).glob('*_test.image')):
        rows.append(dict(input_role='sma_discovery_input',
            relative_path=path.relative_to(archive).as_posix(),
            expected_sha256=sha256_workflow_tree(path)))
    for relative in selected:
        rows.append(dict(input_role='sma_selected_clean_support',
            relative_path=relative.as_posix(),
            expected_sha256=sha256_workflow_tree(archive/relative)))
    path = analysis/'frozen_fitter_inputs.csv'
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, ('input_role','relative_path','expected_sha256'))
        writer.writeheader(); writer.writerows(rows)


def fixture_load(manifest):
    with Path(manifest).open(newline='') as stream:
        return list(csv.DictReader(stream))


def fixture_select(archive_root, manifest, **_kwargs):
    root = Path(archive_root)
    rows = fixture_load(manifest)
    paths = [(r, root/('analysis/multi_freq_from_archive' if
              r['input_role'] == 'archive_fitter_input' else '')/r['relative_path'])
             for r in rows]
    return FrozenInputSelection(sha256_file(manifest), root,
        tuple(p for r,p in paths if r['input_role'] != 'sma_selected_clean_support'),
        tuple(p for r,p in paths if r['input_role'] == 'sma_selected_clean_support'))
