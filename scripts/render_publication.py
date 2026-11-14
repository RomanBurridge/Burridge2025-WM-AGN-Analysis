#!/usr/bin/env python3
"""Render publication figures in an explicitly supplied read-only Linux rootfs.

The external runtime is specified in RENDERING_ENVIRONMENT.txt. This program
does not install a host OS, replace scientific inputs or overwrite prior runs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def command_prefix(rootfs: Path, release: Path, output: Path) -> list[str]:
    return ['bwrap', '--die-with-parent', '--unshare-all', '--new-session',
            '--ro-bind', str(rootfs), '/', '--proc', '/proc', '--dev', '/dev',
            '--tmpfs', '/tmp', '--tmpfs', '/mnt', '--dir', '/mnt/release',
            '--dir', '/mnt/output', '--ro-bind', str(release), '/mnt/release',
            '--bind', str(output), '/mnt/output', '--clearenv',
            '--setenv', 'PATH', '/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin',
            '--setenv', 'LC_ALL', 'C.UTF-8', '--setenv', 'TZ', 'UTC',
            '--setenv', 'PYTHONNOUSERSITE', '1', '--setenv', 'PYTHONDONTWRITEBYTECODE', '1',
            '--setenv', 'XDG_CACHE_HOME', '/mnt/output/runtime/cache',
            '--setenv', 'XDG_CONFIG_HOME', '/mnt/output/runtime/config',
            '--chdir', '/mnt/release', '/usr/bin/python3', '-B']


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rootfs', type=Path, required=True)
    parser.add_argument('--release', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rootfs, release, output = (p.resolve() for p in (args.rootfs, args.release, args.output))
    if shutil.which('bwrap') is None:
        parser.error('bubblewrap (bwrap) is required for the isolated renderer')
    if not (rootfs / 'usr/bin/python3').exists() or not (release / 'analysis').is_dir():
        parser.error('Invalid renderer rootfs or release directory')
    if output.exists() or args.output.is_symlink():
        parser.error('Output already exists; preserve it and choose a new run directory')
    if output.is_relative_to(release) or release.is_relative_to(output):
        parser.error('Use a new output directory separate from the source release')
    output.mkdir(parents=True)
    for relative in ('runtime/cache', 'runtime/config', 'logs', 'flux', 'SEDs'):
        (output / relative).mkdir(parents=True, exist_ok=True)
    prefix = command_prefix(rootfs, release, output)
    probe = ('import json,sys,matplotlib,numpy,scipy,PIL; '
             'from matplotlib import ft2font; '
             'v={"python":".".join(map(str,sys.version_info[:3])), '
             '"matplotlib":matplotlib.__version__,"numpy":numpy.__version__, '
             '"scipy":scipy.__version__,"pillow":PIL.__version__, '
             '"freetype":ft2font.__freetype_version__}; '
             'print(json.dumps(v,sort_keys=True)); '
             'assert v=={"python":"3.12.3","matplotlib":"3.6.3", '
             '"numpy":"1.26.4","scipy":"1.11.4","pillow":"10.2.0", '
             '"freetype":"2.13.2"},"Renderer version mismatch"')
    runtime = subprocess.run([*prefix, '-c', probe], check=True, text=True,
                             capture_output=True)
    (output / 'renderer_versions.json').write_text(runtime.stdout, encoding='utf-8')
    analysis = 'analysis/multi_freq_from_archive/'
    inputs = [release / 'analysis/CountingBHS_and_Spectra/BHS_params.csv',
              release / (analysis + 'machinetables/S230_extrapolation_comparison.csv'),
              release / (analysis + 'fitsumfiles/final/fitsummary_final_withextras.txt')]
    # The released final table's canonical name is resolved explicitly.
    if not inputs[2].is_file():
        raise FileNotFoundError(inputs[2])
    before = {str(p.relative_to(release)): digest(p) for p in inputs}
    commands = [
        ('bhs', [analysis + 'build_bhs_flux_figure.py', '--output-dir', '/mnt/output/flux']),
        ('sed', [analysis + 'build_sed_figures.py', '--output-dir', '/mnt/output/SEDs']),
    ]
    for label, argv in commands:
        with (output / 'logs' / f'{label}.stdout.log').open('wb') as stdout, \
             (output / 'logs' / f'{label}.stderr.log').open('wb') as stderr:
            subprocess.run([*prefix, '/mnt/release/' + argv[0], *argv[1:]],
                           check=True, stdout=stdout, stderr=stderr)
    files = sorted([*(output / 'flux').glob('*.png'), *(output / 'flux').glob('*.pdf'),
                    *(output / 'SEDs').glob('*.png')])
    if len(files) != 8:
        raise ValueError(f'Expected seven PNGs and one PDF, got {len(files)} artifacts')
    if before != {str(p.relative_to(release)): digest(p) for p in inputs}:
        raise ValueError('Scientific inputs changed while rendering')
    report = {'status': 'rendered_not_yet_visually_reviewed', 'scientific_inputs': before,
              'artifacts': {str(p.relative_to(output)): digest(p) for p in files},
              'runtime': json.loads(runtime.stdout)}
    (output / 'RENDER_REPORT.json').write_text(json.dumps(report, indent=2) + '\n',
                                            encoding='utf-8')
    print(json.dumps(report, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
