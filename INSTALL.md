# Installing the fitting environment

These commands describe the tested **CPython 3.12 / Linux x86_64** route.
The clean installation and runtime-import probe were exercised with Python
3.12.3. They do not establish support for other Python versions or platforms.
Use an independent working copy of the source and new environment/runtime
directories. Do not install into the author's preserved research environment.

## 1. Install the pinned packages

Run from the extracted release root. Replace both example paths with your own
new absolute paths; the virtual environment must not already exist.

```bash
python3.12 -m venv /absolute/path/to/new-wmagn-venv

env -u PYTHONPATH PYTHONNOUSERSITE=1 \
  /absolute/path/to/new-wmagn-venv/bin/python -m pip install \
  -c constraints-fitting.txt \
  -r requirements-casa-linux-cp312.txt \
  -r requirements-mir.txt

/absolute/path/to/new-wmagn-venv/bin/python -m pip check
```

The full installation shown above was tested. `requirements.txt` pins the
analysis packages, including `reproject`. `constraints-fitting.txt` freezes
the resolved fitting stack. All 152 dependency pins match the installed
validation environment. The installation used pip 25.1.1; constraints do not
install or bootstrap pip itself. The unused `wheel` distribution was not
installed and is not a fitting requirement. The CASA requirements select the official
CASA 6.7.6.14 CPython-3.12 `manylinux_2_28` wheel by URL and SHA256; the
default newer-platform wheel failed an OpenSSL compatibility check on the
validation host. No system-library replacement is needed for this route.

`requirements-mir.txt` adds the optional `python-casacore` MS writer through
the PyUVData `casa` extra. It is needed for MIR-to-MS conversion, not for the
supported route that starts with a delivered CASA Measurement Set. It may
be omitted if that optional conversion will not be used; the combined command
above is the full tested package installation.

## 2. Create the isolated CASA runtime

CASA also needs reference catalogs and `analysisUtils`, separately from the
Python packages. Choose a new runtime directory whose parent already exists:

```bash
env -u PYTHONPATH \
  /absolute/path/to/new-wmagn-venv/bin/python -B scripts/setup_runtime.py \
  --runtime /absolute/path/to/new-wmagn-runtime \
  --pull-data
```

This downloads the pinned official `analysisUtils` 2.398 archive and verifies
its hash before safe extraction. An already downloaded, identical archive may
be supplied with `--analysis-utils-archive /path/to/analysis_scripts.tar`.
The `--pull-data` flag runs the official `casaconfig` reference-data download
in the chosen runtime, with its output captured in `logs/pull_data.log`.
It records the downloaded data inventory; subsequent automatic updates are
disabled. A setup without `--pull-data` does not supply fitting-ready catalogs.

The setup refuses existing directories, including partial failed setups.
Keep a failed directory for diagnosis and use a different new path for another
attempt. The runtime's config, CASA data, logs, caches and plotting/Jupyter
settings are isolated from existing user settings; `HOME` is not overridden.
`runtime_manifest.json` records the selected interpreter, dependency bundle,
configuration hashes and reference-data inventory. Moving or replacing the
virtual environment requires a new runtime setup.

## 3. Verify actual imports and a Measures lookup

Use the generated executable, not a bare Python interpreter, for CASA fitting:

```bash
/absolute/path/to/new-wmagn-runtime/wm-agn-casa-python -c \
  'import numpy,scipy,pandas,matplotlib,astropy,reproject,pyuvdata,analysisUtils,mpfit; from casatools import image,measures; from casatasks import imfit,imhead; from casaconfig import config; m=measures(); print(m.observatory("ALMA")); print(config.measurespath, config.data_auto_update, config.measures_auto_update)'
```

The lookup must return an ALMA position, the data path must be inside the new
runtime, and both update flags must be `False`. The actual clean-runtime
import/lookup probe passed during candidate testing. `analysisUtils` may print
its ordinary non-CASA-shell banner or a warning about its unused optional
`pytz` helper; these messages do not themselves indicate a fitter-import
failure. They are not justification for ignoring a traceback or failed lookup.

The launcher accepts a script plus arguments, `-c CODE`, or `-m MODULE`.
Use it as `--python-bin` for `run_continuum_four_pass.py` and set
`--analysis-utils-dir` to:

```text
/absolute/path/to/new-wmagn-runtime/analysis_utils/analysis_scripts
```

The four-pass runner freezes that full bundle. The launcher uses the frozen
copy only if its complete file inventory matches the pinned bundle. An
explicit `WM_AGN_ANALYSIS_UTILS_DIR` override must likewise match it; unrelated
ambient `PYTHONPATH` code is not imported.

## 4. Retrieve and validate the observations

Installation does not provide the external telescope data. Follow
[DATA_ACCESS.md](DATA_ACCESS.md) and the four-pass workflow documentation.
The frozen discovery manifest is distinct from the final 501-row catalog:
it requires 544 discovery inputs plus six selected-clean support trees.

```bash
/absolute/path/to/new-wmagn-venv/bin/python \
  analysis/multi_freq_from_archive/frozen_fitter_inputs.py check-manifest

/absolute/path/to/new-wmagn-venv/bin/python \
  analysis/multi_freq_from_archive/frozen_fitter_inputs.py validate \
  --archive-root /absolute/path/to/recovered-project
```

The public recovery command can download direct ALMA/VLA files and rebuild
the five reviewed Stokes-I products with exact hash checks. It does not supply
the historical SMA image trees. A fresh SMA archive reduction is a separately
identified route: an independent recent delivery gave different fluxes,
including the unresolved NGC 4258 result. Do not claim that installing the
software resolves that difference or guarantees exact fresh-SMA reproduction.

Optional MIR conversion uses the ordinary environment Python, not the CASA
launcher. The input must already have the archive's requested averaging.
Run it only in a new working copy; existing outputs and partial failed MSs
are not overwritten. The converter's path/file-safety tests do not constitute
a completed fresh-MIR conversion test.

## Figure rendering is a separate environment

The fitting stack above uses Matplotlib 3.11.2 and NumPy 2.3.5. It is not the
environment used to verify exact publication-figure bytes. That rendering
stack used Ubuntu 24.04, Python 3.12.3, Matplotlib 3.6.3, NumPy 1.26.4,
SciPy 1.11.4, Pillow 10.2.0, system FreeType 2.13.2 and the Agg backend.
Fonts, rendering configuration and empty-cache tests also matter. Follow
[RENDERING_ENVIRONMENT.txt](RENDERING_ENVIRONMENT.txt) and
`scripts/render_publication.py` for that separate route; do not replace CASA fitting dependencies to
make a figure hash pass or infer a scientific numerical error from a graphics
library difference alone.

## Run the complete source-test inventory

With the fitting runtime and independently prepared locked Ubuntu renderer
rootfs available, run from the release root:

```bash
python3 scripts/test_release.py \
  --fitting-python /absolute/path/to/new-wmagn-runtime/wm-agn-casa-python \
  --renderer-rootfs /absolute/path/to/locked-ubuntu2404-rootfs \
  --output /absolute/path/to/new-external-test-output \
  --release .
```

The output directory must not already exist and must be outside the release
and renderer rootfs. This command inventories both public test suites, runs
scientific/table/workflow tests with the supplied fitting Python, and routes
the two figure-artifact test modules through the locked renderer. A nonzero
invocation, missing required module, or source-tree change makes the command
fail; inspect the generated logs, test counts and skips rather than treating
an invocation alone as a pass. It does not download observations or substitute
for actually running the four-pass fitter and downstream regression.

Successful installation/imports are only the first checkpoint. The original
`CANDIDATE_TEST_REPORT.txt` records the rc1 actual fitting, workflow gates,
downstream numerical comparison and figure/table verification inherited by
v1.0.11. `STABLE_RELEASE_PROMOTION.txt` records the separate metadata-only
promotion checks; it does not claim a repeated scientific run.
