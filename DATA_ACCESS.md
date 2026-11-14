# External Data Access

This Zenodo release intentionally excludes archival data files and large or intermediate generated CASA image products. The lightweight publication figures are included, while the external observatory data can be retrieved from the original archives.

## Canonical recovery manifest

Start with
`analysis/multi_freq_from_archive/machinetables/continuum_recovery_manifest.csv`.
The final-release manifest has one row for each of the 501 canonical
continuum measurements and records the exact release-relative product path or
literature citation, a stable record ID, the observatory key, and the
authoritative retrieval service.  Its 485 locally fitted products comprise 222
ALMA FITS files, 248 NVAS images, and 15 SMA CASA image directories.  They have
recorded byte sizes and SHA-256 digests totaling 6,040,874,990 payload bytes;
the other 16 rows are explicitly
marked `literature_only_no_local_product` and link to their cited papers.
These counts were regenerated from the exactly adjudicated authoritative Task
28 run.

The audit-time integrity values are locked in
`analysis/multi_freq_from_archive/continuum_product_checksums.csv`.  A normal
structural check is lightweight and does not need the external products:

```bash
python3 analysis/multi_freq_from_archive/build_continuum_recovery_manifest.py --check
```

For an augmented local project tree laid out as described below, verify the
exact files and CASA image directories against every recorded checksum:

```bash
python3 analysis/multi_freq_from_archive/build_continuum_recovery_manifest.py \
  --check --archive-root /path/to/augmented/Burridge2025-WM-AGN-Analysis
```

Regular products use ordinary file SHA-256.  A CASA image is a directory, so
`sha256-relpath-tree-v1` hashes its sorted relative directory names, symlink
targets, and regular-file sizes and contents.  Absolute locations, timestamps,
owners, permissions, and inode numbers are excluded.  The manifest's
directory `local_size_bytes` is the sum of its regular-file payload sizes;
directory entries and link-target text do not add to that size.  The manifest's
`local_product_hash_verified` status records the audit-time verification; it
does not imply that the large product is bundled in this Zenodo release.

## Frozen fitter-input membership

The 501-row recovery manifest describes the final published measurements; it
is not the complete discovery-pass input list.  The exact Task 28 fitter
membership is frozen separately in
`analysis/multi_freq_from_archive/frozen_fitter_inputs.csv`.  That public-safe
file contains 550 records: 529 archive FITS files (222 ALMA and the exact 307
historical NVAS products), 15 SMA dirty-image discovery trees, and six
selected-clean SMA image/residual support trees.  Only the first 544 records
are discovery inputs.  The six support trees are three image/residual pairs
for NGC 1068, NGC 3079, and NGC 4258.

The frozen file prevents a fresh archive query from silently changing the
historical selection.  In particular, the reproducible public NVAS query now
returns 341 products, but only the hash-locked historical 307 are admitted.
Regular files use ordinary whole-file SHA-256.  The SMA trees use the original
four-pass runner's `sha256-workflow-tree-v1` algorithm, which is intentionally
different from the final-measurement manifest's
`sha256-relpath-tree-v1` algorithm.

Check the lightweight manifest itself with:

```bash
python3 analysis/multi_freq_from_archive/frozen_fitter_inputs.py check-manifest
```

For a recovered project tree, hash all 550 required products and return the
exact 544-input discovery whitelist plus the six selected-clean paths with:

```bash
python3 analysis/multi_freq_from_archive/frozen_fitter_inputs.py validate \
  --archive-root /path/to/augmented/Burridge2025-WM-AGN-Analysis
```

The optional recovery command downloads the 524 directly retrievable archive
files and the five public IQUV parent products.  Add `--derive-stokes-i` to
load the reviewed deterministic Task 28 extractor, rebuild the five historical
Stokes-I products through temporary files, and require their exact frozen
whole-file hashes before installation:

```bash
python3 analysis/multi_freq_from_archive/frozen_fitter_inputs.py recover-public \
  --archive-root /path/to/augmented/Burridge2025-WM-AGN-Analysis \
  --derive-stokes-i
```

Without that flag, the command only caches and hash-verifies the five parents.
Existing matching outputs are reused; nonmatching outputs are never
overwritten.  The recovery rows point to the real parent products rather than
to nonexistent derived-output URLs.  The 21 SMA products remain a separate
archive/calibration step.  A final `validate` call is the required exact-byte
gate before fitting.

## Inclusive ALMA archive census and additions

The first inclusive public-release census through 2025-01-27 found 112 ALMA
Member OUS records absent from the historical input inventory.  Seven products
satisfied the then-current continuum rules and 105 received explicit
exclusions.  Six of the seven were ESO products hidden by a raw `SO` substring
test, and the seventh was an independently omitted NGC 1068 product.

Task 28 then audited the filename-dependent selection rules themselves.  It
approved 49 pre-cutoff products omitted only by filename, 26 straightforward
post-cutoff target continuum datasets, and one Stokes-I measurement from each
of seven MFS/polarization deliveries.  Together with the seven earlier
admissions, `continuum_product_addition_policy.csv` therefore contains 89
products relative to the immutable 409-row historical baseline.  The prepared
four-pass inventory now contains 544 inputs (222 ALMA, 307 VLA, and 15 SMA),
and the exactly reviewed raw publication table contains 495 rows (406 retained
historical-common rows plus 89 admissions).

The shipped discovery code evaluates scientific eligibility from reviewed
product metadata and basename tokens rather than treating `tt0`, `mfs`,
`manual`, or `split.cal` as automatic rejection terms.  Alternate reductions
and per-spectral-window products are not counted as independent full-band
measurements.  The archive inputs remain external: the release records their
retrieval identifiers and hashes but does not bundle their FITS bytes.  The
census, approval manifests, query/provenance records, and preserved pre-fix
materials are included under
`analysis/multi_freq_from_archive/archive_inventory_audit/`; that directory
must contain no FITS products.

## Selected S230 image WCS

`analysis/multi_freq_from_archive/s230_wcs_evidence.csv` is the portable WCS
snapshot for the exact 19 products selected in
`machinetables/S230_extrapolation_comparison.csv`.  It records the matching
recovery ID and product hash plus the image's celestial reference code, raw
FITS `RADESYS`/`EQUINOX`/`CTYPE` when applicable, projection, reference pixel
and value, increments, PC matrix, and transform convention.

With the recovered project tree and CASA Python, re-read the actual products
and compare the snapshot with:

```bash
/path/to/casa-python analysis/multi_freq_from_archive/extract_s230_wcs_evidence.py \
  --archive-root /path/to/augmented/Burridge2025-WM-AGN-Analysis --check
```

The lightweight component/reference-coordinate join can then be checked
without the external images:

```bash
python3 analysis/multi_freq_from_archive/build_s230_astrometry_audit.py --check
```

See `analysis/multi_freq_from_archive/machinetables/s230_anchor_astrometry_schema.md`
for the common-frame transform and evidence boundaries.

## ALMA

ALMA continuum products can be retrieved from the ALMA Science Archive:

- Web archive: https://almascience.org/aq/
- Regional archive mirrors: https://almascience.eso.org/aq/ , https://almascience.nao.ac.jp/aq , and https://almascience.nrao.edu/aq/
- Programmatic access: `astroquery.alma`

The recovery manifest lists the full ALMA product filename and Member OUS UID
for every ALMA row. These identifiers can be used in the ALMA Science Archive
query interface or through `astroquery.alma` to locate the corresponding
products.  The early Cycle-1 product
`uid___A002_X5b1929_X57d.continuum.fits` is explicitly cross-referenced to
project `2012.1.00912.S`, Member OUS `uid://A002/X5a9a13/X528`, and execution
block `uid://A002/X5b1929/X57d`.  Seven additional Task 28 deliveries retain
`uid___A002...` execution-block filenames; their Member OUS and proposal
identities are taken from the reviewed, hash-validated Task 28 input manifest.
ALMA also provides a Request Handler download script for command-line retrieval
after selecting products in the archive.

For the NGC 4945 fixed-nucleus correction,
`analysis/multi_freq_from_archive/ngc4945_nuclear_photometry.csv` lists all 45
required product paths and SHA-256 hashes: 27 historical products plus 18
Task 28 additions.  Arrange the downloaded products
under the same relative `ALMA/NGC4945/...` paths and pass their common parent
to `rerun_ngc4945_nuclear.py --archive-root`; the script transforms the
published source-18 ICRS coordinate through each image WCS.

Minimal programmatic sketch:

```python
from astroquery.alma import Alma

alma = Alma()
result = alma.query_object("NGC 1068", public=True)
print(result.colnames)
```

## VLA / NVAS

VLA data and products can be retrieved through NRAO:

- NRAO Data Archive / Archive Access Tool: https://data.nrao.edu
- NRAO VLA archive documentation: https://science.nrao.edu/facilities/vla/archive/index
- NVAS historical VLA image products: https://www.vla.nrao.edu/astro/nvas/
- Programmatic NVAS image-list access: `astroquery.nvas`

The recovery manifest lists both the NVAS program token and exact image
filename for every fitted historical VLA row. These filenames can be matched
through the NVAS service. Raw and modern VLA data can be located through the
NRAO Archive Access Tool; NRAO also provides TAP-based scripted metadata
queries, while downloads are completed through the archive web/request
workflow. The generated `analysis/multi_freq_from_archive/machinetables/continuum_product_aliases.csv` sidecar preserves both NVAS filenames and SHA-256 hashes for the two reductions of the 1998-01-06 IC 2560 BB091 observation; only the later 2009 `BB0091` reduction is retained as a measurement.

Three NGC 4258 NVAS products are intentionally not used for nuclear photometry:
the 1997-01-04 BR0043 15.36 GHz image and the 1996-12-31 AH594 22.46 GHz
image, plus the 2004-05-22 AH0847 22.46 GHz image.  The first two free fits
land 10.4 and 6.5 restoring-beam major axes from the maser nucleus; bright
diagonal stripe artifacts cross the nucleus in the AH0847 image and make its
legacy 79.89 mJy row unreliable.  `continuum_product_exclusion_policy.csv` and
the generated `machinetables/continuum_product_exclusions.csv` preserve all
three filenames, SHA-256 hashes, FITS-header evidence, diagnostics, and former
table rows.  The continuum table instead identifies the first two
same-observation reprocessed detections from Doi et al. (2013, ApJ 765, 63;
doi:10.1088/0004-637X/765/1/63) as literature measurements; those values do not
claim provenance from the failed NVAS images.  The AH0847 image has no
replacement row.

Minimal NVAS sketch:

```python
from astroquery.nvas import Nvas

urls = Nvas.get_image_list("12h18m57.5s +47d18m14s", radius="0d0m30s")
for url in urls[:5]:
    print(url)
```

## SMA

SMA data can be retrieved from the SMA public archive:

- SMA Data Archive: https://lweb.cfa.harvard.edu/rtdc/SMAdata/archives/
- SMA archive instructions: https://lweb.cfa.harvard.edu/rtdc/SMAdata/instructions/archive_use.html

The SMA observations are associated with programs `2022B-H002` and `2021B-H004`. Search the SMA archive by project code and observing date. The archive provides public data after the proprietary period and can provide raw MIR data or CASA measurement-set formats when available.  All 15 SMA continuum-table rows use the `2022B-H002` track observed on 2023-03-16; the recovery manifest records each selected `nitercal_0` or clean CASA image path and its deterministic directory checksum.

This release includes the observing-script text files used by the calibration scripts:

- `analysis/sma_calibration/txt_230316.txt`
- `analysis/sma_calibration/txt_220131.txt`

Raw SMA interferometry directories contain the science targets and calibration sources observed during the night; keep that structure intact when placing data under `analysis/sma_calibration/sma_data/`.
