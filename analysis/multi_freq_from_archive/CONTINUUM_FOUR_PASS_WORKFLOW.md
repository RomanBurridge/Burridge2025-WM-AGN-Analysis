# Reproducible four-pass continuum workflow

## What was restored

The historical fitter used these modes in order:

1. `prelimchecker`
2. `fit_sizenoloc`
3. `getchoosercoords`
4. `fixtochosenfreq`

The missing `get_chosenfreq.ipynb` made the transition between modes 3 and 4
non-reproducible.  It is replaced by `chosen_frequency_anchors.py`, which
writes a compact candidate manifest, applies a deterministic selection rule,
checks the winner against the reviewed product, reference-frame, coordinate,
and tolerance rules in `chosen_frequency_component_policy.csv`, and copies
only the winning CASA model trees into the generated `images/chosenfreq/`
directory.  That compact nine-source policy is derived from the broader
19-source `s230_anchor_identity_audit.csv`.  It protects the historical-common
anchors; a newly admitted product that becomes an anchor must also pass the
separate exact addition/component-identity review.

A filename inventory of the working Documents trees, the live research tree,
and the pre-FEL hardware archive found no copy of `get_chosenfreq.ipynb`.
The preserved legacy fitter retains only its execution call, immediately after
the per-source minimum-geometric-beam model cleanup.  The replacement therefore
implements that visible rule and the missing model-copy handoff directly,
rather than trying to reconstruct unreviewable notebook state.

The deterministic generic rule is: detected candidates from 200--400 GHz with
geometric-mean restoring beam below 1 arcsec, then the smallest geometric-mean
beam.  Exact ties are resolved by portable input product, telescope, model
name, and run context.  A generic winner is not automatically accepted as a
nucleus.  The fitted coordinate must remain within the explicitly recorded
tolerance of the reviewed component, even when the product path matches.  A
product, frame, coordinate, or disposition mismatch stops consolidation.

For the eight retained generic anchors in the frozen historical-common set, the reviewed coordinate is the
preserved historical model peak after the independent identity disposition in
the S230 audit, and the guard radius is one historical restoring-beam major
axis.  That radius is a component-identity tolerance, not an astrometric
uncertainty.  NGC 4945 instead records the Emig et al. (2020) source-18
coordinate; its generic candidate is rejected regardless of separation.

## NGC 4945 is an explicit exception

`chosen_frequency_anchor_reference.csv` is an historical identity reference,
not a list of nine scientifically accepted nuclei.  Its NGC 4945 row records
the rejected source-17 SSC/SNR-associated seed at
`13:05:27.4933 -49:28:05.1422` because that is what the preserved 409-row
historical final-stage run consumed.

The current workflow never treats that seed as final nuclear photometry.
In schema 5, `prepare` authenticates the shipped
`source18_initial_history/` assets: 27 historical photometry records and 66
historical correction records, with public-release derivation provenance.
It installs them only in the new copied workspace and preserves immutable
`source18_bootstrap/` snapshots. The 409-row reference supplies product keys,
not historical OLD flux/SNR values. Preparation neither fits sources nor
replays the newly adopted results. It is independent of legitimate changes
to final diagnostic fields after a reviewed run, while preserving the exact
historical lineage and rejecting reused or already-final working metadata.
After the four generic passes, `run_continuum_four_pass.py` must run
`rerun_ngc4945_nuclear.py --apply` at the Emig et al. (2020) source-18
AGN/maser coordinate and then run its guard check.  Failure of either command
fails the workflow; the generic fixed-to-chosen table is only an intermediate
product.

## Search-region geometry

The search region is intentionally an isotropic square with both full side
widths derived from the restoring-beam major-axis FWHM (`Bmaj`).  It is not a
rotated beam ellipse, so `Bmin` is used for the geometric-mean resolution
eligibility test but not for the square search width.

- `prelimchecker` and `fit_sizenoloc`: `10 Bmaj` square.
- `getchoosercoords`: `max(10 Bmaj, 1 arcsec)` square, with the reviewed
  NGC 3393 one-arcsecond override.
- `fixtochosenfreq` with an anchor: `3/(2 sqrt(2 ln 2)) Bmaj` square.  This
  is a legacy full width of three Gaussian standard deviations (`1.274 Bmaj`),
  and therefore a half-width of `1.5 sigma`; it is not a `3 sigma` radius.
- `fixtochosenfreq` without an anchor: `min(20 Bmaj, 3 arcsec)` square; the
  identified NGC 3079 product explicitly uses the three-arcsecond width.
- A fit explicitly classified by `write_region` as extended or
  multiple-component retains the historical 10-arcsecond square override.

`square_bmaj_region_width_arcsec()` is the single tested implementation of
these full side widths.  The old unused minor-axis variable has been removed.

## Final publication-noise geometry

The final publication classification deliberately decouples its noise aperture
from the smaller fit aperture.  For each ordinary final image fit, the code
uses an unrotated square centered on the fitted or fixed component with full
width `10 Bmaj`, intersects
it with the valid image footprint, and excludes the fixed-beam `3 sigma`
ellipse.  It adopts `1.482602218505602 * medabsdevmed` exactly once as the
Gaussian-equivalent local scatter, requires at least 50 effective independent
beam samples, and uses fitted flux divided by this scaled MAD for the reported
S/N and the `S/N >= 3` classification.  An ordinary upper limit is the larger
of the fitted flux and three times the same scaled MAD.  Classic RMS, standard
deviation, mean, edge truncation, and RMS/scaled-MAD disagreement remain
audited diagnostics rather than row-dependent alternate estimators.

The 12 predesignated SMA `nitercal_0` dirty-image failures bypass this ordinary
classifier and retain their reviewed one-beam RMS upper-limit policy.  Their
Gaussian-fit flux is diagnostic only: a converged fit records the fitted
component flux, while a nonconverged no-component fit records that diagnostic
as blank without changing the independent three-RMS limit.  No widened
final-fit retry is used: the preserved retry block was unreachable in the
ordinary workflow and was removed rather than activated without a new
scientific review.

## Safe execution

The fitter creates and removes generated directories, and CASA can update lock
metadata in CASA table directories.  Do not run the four passes in the release
Git worktree or in the preserved hardware archive.

The deterministic `prepare` command builds this layout:

- a lightweight copy of the release without `.git`;
- real writable `ALMA/<source>/` and `VLA/<source>/` directories containing
  independent copies of the exact 529-file Task 28 whitelist in
  `frozen_fitter_inputs.csv`: 222 ALMA and the historical 307 VLA products,
  not the larger set returned by a new archive query;
- real copied SMA CASA image/residual directories, never directory symlinks;
  the 15 hash-locked SMA discovery images; and
- the six hash-locked selected-clean support trees (three clean image/residual
  pairs) together with their date/flag-filtered manifest.

The default copy mode gives OS-level separation from the preserved archive.
Before writing the workspace, `prepare` validates the frozen manifest's own
SHA-256 and every declared product hash, then uses only its exact paths.  Each
of the 89 addition-policy rows separately binds the product bytes, expected
two-decimal output frequency, UTC date label, and expected run disposition.
`prepare` reads
`DATE-OBS` and the FITS frequency-axis `CRVAL` directly and rejects any mismatch.
It also stops before CASA if an addition shares source, telescope, rounded
frequency, and UT date with a baseline or another addition, because suffix
stability must be implemented explicitly rather than silently renumbering a
historical common row.
An explicit `--archive-input-mode symlink` exists for low-disk diagnostics,
but is not recommended for the production rebuild because the process could
technically write through a file link.  In either mode, directory links are
rejected.

Set these paths for your environment; the audited run is documented in the
validation evidence.  `prepare` refuses an existing destination and does not
import or run CASA:

```bash
RELEASE=/path/to/wmagn_candidate
ARCHIVE=/path/to/augmented_archive/Burridge2025-WM-AGN-Analysis
WORKSPACE=/path/to/validation_workspace/release
CASA_PY=/path/to/new-wmagn-runtime/wm-agn-casa-python
ANALYSIS_UTILS=/path/to/new-wmagn-runtime/analysis_utils/analysis_scripts

python3 "$RELEASE/analysis/multi_freq_from_archive/run_continuum_four_pass.py" prepare \
  --release "$RELEASE" \
  --archive-project-root "$ARCHIVE" \
  --workspace "$WORKSPACE" \
  --date 230316 --flagnum 2 \
  --archive-input-mode copy \
  --sma-discovery-relative analysis/sma_calibration/230316_128/images/nitercal_0/flag2_cal1/fields/science

python3 "$WORKSPACE/analysis/multi_freq_from_archive/run_continuum_four_pass.py" plan \
  --workspace "$WORKSPACE" \
  --python-bin "$CASA_PY" \
  --analysis-utils-dir "$ANALYSIS_UTILS" \
  --date 230316 --flagnum 2 --calnum 1 \
  --smapath analysis/sma_calibration/230316_128/images/nitercal_0/flag2_cal1/fields/science \
  --selected-clean-manifest analysis/sma_calibration/selected_clean_iterations.csv \
  --use-selected-clean \
  --canonical-reference analysis/multi_freq_from_archive/fitsumfiles/historical/fitsummary_fixtochosen_sorted_409_reference.txt \
  --exclusion-policy analysis/multi_freq_from_archive/continuum_product_exclusion_policy.csv \
  --run-id augmented-local-scaled-mad-230316-flag2

python3 "$WORKSPACE/analysis/multi_freq_from_archive/run_continuum_four_pass.py" run \
  --workspace "$WORKSPACE" \
  --python-bin "$CASA_PY" \
  --analysis-utils-dir "$ANALYSIS_UTILS" \
  --date 230316 --flagnum 2 --calnum 1 \
  --smapath analysis/sma_calibration/230316_128/images/nitercal_0/flag2_cal1/fields/science \
  --selected-clean-manifest analysis/sma_calibration/selected_clean_iterations.csv \
  --use-selected-clean \
  --canonical-reference analysis/multi_freq_from_archive/fitsumfiles/historical/fitsummary_fixtochosen_sorted_409_reference.txt \
  --exclusion-policy analysis/multi_freq_from_archive/continuum_product_exclusion_policy.csv \
  --run-id augmented-local-scaled-mad-230316-flag2
```

`plan` is read-only.  It verifies all 544 discovery inputs (222 ALMA, 307 VLA,
and 15 SMA), the six selected-clean support trees, the exact
`frozen_fitter_inputs.csv` SHA-256, all 529 ALMA/VLA file hashes, the workspace
marker and workflow
input hashes, the complete frozen analysis-scripts bundle, required CASA and
Python imports, empty output targets, and the current comparison inputs before
printing the six exact commands (four generic passes, source-18 apply, and
source-18 check).  The 89-row Task 28 admission policy is scientific review
metadata within that already frozen 222-product ALMA membership; it is not an
extra set added to a smaller preparation baseline.

The runner refuses Git worktrees, missing/mismatched workspace markers,
fitters changed after workspace initialization, directory-symlinked SMA
inputs, and pre-existing generated outputs.  It gives every legacy stage
closed standard input, so a forgotten `input()` call fails instead of hanging,
and external notifications remain disabled.  The plan records SHA-256 digests
for every declared workflow input, including the fitter, consolidator,
coordinate/geometry/noise helpers, identity policies, the exact boundary-
rejection policy and evidence, source coordinate input,
historical reference, source-18 correction, selected-clean manifest, canonical
comparison table, and exclusion policy.  It also copies the complete required
analysis-scripts bundle into the evidence directory, executes that frozen copy,
and verifies both the bundle and workflow-input hashes immediately before and
after every stage.  Each stage's command, log, and output-tree hashes are preserved under
`workflow_runs/<run-id>/`.  The runner also checks that all prepared input-file
inventory evidence remains unchanged: it rehashes all 529 ALMA/VLA FITS files,
the 15 copied SMA discovery trees, and the six selected-clean support trees,
while rejecting directory symlinks and any missing, extra, or substituted
frozen identity.

The 529 preserved ALMA/VLA FITS inputs total 6,235,773,120 bytes (about
5.81 GiB); preparation therefore needs at least 6 GiB for those copies alone.
A previous isolated final pass used about 8.2 GB of
generated CASA model trees and took 67 minutes.  Because four mode-specific
fit products are retained, budget 40 GB and preferably 50 GB free, and allow
about 4--6 hours for the full run plus the source-18 correction and hashing.

## Long runs over SSH

Do not rely on an interactive terminal connection to preserve a multi-hour
fit. On a Linux host with an active systemd user manager, the same reviewed
`run` command can be launched as a separate user service after preparation
and a successful `plan`:

```bash
systemd-run --user --unit=wm-agn-continuum --service-type=exec \
  --property="WorkingDirectory=$WORKSPACE" \
  /usr/bin/python3 "$WORKSPACE/analysis/multi_freq_from_archive/run_continuum_four_pass.py" run \
  --workspace "$WORKSPACE" \
  --python-bin "$CASA_PY" \
  --analysis-utils-dir "$ANALYSIS_UTILS" \
  --date 230316 --flagnum 2 --calnum 1 \
  --smapath analysis/sma_calibration/230316_128/images/nitercal_0/flag2_cal1/fields/science \
  --selected-clean-manifest analysis/sma_calibration/selected_clean_iterations.csv \
  --use-selected-clean \
  --canonical-reference analysis/multi_freq_from_archive/fitsumfiles/historical/fitsummary_fixtochosen_sorted_409_reference.txt \
  --exclusion-policy analysis/multi_freq_from_archive/continuum_product_exclusion_policy.csv \
  --run-id augmented-local-scaled-mad-230316-flag2

systemctl --user status wm-agn-continuum.service
journalctl --user -u wm-agn-continuum.service --no-pager
```

Use new unit/run names and a new prepared workspace for a new attempt; do not
reset partial outputs or bypass the no-reuse guards. This separates the job
from the client terminal, not from the workstation: the host and user manager
must remain available. It changes no fitting arguments or scientific settings.
An active service is not a successful completion. After it exits, inspect all
six stage snapshots and final validation records. The intentional
`review_required` stop also has a nonzero exit; distinguish it from a failed
stage before accepting or promoting anything. If a connection is lost, check
the service and actual processes before deciding the computation stopped.

## Validation evidence and limit

The checked-in historical reference records nine anchor product identities,
input SHA-256 digests, model-tree SHA-256 digests, WCS frames, beam values, and
model-peak coordinates.  All nine preserved generic identities are ALMA
products.  The completed isolated historical final-stage validation separately
reproduced all 409 accepted rows and all 46 exclusions exactly.  Its counts,
CASA version, output digests, code provenance, source-evidence digests, and
environment limitation are preserved compactly in
`historical_final_stage_validation.json`.

That 409/46 result validates the historical fourth stage and the preserved
historical model trees.  It does not substitute for executing the corrected
four-pass workflow: current WCS routing, component-policy checks, NGC 4258
product exclusions, and the mandatory NGC 4945 source-18 correction are
intentional changes.  A new full run must be reviewed as a corrected rebuild,
not required to be byte-identical to the rejected historical component seed.

The immutable historical raw baseline has 409 rows because it retains three
NGC 4258 NVAS products now explicitly rejected by
`continuum_product_exclusion_policy.csv`.  Two are failed reductions replaced
by the same-observation Doi et al. (2013) measurements; the AH0847 2004-05-22
image is independently rejected because bright diagonal stripe artifacts cross
the nucleus, and it has no replacement row.  The corrected fitter loads this
hash-locked policy, verifies each rejected FITS artifact before skipping it,
and fails closed if the policy is missing or malformed.  The immutable
common-product expectation therefore contains 406 rows.  Product admission is
separate: all `N` recovered inputs are byte/header checked before CASA, while
only the policy subset `P` with `publication_row_required` may appear in the
fresh raw table.  The exact generated count is `406 + P`; with all 89 Task 28
admissions marked `publication_row_required`, that count is 495.  An input-only
admission must remain absent.  After the mandatory NGC 4945 correction, the
runner compares product identity and all seven reported fields against the
frozen 409-row pre-scaled-MAD baseline in
`fitsumfiles/historical/fitsummary_fixtochosen_sorted_409_reference.txt`
after removing exactly those three reviewed products. It writes
`current_corrected_table_comparison.json`.  Product identity, order, and every
non-flux/SNR field must match exactly; flux/SNR differences caused by the final
scaled-MAD policy require an exact, product-complete reviewed-delta ledger and a
hash-frozen adjudication before the rebuild can be accepted.  Every admitted
input also has a distinct post-run row in
`reviewed_new_product_additions.csv`, including explicit no-publication
outcomes.  Detected additions require component-identity review; added upper
limits retain admission provenance without claiming an adopted component.  The
CASA-free, no-clobber `render-review-evidence` command derives both ledgers and
the exact freeze JSON before `adjudicate` recomputes them.  The later machine-table builder
adds the two reviewed Doi et al. replacements; those literature rows are not
raw fitter outputs and are deliberately outside this comparison.

The cumulative `validated_continuum_corrections.csv` contains 84 records from
the formatter, frame-routing, RMS-era noise, and NGC 4945 source-18 correction
lineage. Its run-specific preapplication snapshot contains the 66 historical
rows, and Task 28 adds the 18 newly admitted NGC 4945 products. It is correction
lineage rather than a complete account
of the source-excluded scaled-MAD changes.  The product-complete
run-specific `reviewed_current_to_generated_deltas.csv` and
`reviewed_new_product_additions.csv` ledgers, bound by
`reviewed_local_scaled_mad_freeze.json` and certified by
`reviewed_delta_adjudication.json`, are the authoritative final change scope.

For a release, copy those four files byte-for-byte from the adjudicated
`workflow_runs/<run-id>/` directory to exactly
`analysis/multi_freq_from_archive/reviewed_current_to_generated_deltas.csv`,
`analysis/multi_freq_from_archive/reviewed_new_product_additions.csv`,
`analysis/multi_freq_from_archive/reviewed_local_scaled_mad_freeze.json`, and
`analysis/multi_freq_from_archive/reviewed_delta_adjudication.json`.  Do not
rewrite paths or hashes inside them.  These root-level copies are the stable
canonical release evidence; the bulky run directory is execution evidence and
is not part of the software release.  After promoting the sorted table, noise
audit, identity policy, correction manifest, and NGC 4945 sidecars, run:

```bash
python3 analysis/multi_freq_from_archive/reviewed_release_evidence.py --check
```

The checker fails closed on an absent or partial release bundle, duplicate JSON
fields or exact JSON type/schema drift, any noncanonical path in the complete
22-artifact freeze inventory, a ledger/freeze/adjudication hash mismatch, a
released-runner mismatch, wrong final counts, a frozen release-side artifact
mismatch, or any difference from a fresh strict ledger recomputation.  The
test suite may skip its canonical integration check only while all four
root-level files are genuinely absent; the command above never skips and is a
mandatory release gate.

## Historical authoritative Task 28 v1.0.10 execution

This section records the original accepted scientific baseline. A repaired
candidate must generate its own current-run evidence with its current code;
the original runner hashes and old adjudication must not be rewritten to
make new code pass. See the candidate test report for the new execution.

Task 28 audited the filename-dependent ALMA selection rules and admitted 82
additional products beyond the seven admitted at the A8 checkpoint: 49
pre-cutoff filename-only omissions, 26 post-cutoff target continuum datasets,
and seven Stokes-I MFS products.  Relative to the immutable 409-row historical
baseline, the complete addition policy therefore contains 89 products.

The isolated run `task28-authoritative-v2-230316-flag2-20260825` prepares 544
inputs (222 ALMA, 307 VLA, and 15 SMA), retains 406 historical-common products,
and adds all 89 reviewed ALMA products.  Its raw publication output contains
495 rows and its complete diagnostic-noise union contains 498 rows, including
the three exact policy-bound boundary rejections.  The NGC 4945 fixed-nucleus
stage covers 45 products while preserving the exact pre-Task-28 lineage of the
27 previously released products.

After downstream bandwidth-smearing and alias exclusions, 485 local products
remain (222 ALMA, 248 NVAS, and 15 SMA).  Sixteen literature measurements bring
the release table to 501 rows.  The 498 rows belonging to the 21-source sample
produce 14,193 unordered within-source pairs and a 105-pair strict screen,
including 103 different-day and two same-day comparisons.  The selected
$S_{230}$ table contains 12 detections and seven upper limits drawn from 11
ALMA, seven SMA, and one NOEMA measurement.

The original v1.0.10 release promoted only the hash-frozen Task 28 reviewed-delta, addition,
freeze, and adjudication artifacts.  The A8 execution below remains a fully
validated historical checkpoint and must not be substituted for the Task 28
release evidence.

## Historical augmented execution checkpoint (A8)

The 462-input and seven-admission counts in this subsection are retained only
as the superseded A8 checkpoint; the current Task 28 membership is the frozen
544-input inventory described above.

The inclusive public-release census through 2025-01-27 found 112 ALMA Member
OUS records absent from the historical inventory.  Seven are eligible
additions and 105 have documented exclusions.  Six additions are ESO products
that a legacy raw `SO` substring test mistakenly rejected by matching the
letters within `ESO`; the seventh is an independently omitted NGC 1068 product.
The shipped `get_urls.ipynb` now matches chemical tokens against the product
basename with token boundaries.
The lightweight 112-row census, seven-row admission manifest,
query/provenance record, and preserved pre-fix notebook are packaged under
`archive_inventory_audit/`; that directory intentionally contains no FITS.

The final external inventory is consequently 462 inputs: 140 ALMA FITS files,
307 VLA FITS files, and 15 SMA images.  The isolated schema-4 run
`cgcg-augmented7-local-scaled-mad-230316-flag2-authoritative8` completed all four
CASA passes, the NGC 4945 source-18 apply/check, final inventory and identity
validation, no-clobber review rendering, and exact adjudication on 2026-08-20
(HST).  Exact review accepted a 413-row publication table (`406 + 7`) and a
416-row generic diagnostic audit: the same 413 publication products plus the
exact three policy-bound converged fits rejected by the historical
half-maximum-contour boundary gate.  The immutable
`fitsumfiles/historical/fitsummary_fixtochosen_sorted_409_reference.txt` remains
unchanged and distinct from both the historical A8 output and current Task 28 output.

The A8 plan has SHA-256
`00b0b3dff3c5cb00920934eaef87f503a4c99172e0fd7e2074d1eb0eff40c9eb`.
The promoted 413-row sorted table and 416-row noise audit have SHA-256 values
`2756c1b5abc1d548dea76fba826602a8ddbe04b18424aded4cdf7c68bdc839c9`
and `f4c4a1e2e3170e6b3a2db222305bfdd86f70a2166ddd411c5141849b57f6c93c`.
The exact reviewed delta ledger, seven-row addition ledger, schema-3 freeze,
and schema-2 adjudication have SHA-256 values
`0da4af0439d549c4a62759d20d2f80291e5eafd15c20288801406a97c810a2da`,
`479a34f7bfb9d61d2a835d96542f8409d44011158f9aa2ed32a4d96432b668fe`,
`5aa6b7c2fdefaf20b9fa751424a58ac1c5d8fe807b31c7f75fe3c56f0670322d`,
and `f29c75068f4c1329156586f4bf6a68daa6e2bc331709a79e8e11464c62eb6234`,
respectively.  The adjudication reports 393 reviewed delta rows, seven reviewed
additions, three boundary rejections, 19 classification identities, 27 NGC 4945
photometry rows, and seven NGC 4945 component records.

Downstream bandwidth-smearing, alias, and literature assembly produced 416
rows: 403 retained local measurements (140 ALMA, 248 NVAS, and 15 SMA) plus 13
literature-only rows.  The recovery manifest totals 3,916,034,030 local payload
bytes.  Its 413 analysis-sample measurements produce 10,137 unordered
within-source pairs and a 95-pair strict screen, including two same-day pairs.
The 19-row selected-anchor
audit contains 12 detections and seven upper limits drawn from ten ALMA and nine
SMA observations.  The admitted ESO 558-G009 X600 row is a 313.11 GHz,
0.66 mJy detection in a 0.71817-arcsec beam and gives
`S_230 = 0.8190724285773601 +/- 0.09125138464542781 mJy`.

## Historical pre-scaled-MAD corrected execution

The 455-input counts below likewise describe the earlier pre-Task-28 run, not
the current fitter membership.

The RMS-based corrected workflow completed successfully on 2026-08-18 (HST) in
the isolated copied workspace under run ID
`current-corrected-230316-flag2-retry3`.  The evidence directory is retained in
the companion correction audit at
`end_to_end_validation/current_four_pass_20260818_retry3/release/analysis/multi_freq_from_archive/workflow_runs/current-corrected-230316-flag2-retry3/`.

Independent read-only checks verified all of the following:

- all four generic stages completed and each has a stage-artifact snapshot;
- all nine chosen-frequency anchors passed the product, coordinate, frame,
  component-policy, and model-provenance check;
- the mandatory NGC 4945 source-18 apply and check stages passed for all 27
  measurements and all seven component-audit records;
- the prepared 455-input inventory remained unchanged: 133 ALMA and 307 VLA
  FITS files passed the content-hash recheck, while the 15 copied SMA discovery
  images passed the runner's count and no-directory-symlink guards; an
  independent post-run deterministic tree-hash comparison additionally
  matched all 15/15 SMA discovery trees and all 6/6 selected-clean
  image/residual trees to the original archive;
- the final table contains 406 unique rows, all product paths are
  release-relative and resolve inside the isolated release, and all 15 SMA
  paths use `analysis/sma_calibration/...`; and
- the generated table is exactly equal, in order and in all seven reported
  fields plus product identity, to the frozen 409-row historical reference after
  removing only the three reviewed NGC 4258 products.

The generated 406-row table has SHA-256
`da0c6aa568b70cbba0899926538a9e2ae16de26d44ec4b76823f628d5f5f6778`.
The exact-comparison and input-integrity records have SHA-256 values
`cc04f7c735f7ea1203d621a5e25a04a6de9d57bee554a236044cf24472b50e7a`
and `ea6928f95fe943cb9ca928b69c32a9ccf12e02c04f0a660da98da4911518edd4`,
respectively.  The aggregate hashes of the sorted SMA tree-name/hash mappings
are `241b550ab23cbd9f5560248c0ae731fc37e81b339fbcbd5ee8b6aa7364c120c1`
for the 15 discovery trees and
`4cab219804e4d2b6f22f25a6dd5cbfc6d683e43fdd0f887ae3631bdb8ebfc4a8`
for the six selected-clean trees. Those historical SMA hashes were
independent post-run evidence; that earlier runner did not yet record them in
its prepare/integrity JSON. The current runner records and checks both sets.
This completed the earlier RMS-based four-pass execution gate and is retained
as audit history; it does not reproduce the final source-excluded scaled-MAD
publication classifier.  The separate IC 2560
author-decision gate was subsequently closed on 2026-08-19 with an explicit
formal-only treatment.  The v1.0.9 DOI reservation and insertion gate closed
on 2026-08-21 HST with Zenodo version DOI
`10.5281/zenodo.22051996`; that version was subsequently published.
Version v1.0.10 was subsequently published at `10.5281/zenodo.22121490`.
Publication of this repaired candidate remains an external authorization gate.
