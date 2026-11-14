# WM-AGN Analysis Code

## Software-repair release: v1.0.11

This release accompanies the published **ApJ article** cited below
(manuscript AAS72796R3 / AE8A42).
The published Zenodo record for v1.0.11 is
[`10.5281/zenodo.23759601`](https://doi.org/10.5281/zenodo.23759601).
Use this version-specific identifier when citing v1.0.11.
The original v1.0.10 archive remains unchanged at
`10.5281/zenodo.22121490`; that DOI identifies the older version.

v1.0.11 promotes the validated v1.0.11-rc1 package with version metadata and
documentation changes only. Scientific source, tests, dependencies, numerical
products and figures are byte-identical to that candidate. The original
`CANDIDATE_TEST_REPORT.txt` and `RELEASE_NOTES_v1.0.11-rc1.md` are retained as
historical validation evidence. Their 467 passed / seven skipped test results
were obtained for rc1; the scientific tests were not rerun for this metadata-only
promotion. See `RELEASE_NOTES_v1.0.11.md` and `STABLE_RELEASE_PROMOTION.txt`
for the current release description and independent packaging checks.

Start with `INSTALL.md` for the tested clean Python/CASA installation and
`DATA_ACCESS.md` for the external observations. The frozen fitter manifest
requires exactly 544 discovery products plus six selected-clean supports.
Archive discovery is not allowed to silently add later query results.
Preparation creates an independent workspace and initializes the authenticated
NGC 4945 historical metadata there; the actual fitting and correction must
still run. Published measurements are not replayed as a substitute for fits.

This is a code-and-lightweight-products archive, not a bundle of raw telescope
data or a preinstalled CASA environment. Exact historical fitting requires the
hash-matching preserved input products. Fresh SMA archive processing remains
a separate reproducibility limitation: an independently delivered/reprocessed
NGC 4258 track differed from the adopted result and has not replaced it.
See the original `CANDIDATE_TEST_REPORT.txt` for the completed rc1 checks
and remaining limits inherited by this release. Installation, data retrieval, fitting, review,
and rendering are explicit stages; this is not a claim of unattended recovery
of every historical SMA product from today's archive.

This archive contains the analysis code and lightweight machine-readable tables for the manuscript

**When the Shadow Meets Its Measure: Assessing the Feasibility of Submillimeter Black Hole Shadow Imaging in Megamaser-disk Active Galactic Nuclei**

Roman N. Burridge and Geoffrey C. Bower (2026), *The Astrophysical Journal*,
**1011**, 273. Journal DOI:
[`10.3847/1538-4357/ae8a42`](https://doi.org/10.3847/1538-4357/ae8a42).

The code calibrates and images SMA observations, gathers and fits archival ALMA/VLA continuum products, estimates beam-matched continuum screening flux densities, converts the selected 200--400 GHz reference measurements to a common 230 GHz frequency, computes black-hole-shadow angular sizes, estimates contamination from thermal dust and extended jets, and produces the tables/figures used in the paper.

Large telescope data products are intentionally not included in this release. See `DATA_ACCESS.md` for instructions on retrieving ALMA, VLA/NVAS, and SMA data from their public archives.

This release is distributed under the MIT License; see `LICENSE`.

The original version v1.0.10 superseded the published v1.0.9 archive for reproduction of
the final revised manuscript.  It incorporates the authoritative Task 28
continuum closeout, expands the machine-readable continuum table from 425 to
501 rows, and regenerates the dependent 230 GHz, contamination, table, and
figure products.  The black-hole parameter and uncertainty corrections
introduced in v1.0.9 are retained and tested here.  The v1.0.10 Zenodo version
DOI is `10.5281/zenodo.22121490`.  IC 2560 is explicitly retained as a
formal-only result: the shared
distance cancels, no separately quantified model systematic is available, and the
unquantified model contribution is not treated as zero.  Its physical-mass
bounds likewise contain only the available fitting and distance components,
not a complete uncertainty total.  A numerical IC 2560 fit that appeared only
in arXiv:2409.00772v1 is audit context rather than a release comparison: v2
removed that analysis because the input data were not officially published,
and the peer-reviewed 2025 article contains no IC 2560 fit.  The exact scope and
verification record are in `CORRECTION_CHANGELOG.md`.  For NGC 4258, the
adopted paired-chain shadow-size result and its separate 1 percent
non-Keplerian allowance are retained.  A later component-selection check using
a subset of the public VLBI epochs is recorded separately as a non-probabilistic sensitivity and is not
folded into the adopted calculation or figure error bar.

The historical v1.0.10 continuum-release architecture was completed by the isolated
schema-4 Task 28 run
`task28-authoritative-v2-230316-flag2-20260825` and its exact review
adjudication. The repaired runner uses schema 5, frozen input membership,
and separately authenticated initial-history assets. Its current-run evidence
must be generated by actual execution, not copied from the historical run.
The prepared external inventory contains 544 inputs (222 ALMA,
307 VLA, and 15 SMA).  Relative to the immutable 409-row historical baseline,
the reviewed publication fit table contains 406 retained historical-common
rows and 89 admitted ALMA products, for 495 rows.  Its complete diagnostic
noise audit contains those 495 rows plus three exact, hash-bound, converged
fits rejected by the historical half-maximum-contour boundary gate.  After
the nine-row VLA bandwidth-smearing policy and the one-row IC 2560 alias policy
are applied and 16 literature measurements are added, the final machine table
and recovery manifest contain 501 measurements: 485 local products (222 ALMA,
248 NVAS, and 15 SMA) plus 16 literature-only records.  Of these, 498 belong
to the 21-source analysis sample and three are the additional NGC 5495 catalog
measurements.  The regenerated 19-source $S_{230}$ audit contains 12 detections
and seven upper limits selected from 11 ALMA, seven SMA, and one NOEMA
measurement.  Counts below that describe older states are labeled historical
and must not be used as final-release checks.

The inclusive archive census, admission manifests, query/provenance records,
and preserved pre-fix discovery materials are packaged under
`analysis/multi_freq_from_archive/archive_inventory_audit/`; that directory
contains no archival FITS bytes.  The canonical review ledgers, freeze,
adjudication, product hashes, and recovery manifest bind the final release to
the authoritative Task 28 run.

## Contents

- `BH_Parameters.ipynb`: readable diagnostic companion for the audited black-hole parameter workflow; it is not a hidden table-building engine.
- `VLBI_sensitivity.ipynb`: computes representative VLBI sensitivity and angular-resolution quantities.
- `analysis/CountingBHS_and_Spectra/bhs_source_inputs.csv`: canonical 21-source mass, distance, treatment, uncertainty-component, and citation inputs, including a separate DOI (or `NA`) for each adopted mass uncertainty.  The UGC 6093 record retains Zhao et al.'s asymmetric fitted-velocity interval and separately typed, source-reported absolute-reference-position systematic.  The ESO 558-G009 record retains the exact `27/73.9 Mpc` velocity-derived measurement term and a separately typed midline/inclination geometry systematic.  The no-primary Mrk 1029 record uses the selected Huchra velocity `9076 +/- 32 km/s`, hence the exact `32/73.9 Mpc` measurement component, while preserving the second published velocity as an ambiguity in provenance.  The no-primary NGC 1320 record likewise preserves Gao et al.'s selected `2649 +/- 16 km/s` value and exact `16/73.9 Mpc` measurement component.
- `analysis/CountingBHS_and_Spectra/bhs_latex_metadata.csv`: source typography, inclination text, and BibTeX keys used to render manuscript Table 1 without duplicating scientific numbers.  NGC 1320 is reported as `unconstr.`; the `i=90 degrees` used by later restricted fits is a fixed assumption, not a measured inclination.
- `analysis/CountingBHS_and_Spectra/pesce_theta_g_quantiles.csv`: sanitized six-source paired-chain angular-gravitational-radius quantiles supplied for this correction.
- `analysis/CountingBHS_and_Spectra/bhs_model_sensitivity_checks.csv`: ten attached non-probabilistic checks for NGC 4258, NGC 2960, NGC 6323, J0437+2456, Mrk 1029, NGC 1320, Circinus, NGC 3079, NGC 4388, and NGC 4945.  The NGC 4258 branch records the VLBI data-selection configurations and endpoint masses and distances from van der Boom & Kushnir (2026; DOI `10.1093/mnras/stag1327`); NGC 2960 and NGC 6323 record the three physically allowed Kuo et al. (2018; DOI `10.3847/1538-4357/aabff1`) disk configurations and exact one-sided endpoint ranges; J0437+2456, Mrk 1029, and NGC 1320 record separate edge-on Mestel-disk mass spans from Deepshikha et al. (2024; DOI `10.1017/S1743921323002491`), with the latter two standalone diagnostics because they have no adopted primaries.  The remaining four branches record the one-sided 3, 34, 0.1, and 24 percent turbulent-support mass-bias estimates for Circinus, NGC 3079, NGC 4388, and NGC 4945 from Kuo et al. (2020, Appendix A; DOI `10.1093/mnras/staa2260`).
- `analysis/CountingBHS_and_Spectra/bhs_audit_only_model_sensitivity_checks.csv`: one separately validated, audit-only, unattached NGC 1320 record preserving Kuo et al.'s (2020, Appendix A) conditional analytic 1 percent turbulent-support result.  With no robust primary mass, it cannot attach to a generated row or become an adopted uncertainty, shadow size, or figure point.
- `analysis/CountingBHS_and_Spectra/bhs_nonadopted_ratio_fits.csv`: thirteen conditional literature $M/D$ fits from Villaraos et al. (2022) and González-Juárez et al. (2024), including retained-primary comparisons for NGC 1194, NGC 2273, NGC 2960, NGC 4388, NGC 5765b, CGCG 074-064, NGC 6264, UGC 3789, NGC 6323, ESO 558-G009, and J0437+2456 and the non-adopted Mrk 1029 and NGC 1320 fits, recorded separately from the adopted results with their model assumptions and non-adoption bases.
- `analysis/CountingBHS_and_Spectra/build_bhs_params.py`: deterministic black-hole parameter and uncertainty builder; it validates all three non-adopted-check sidecars and attaches eligible fields only after constructing each generated result, while validating but never attaching the separate audit-only record.
- `analysis/CountingBHS_and_Spectra/BHS_params.csv`: generated machine-readable black-hole-shadow parameter table, including separate non-adopted model-check and literature-fit fields.
- `analysis/CountingBHS_and_Spectra/build_bhs_latex_table.py`: deterministic manuscript Table 1 renderer driven by `BHS_params.csv` and `bhs_latex_metadata.csv`.
- `analysis/multi_freq_from_archive/build_bhs_flux_figure.py`: deterministic manuscript BHS-size versus 230 GHz flux figure builder.
- `analysis/multi_freq_from_archive/build_sed_figures.py`: deterministic six-panel SED builder; its angular-to-linear insets are derived from the canonical adopted-distance table.
- `analysis/multi_freq_from_archive/build_continuum_latex_excerpt.py`: deterministic renderer for the ten-row printed excerpt of the full continuum machine table.
- `analysis/multi_freq_from_archive/SpectraGraphs/sed_figure_manifest.json`: row counts, distance-derived inset scales, and hashes for the six SED PNGs.
- `analysis/sma_calibration/`: SMA calibration and imaging scripts, observing-script text files, molecular-line flag list, and the SMA MIR-to-measurement-set conversion helper.
- `analysis/multi_freq_from_archive/`: continuum fitting, archival-flux compilation, 230 GHz extrapolation, contamination estimates, variability-screening notebooks, and plotting notebooks.
- `analysis/multi_freq_from_archive/continuum_product_alias_policy.csv`: audited retained/excluded archive-product policy, including hashes and matching FITS-header evidence for the two IC 2560 NVAS reductions of one observation.
- `analysis/multi_freq_from_archive/continuum_product_exclusion_policy.csv`: audited failed-product policy for three NGC 4258 NVAS reductions, including immutable hashes, headers, fit and image-artifact evidence, and reasons for exclusion. Two observations have exact Doi et al. (2013) literature replacements; the artifact-dominated 2004-05-22 AH0847 image has no replacement row.
- `analysis/multi_freq_from_archive/continuum_bandwidth_smearing_exclusion_policy.csv`: audited nine-product VLA exclusion policy implementing `x <= 0.5`; `machinetables/vla_bandwidth_smearing_exclusions.csv` preserves every numerical input and generated disposition.
- `analysis/multi_freq_from_archive/continuum_boundary_rejection_policy.csv` and `continuum_boundary_rejection_evidence.csv`: exact three-product contract for converged fits excluded by the historical half-maximum-contour boundary gate.  The products were retained in the historical A8 416-row diagnostic audit and remain in the final Task 28 498-row diagnostic-noise union while staying outside the corresponding 413- and 495-row publication tables; input bytes, independent rerun evidence, boundary metrics, and reproduced CASA model payloads are hash-bound.
- `analysis/multi_freq_from_archive/build_continuum_recovery_manifest.py`: deterministic builder and verifier for the final 501-row continuum recovery inventory.
- `analysis/multi_freq_from_archive/continuum_product_checksums.csv`: the final audit-time size and SHA-256 lock covers all 485 locally fitted products; CASA image directories use a relative-tree checksum that excludes workstation paths and filesystem metadata.
- `analysis/multi_freq_from_archive/machinetables/continuum_recovery_manifest.csv` and `continuum_recovery_manifest_summary.json`: one-to-one archive/literature provenance for every canonical continuum row, including stable record IDs, ALMA Member OUS UIDs, NVAS programs and filenames, SMA program/date/image paths, retrieval URLs, and explicit literature-only status.
- `analysis/multi_freq_from_archive/continuum_source_roles.csv`: machine-verifiable inventory of the 21-source analysis sample and the additional NGC 5495 catalog entry.
- `analysis/multi_freq_from_archive/validated_continuum_corrections.csv`: cumulative 84-row formatter/frame/earlier-noise and NGC 4945 source-18 lineage.  Its run-specific preapplication snapshot contains the 66 historical rows; Task 28 adds the 18 newly admitted NGC 4945 products.  It is not the complete final scaled-MAD change inventory; the canonical `analysis/multi_freq_from_archive/reviewed_current_to_generated_deltas.csv` and `analysis/multi_freq_from_archive/reviewed_new_product_additions.csv` ledgers own that scope.
- `analysis/multi_freq_from_archive/reviewed_release_evidence.py`: fail-closed validator for the four byte-preserved canonical review artifacts: the two ledgers above, `reviewed_local_scaled_mad_freeze.json`, and `reviewed_delta_adjudication.json`.  The checker binds their exact schemas and hashes to one authoritative run and to the promoted release-side table, policies, noise audit, identity policy, and NGC 4945 sidecars.
- `analysis/multi_freq_from_archive/continuum_product_addition_policy.csv`: 89 hash- and header-locked ALMA admissions derived from the inclusive archive census; all 89 require publication rows.
- `analysis/multi_freq_from_archive/archive_inventory_audit/`: lightweight census and admission evidence, query/provenance materials, and the preserved pre-fix discovery notebook; no FITS products.
- `analysis/multi_freq_from_archive/fitsumfiles/historical/fitsummary_fixtochosen_sorted_409_reference.txt`: immutable 409-row pre-scaled-MAD baseline, kept distinct from the live 495-row sorted fitter output.
- `analysis/multi_freq_from_archive/ngc2960_frame_correction_audit.json`: immutable-input and coordinate-frame evidence for the historical pre-scaled-MAD `0.47*` mJy decision.  Its old RMS classification is superseded by the exact final scaled-MAD ledger and must not be cited as the current classification.
- `analysis/multi_freq_from_archive/rerun_ngc4945_nuclear.py`: deterministic fixed-centroid remeasurement of all 45 NGC 4945 continuum images at the Emig et al. (2020) source-18 nucleus; the 27 previously released rows retain their exact pre-Task-28 lineage.
- `analysis/multi_freq_from_archive/ngc4945_nuclear_photometry.csv` and `ngc4945_component_audit.csv`: per-image corrected photometry and the accepted/rejected component audit, including the archived source-17 result.
- `analysis/multi_freq_from_archive/s230_anchor_identity_audit.csv`: exact 19-source selected-observation inventory with nuclear-identity dispositions and caveats.
- `analysis/multi_freq_from_archive/build_s230_astrometry_audit.py`, `s230_astrometry_policy.csv`, and `s230_wcs_evidence.csv`: fail-closed builder inputs for the 19-row selected-anchor component-coordinate and hash-verified image-WCS audit.
- `analysis/multi_freq_from_archive/machinetables/s230_anchor_astrometry_audit.csv` and `s230_anchor_astrometry_schema.md`: machine-readable fitted/reference coordinates, ICRS transforms, signed offsets, beam separations, uncertainty/calibrator/QA availability, and accepted/ambiguous/blended/upper-limit dispositions.  The final 19-row audit is required to contain 12 detections and seven upper limits.  The NGC 5765B reference is explicitly classified as ICRS from its VLBA Calibrator Survey provenance (Beasley et al. 2002; DOI `10.1086/339806`), rather than transformed from FK5.
- `analysis/multi_freq_from_archive/chosen_frequency_component_policy.csv`: reviewed nine-source historical chosen-frequency product, coordinate, tolerance, and generic-seed actions derived from that S230 identity audit; admitted additions have their own reviewed identity ledger.
- `analysis/multi_freq_from_archive/chosen_frequency_anchors.py`: deterministic replacement for the missing chosen-frequency consolidation notebook, with reviewed component-policy enforcement and hashed CASA-model provenance.
- `analysis/multi_freq_from_archive/chosen_frequency_anchor_reference.csv`: nine-row historical anchor-identity reference; its NGC 4945 source-17 seed is explicitly rejected and retained only to validate historical provenance.
- `analysis/multi_freq_from_archive/historical_final_stage_validation.json`: compact, hash-locked record of the isolated historical 409-row/46-exclusion final-stage validation and its limits.
- `analysis/multi_freq_from_archive/run_continuum_four_pass.py` and `CONTINUUM_FOUR_PASS_WORKFLOW.md`: isolation-gated four-mode orchestrator and its exact geometry, safety, and validation protocol.
- `analysis/multi_freq_from_archive/machinetables/fitsummary_machine.txt`: full machine-readable continuum-fitting table corresponding to the printed table in the paper.
- `analysis/multi_freq_from_archive/machinetables/continuum_product_aliases.csv`: generated provenance sidecar preserving both IC 2560 product identities and the reason the older reduction is excluded from measurement calculations.
- `analysis/multi_freq_from_archive/machinetables/continuum_product_exclusions.csv`: generated audit sidecar preserving all three failed NGC 4258 products and their former table rows without treating any image as nuclear photometry.
- `analysis/multi_freq_from_archive/machinetables/S230_extrapolation_comparison.csv`: 230 GHz reference-flux table, including the adopted spectral index and the alpha range used to assess model sensitivity.
- `analysis/multi_freq_from_archive/build_variability_tables.py`: deterministic variability diagnostic-table generator.
- `analysis/multi_freq_from_archive/build_contamination_tables.py`: deterministic 230 GHz contamination-table generator.
- `analysis/multi_freq_from_archive/contamination_literature_inputs.csv`: literature measurements used by the contamination calculation, with explicit provenance and status.
- `analysis/multi_freq_from_archive/contamination_adoption_policy.csv`: explicit adopted, rejected, unavailable, or pending status for each contamination result.

The archival FITS files, raw SMA data, CASA image products, CASA logs, measurement sets, and local environment files are excluded.

## Python/CASA Dependencies

The analysis uses a CASA-capable fitting environment and a separate locked
publication-rendering environment. Follow `INSTALL.md`: it supplies tested
version pins, the compatible Linux/CPython 3.12 CASA wheel, the pinned
analysisUtils download, and isolated CASA data/configuration. The plotting
stack in the fitting environment is not the byte-reproduction renderer;
see `RENDERING_ENVIRONMENT.txt` and `scripts/render_publication.py`.

Core Python packages:

- `numpy`
- `scipy`
- `pandas`
- `matplotlib`
- `astropy`
- `astroquery`
- `requests`
- `pyuvdata` for SMA MIR conversion
- `jupyter`
- `reproject` (required by the actual fitter import)

CASA-side packages:

- CASA 6 with `casatasks` and `casatools`
- `analysisUtils` for coordinate conversion utilities used during ALMA fitting

The shell scripts accept `PYTHON_BIN=/path/to/casa-python` when a CASA-specific interpreter is needed.

## Reproducing the Main Analysis

Run notebooks from the release root so that relative paths resolve correctly.

1. Retrieve the external telescope data following `DATA_ACCESS.md`.
   The checked-in recovery inventory can be validated without the large data:
   ```bash
   python3 analysis/multi_freq_from_archive/build_continuum_recovery_manifest.py --check
   ```
   After arranging the augmented external archive tree, recompute and compare
   all 485 fitted-product hashes with:
   ```bash
   python3 analysis/multi_freq_from_archive/build_continuum_recovery_manifest.py \
     --check --archive-root /path/to/augmented/Burridge2025-WM-AGN-Analysis
   ```
   Validate the lightweight 19-row S230 astrometry join with:
   ```bash
   python3 analysis/multi_freq_from_archive/build_s230_astrometry_audit.py --check
   ```
   Re-extract and compare the actual selected-image WCS with CASA Python when
   the preserved project tree is available:
   ```bash
   PYTHON_BIN=/path/to/casa-python
   "$PYTHON_BIN" analysis/multi_freq_from_archive/extract_s230_wcs_evidence.py \
     --archive-root /path/to/augmented/Burridge2025-WM-AGN-Analysis --check
   ```
2. Place ALMA, VLA/NVAS, and SMA products under:
   - `analysis/multi_freq_from_archive/ALMA/`
   - `analysis/multi_freq_from_archive/VLA/`
   - `analysis/multi_freq_from_archive/SMA/`
   - `analysis/sma_calibration/sma_data/`
3. Convert SMA MIR data to CASA measurement sets if needed:
   ```bash
   cd analysis/sma_calibration/sma_data
   python3 sma2casa.py 230316_128
   ```
4. Calibrate/image SMA observations using a CASA-capable interpreter:
   ```bash
   cd analysis/sma_calibration
   /path/to/casa-python calibrate_sma_dirty.py --date 230316 --flagnum 2 --sciniter 0
   /path/to/casa-python calibrate_sma_clean.py --date 230316 --flagnum 2 --sciniter 6
   /path/to/casa-python calibrate_sma_clean.py --date 230316 --flagnum 2 --sciniter 7
   /path/to/casa-python calibrate_sma_clean.py --date 230316 --flagnum 2 --sciniter 10
   ```
   The three clean iterations are all required by
   `selected_clean_iterations.csv`: NGC 4258 uses 6, NGC 3079 uses 7, and
   NGC 1068 uses 10.
5. Fit continuum images:
   ```bash
   cd analysis/multi_freq_from_archive
   PYTHON_BIN=/path/to/casa-python bash run_continuum_fits_dirty.sh
   PYTHON_BIN=/path/to/casa-python bash run_continuum_fits_clean.sh
   ```
   The two shell runners above are final-stage convenience runners and require
   an already generated `images/chosenfreq/` anchor set.  A raw-image rebuild
   must first execute the four modes in their historical order through the
   isolation-gated `run_continuum_four_pass.py`; see
   `analysis/multi_freq_from_archive/CONTINUUM_FOUR_PASS_WORKFLOW.md`.  The
   deterministic replacement for the missing consolidation notebook verifies
   every generic minimum-beam winner against the numeric product-and-coordinate
   rules in `chosen_frequency_component_policy.csv`, which are derived from the
   broader `s230_anchor_identity_audit.csv`.  It fails if a smaller-beam nearby
   component does not match the reviewed nuclear/component decision.

   The generic NGC 4945 anchor is the historically selected but rejected
   source-17 component.  The four-pass orchestrator therefore applies and
   checks the mandatory source-18 fixed-nucleus rerun before it reports
   completion; its uncorrected fixed-to-chosen table is not a final product.

   The earlier RMS-based corrected isolated workflow completed on 2026-08-18 (HST): all
   four passes, nine guarded anchors, the NGC 4945 source-18 apply/check, and
   the unchanged 455-input inventory guard passed.  Its 406 output rows exactly
   match the frozen 409-row historical raw reference after the three reviewed NGC 4258
   exclusions, and every output product path is release-relative.  The compact
   evidence and checksums are retained as pre-scaled-MAD audit history in
   `analysis/multi_freq_from_archive/CONTINUUM_FOUR_PASS_WORKFLOW.md` and
   `CORRECTION_CHANGELOG.md`.

   The active fitter uses ordinary two-decimal formatting and returns coherent
   ICRS and J2000 coordinate pairs.  For each ordinary final fit, it measures
   source-excluded local scatter in a full-width `10 Bmaj` square, adopts
   `1.482602218505602 * medabsdevmed` exactly once, requires at least 50
   effective independent beam samples, and uses fitted flux divided by that
   scaled MAD for S/N and classification.  Classic RMS, standard deviation,
   mean, and RMS/scaled-MAD disagreement remain diagnostics rather than
   row-dependent alternate estimators; 12 predesignated dirty-SMA failures keep
   their separate one-beam RMS upper-limit policy.  A nonconverged forced-SMA
   component fit is recorded with a blank, unused fitted-flux diagnostic rather
   than a fabricated zero, while its independent three-RMS limit is retained.
   Search boxes are intentional full-width squares derived from `Bmaj`, not
   beam ellipses; the exact stage-specific widths are recorded in the four-pass
   workflow document and implemented by the tested
   `square_bmaj_region_width_arcsec()` helper.
6. Adjudicate and promote the final fitter table and NGC 4945 correction:
   the isolation-gated four-pass orchestrator runs the fixed-source-18 apply and
   check itself, validates the exact 45-row photometry and seven-row component
   audits, and then stops at `review_required` when the scaled-MAD layer changes
   flux or S/N fields.  Accept only the exact product-complete reviewed-delta
   ledger through the evidence-only `adjudicate` command.  Promote the
   adjudicated 495-row sorted table and its noise/source-18/audit sidecars; do
   not rerun `combine_fitsummaries.py`, which would reconstruct a legacy path
   rather than preserve the reviewed table.  The source-specific fit fixes the
   centroid at the source-18 AGN/maser nucleus from Emig et al. (2020), fixes
   the Gaussian axes and position angle to each restoring beam, and uses the
   same local scaled-MAD publication-noise policy as other ordinary final fits.
   Copy the following four run-evidence files byte-for-byte, without editing or
   re-rendering them, to these exact canonical release paths:
   `analysis/multi_freq_from_archive/reviewed_current_to_generated_deltas.csv`,
   `analysis/multi_freq_from_archive/reviewed_new_product_additions.csv`,
   `analysis/multi_freq_from_archive/reviewed_local_scaled_mad_freeze.json`, and
   `analysis/multi_freq_from_archive/reviewed_delta_adjudication.json`.  After
   the sorted table and named sidecars are promoted, require:
   ```bash
   python3 analysis/multi_freq_from_archive/reviewed_release_evidence.py --check
   ```
   This command deliberately fails when all four files are absent, when only a
   subset is present, or when any schema, hash binding, final count, or strict
   ledger recomputation differs.
7. Build the corrected machine-readable continuum table:
   ```bash
   cd analysis/multi_freq_from_archive
   python3 build_fitsummary_machine.py
   ```
   This retains the later 2009 `BB0091` reduction as the sole IC 2560
   measurement, preserves the older `BB091` reduction in the generated alias
   sidecar and excludes three failed NGC 4258 NVAS products.  The 15.36 and
   22.46 GHz failed reductions are replaced by the 15.365 and 22.460 GHz
   reprocessed detections from Doi et al. (2013), Tables 1--2; the
   artifact-dominated 2004-05-22 AH0847 image has no replacement row.  Their
   filenames, SHA-256 hashes, headers, fit/image diagnostics, former rows, and
   decisions remain in the generated exclusion sidecar.  The builder also
   verifies the exact 21-source analysis sample and the additional NGC 5495
   catalog entry.
   The 84 records in `validated_continuum_corrections.csv` preserve cumulative
   formatter/frame/earlier-noise and source-18 lineage, including the 66-row
   preapplication history and 18 newly admitted NGC 4945 products; they are not the complete
   final scaled-MAD change inventory.  The promoted reviewed-delta and addition
   ledgers must exactly own the final change scope.  SMA product
   provenance is stored as repository-relative `analysis/sma_calibration/...`
   paths in both the sorted and final fit summaries; recognized legacy absolute
   paths are normalized during rebuilding and unknown absolute paths fail.
8. From the release root, build and validate the audited black-hole parameter
   table:
   ```bash
   python3 analysis/CountingBHS_and_Spectra/build_bhs_params.py
   python3 analysis/CountingBHS_and_Spectra/build_bhs_params.py --check
   ```
   The production table contains six paired-posterior treatments, thirteen
   shared-distance cancellations, no independent-error cases, and two sources
   without a robust primary mass.  Statistical, systematic, and model terms
   remain separate; a quadrature total is supplied only where its components
   are numerical.  For NGC 4258, the paired-chain quantiles and 1 percent
   non-Keplerian allowance remain the adopted inputs.  The builder separately
   derives a `0.5368--0.5489 microarcsec` data-selection range from the 15
   accepted van der Boom & Kushnir (2026) configurations recorded in
   `bhs_model_sensitivity_checks.csv`.  That range is explicitly
   non-probabilistic and not adopted or combined with the statistical interval,
   the non-Keplerian allowance, or their quadrature total.

   For NGC 2960, the sidecar records the Kuo et al. (2018) fiducial and
   density-constrained linear and Herrnstein disk configurations.  Their exact
   `0.017504--0.017575 microarcsec` one-sided envelope is a maximum 0.404
   percent downward model sensitivity, not a probability interval.  The two
   unconstrained configurations are explicitly excluded because their inferred
   densities would quench masing.  This range is not combined with the adopted
   fixed-distance formal interval, and the latter remains incomplete because
   the paper does not document whether the earlier absolute-coordinate term was
   resolved.

   For NGC 6323, the same sidecar records the Kuo et al. (2018) fiducial
   configuration and the two physically allowed density-constrained disk
   models from the predecessor fit.  Configurations 1, 3, and 5 define the
   separate non-probabilistic `0.009460--0.009537 microarcsec` envelope;
   no-density configurations 2 and 4 are excluded.  The envelope lies within
   the adopted paired result's separate 1 percent allowance and is not
   transferred into, combined with, or plotted as a current-primary error.

   For Mrk 1029, the sidecar records Deepshikha et al.'s (2024; DOI
   `10.1017/S1743921323002491`) one-dimensional, edge-on, coplanar black-hole
   plus Mestel-disk model.  It reports `M_BH=(0.57 +/- 0.16)e6 solar masses`
   and `M_disk=(0.93 +/- 0.65)e6 solar masses` at an imported
   `126.787 +/- 0.059 Mpc`.  Holding that distance fixed maps the reported
   black-hole bounds to the standalone, non-probabilistic
   `0.0003317155--0.0005906155 microarcsec` span.  The exact spot selection,
   source-specific distance provenance, and confidence construction are
   unreported.  With no robust primary, the span has no directional comparison
   to an adopted result and is not adopted, combined, or plotted.

   NGC 1320 likewise remains a no-robust-primary case.  Its attached
   Deepshikha et al. (2024) record is a one-dimensional, edge-on, coplanar
   black-hole plus Mestel-disk diagnostic with `alpha=0.479 +/- 0.078`,
   `M_BH=(4.71 +/- 0.58)e6 solar masses`, and
   `M_disk=(0.10 +/- 1.14)e6 solar masses` at an imported
   `D=38.4 +/- 0.44 Mpc`.  Holding the reported distance fixed gives the
   standalone, non-probabilistic `0.01103254--0.01413127 microarcsec`
   black-hole-only span.  It has no directional comparison to a primary and
   is not adopted, combined, or plotted.  The separate audit-only sidecar
   preserves Kuo et al.'s (2020, Appendix A; DOI
   `10.1093/mnras/staa2260`) conditional analytic 1 percent turbulent-support
   correction as unattached evidence.  It is not an adopted uncertainty and
   cannot create a primary mass, shadow size, or figure point.

   For Circinus, NGC 3079, NGC 4388, and NGC 4945, the same sidecar records the Kuo et al.
   (2020, Appendix A) estimates that turbulent support can make the enclosed
   masses one-sided 3, 34, 0.1, and 24 percent underestimates.  The builder maps
   those conditional factors to `0.0415--0.0428`, `0.0119--0.0159`,
   `0.04535--0.04540`, and `0.0388--0.0481 microarcsec`, respectively.  Circinus retains its adopted
   `0.04152 microarcsec` center and radius-only `0.00078 microarcsec` numerical
   component; NGC 3079 retains its `0.0119 microarcsec` approximate enclosed-mass
   proxy with no formal interval, an explicitly unquantified order-unity model
   caveat, and no complete physical-mass total; NGC 4945 retains its
   `0.0388 microarcsec` center with no formal error; and NGC 4388 retains its publication display
   `0.0454 +/- 0.0011` for the fitting-plus-position result, with exact values
   machine-readable and the sparse-disk model caveat unquantified.
   These ranges are non-probabilistic, are not
   symmetrized or combined, and are not plotted as error bars.  The sidecar
   records the Keplerian Circinus/NGC 4388/NGC 4945 branches and flat-rotation
   NGC 3079 branch, together with their source-specific radii and turbulent
   velocities.

   The attached discriminated-union sidecar loader fails closed on an
   unexpected schema, unknown source or semantic status, fields populated for
   the wrong check type, malformed or unordered endpoint values or factors,
   missing DOI provenance, or multiple attached checks for one source without
   an explicit aggregation policy.  The separate audit-only loader freezes its
   sole NGC 1320 record and rejects any changed identity, attachment status,
   factor, physical assumption, or provenance.  Regression tests freeze the
   adopted NGC 4258
   quantiles, allowance, and total; the adopted Circinus center and radius-only
   component; the adopted NGC 3079 center, no-formal statistical status,
   unquantified-model/incomplete-total statuses, reported legacy Tully--Fisher
   distance label, and blank physical-mass total; the adopted NGC 4945 center
   and no-formal-error status; the adopted NGC 4388 fitting-plus-position component; and the
   adopted NGC 2960 center, formal interval, unquantified-astrometric status,
   and unavailable complete total; and the adopted NGC 6323 paired quantiles,
   one-percent allowance, and total, independently of their sidecar ranges.
   They also freeze Mrk 1029's blank primary mass and shadow fields while
   allowing only its two separately typed, non-adopted checks, and freeze the
   analogous NGC 1320 blanks while attaching only its Deepshikha span and
   proving that the Kuo audit record never reaches generated output.

   The separate `bhs_nonadopted_ratio_fits.csv` sidecar records the reviewed
   conditional $M/D$ fits for NGC 1194, NGC 2273, NGC 2960, NGC 4388, NGC 5765b,
   CGCG 074-064, NGC 6264, UGC 3789, NGC 6323, ESO 558-G009, J0437+2456, Mrk 1029, and NGC 1320 from Villaraos et
   al. (2022) and González-Juárez et al. (2024).  Its loader requires exactly those thirteen
   reviewed rows.  For NGC 1194, NGC 2273, NGC 2960, NGC 4388, NGC 5765b,
   CGCG 074-064, NGC 6264, UGC 3789, NGC 6323, ESO 558-G009, and J0437+2456 it preserves the adopted primary results while
   attaching conditional comparisons; the NGC 6264 row is explicitly an
   older-partial-data reanalysis rather than a same-data refit, while UGC 3789
   uses only the Reid et al. VLBI spot subset and omits the systemic-maser
   acceleration information used by the primary distance fit.  The NGC 2960
   row likewise reuses the older Kuo et al. (2011) spots, fixes the inclination
   and vertical center, selects a 15 degree scatter angle, and records the
   source's additional up-to-4 percent vertical-center sensitivity outside its
   formal interval.  Its adopted Kuo et al. (2018) result instead retains only
   the published fixed-distance formal interval; because that paper does not
   document whether the legacy absolute-coordinate term was resolved, the
   residual astrometric systematic is explicitly unquantified and no complete
   shadow or physical-mass total is claimed; the displayed mass bounds combine
   only the available fit and adopted-distance components.  The NGC 6323 row
   reuses the older Kuo et al. (2011) spots, fixes the Kuo et al. (2015)
   inclination, selects an 11 degree scatter angle, and fits four parameters
   including the vertical center.  Its
   `0.009396 +/- 0.000051 microarcsec` result is 1.61 percent below the adopted
   paired-chain median and remains non-adopted, uncombined, and unplotted; the
   provenance also preserves the Table 2/Section 4.1.2 reduced-chi-squared
   discrepancy.  The J0437+2456 row reuses the Gao et al. (2017) spots, fixes
   the earlier 81-degree inclination and vertical center, uses zero added
   scatter, and gives `0.004360 (-0.000123/+0.000103) microarcsec`, 4.02 percent
   above the retained Pesce et al. (2021) primary.  It remains non-adopted,
   uncombined, and unplotted; the sidecar also records that changing the fixed
   vertical center moves the restricted fit by about 2 percent.  For the other
   two sources it preserves the absence of robust primary masses.  The CGCG row
   also freezes Villaraos et al.'s
   internally inconsistent scatter-angle/$\chi^2$ descriptions as provenance
   without adopting either value.  Regression tests verify that no
   literature-fit field can populate or modify an adopted mass, shadow size,
   or uncertainty.

   Mrk 1029 remains a no-robust-primary case.  Gao et al. (2017) describe two
   plausible velocity systems, `9076 +/- 32` and `9160 +/- 61 km/s`; the audit
   selects the Huchra `9076 +/- 32 km/s` value used in Gao's distance
   construction, giving the exact measurement component `32/73.9 Mpc` and
   `D=120.8*73/73.9=119.3288 Mpc`.  Gao's `i=79 +/- 2 degrees` follows an
   indicative thin-disk classification of the spots and is a fixed geometric
   estimate, not a fitted posterior.  González-Juárez et al. (2024) reuse
   Gao's high-velocity spots in a restricted circular model, fix `i=79 degrees`
   and the vertical center at zero, and tune a 14-degree scattering angle; their
   `M/D=14400 +/- 1100 solar masses/Mpc` result has reduced chi-squared 1.385,
   while changing the position angle changes the ratio by 1.19 percent.  It
   remains a conditional, non-adopted ratio check.  Deepshikha et al.'s
   standalone Mestel-disk check described above remains non-probabilistic.
   Neither check is adopted, combined, or plotted.  No later superseding full
   Mrk 1029 fit was found through 2026-08-20; the later WISDOM tables are
   compilations of Gao's indicative values rather than new fits.

   NGC 1320 also remains a no-robust-primary case.  Gao et al. (2017) report
   two variable redshifted/blueshifted maser groups, no systemic VLBI masers,
   and several spots away from the apparent major axis, allowing a strong warp
   or non-disk/outflow interpretation.  The project distance record uses Gao's
   selected `V_LSR=2649 +/- 16 km/s` value, so its velocity-derived component
   is exactly `16/73.9 Mpc`.  Its inclination is reported as unconstrained;
   `i=90 degrees` is a fixed assumption in the later restricted fit, not a
   measurement.  González-Juárez et al. (2024) reuse the same Gao red/blue
   spots, fit `M/D`, `z_rec`, `x0`, and `y0`, fix `i=90 degrees`, select a
   17-degree scattering angle, and obtain
   `M/D=149700 (-6900/+7100) solar masses/Mpc`, reduced chi-squared 1.127,
   and 2.64 percent position-angle sensitivity.  The corresponding conditional
   shadow is `0.015356 (+0.000728/-0.000708) microarcseconds`.  This restricted
   result does not resolve the missing-systemic, off-axis, warp/outflow
   limitations and remains non-adopted, uncombined, and unplotted.  The
   Deepshikha and Kuo diagnostics described above likewise create no adopted
   value.  No superseding full NGC 1320 fit was found through 2026-08-20.

   The NGC 6323 publication row also uses the same Pesce et al. (2020b) joint
   fit's directed inclination intercept, `(0):91.5 +/- 0.3 degrees`, replacing
   the older Kuo et al. value without changing the adopted mass, distance,
   paired shadow quantiles, one-percent allowance, or BHS-figure point.  The
   restricted `M/D=(0.901 +/- 0.005)e5 solar masses/Mpc` result in
   arXiv:2601.00989v2 is an unrefereed, audit-only comparison with a
   prior-dependent joint table and is not a release sidecar row or a
   superseding primary fit.  No later peer-reviewed full NGC 6323 mass,
   distance, acceleration-plus-VLBI disk, and inclination fit was found through
   2026-08-20.

   J0437+2456 uses the Pesce et al. (2021) nine-parameter fit to the Gao et al.
   VLBI data.  The fixed-distance mass interval gives
   `0.004191 +/- 0.000293 microarcsec`; the exact `40/73.9 Mpc` velocity term
   affects only the distance and physical-mass diagnostics because the shared
   scale cancels from $M/D$.  The fit uses only high-velocity masers, has no
   acceleration information, and leaves the inclination posterior
   unconstrained, so its model contribution and complete shadow and
   physical-mass totals remain unavailable.  Table 1 consequently reports
   `unconstr.` from the adopted Pesce fit rather than carrying forward Gao's
   older geometric estimate.  Villaraos et al. (2022) give the separate
   same-data, fixed-inclination/fixed-center result recorded above.  A short
   Deepshikha et al. (2024; DOI `10.1017/S1743921323002491`) proceedings/poster
   analysis reports a much lower `1.85 +/- 0.23 million solar masses` from a
   one-dimensional, edge-on Mestel-disk decomposition at an external
   `68.27 +/- 0.45 Mpc`.  Because it supplies neither a full astrometric disk
   likelihood nor a documented mass--distance covariance and even yields
   nonphysical negative disk masses for several sample objects, it is retained
   as a separate, non-probabilistic `0.002434--0.003125 microarcsec`
   model-sensitivity check rather than an adopted comparison or probability
   interval.  It does not define a complete model envelope and does not change
   the unquantified-model status.  No published full acceleration-plus-VLBI refit was
   found through 2026-08-20.

   ESO 558-G009 remains a Gao et al. (2017) shared-distance result.  Its
   physical-mass display combines the fixed-distance fitting-plus-position and
   adopted-distance terms in the base interval and shows separately the 6.003
   percent allowance constructed from Gao's generic 5.2 percent midline and 3
   percent inclination terms as `geom. sys.`.  The shadow-size total combines
   that geometry allowance with the separate 2.3 percent fit-plus-position
   component; the shared distance cancels from $M/D$.  Table 1 reports Gao's
   estimated `98 +/- 1 degrees` inclination rather than `e.o.`.  Villaraos et
   al. reuse the same VLBI spots, fix the inclination at 98 degrees, and select
   a 7 degree scatter angle, so their 0.89-percent-higher ratio is retained as
   a conditional same-data check rather than an independent or superseding
   measurement.

   UGC 6093 remains a Zhao et al. (2018) shared-distance result.  Its physical
   mass includes the formal fit, the adopted-distance contribution, and the
   source-reported up-to-5 percent absolute-reference-position term; the Table
   1 renderer displays the position contribution separately as `pos. sys.`.
   Only the shared distance and `H0` factors cancel from `M/D`, so the position
   term remains in the shadow budget.  González-Juárez et al. (2024) report no
   usable UGC 6093 fit because relative spot-position errors rival the angular
   offsets.  That distinct limitation creates no numerical alternate-fit row
   and does not add a second position systematic.

   A final release must also pass the explicit pending-author gate:
   ```bash
   python3 analysis/CountingBHS_and_Spectra/build_bhs_params.py --check --require-final
   python3 analysis/multi_freq_from_archive/export_bhs_params.py
   ```
   Both commands now exit with status 0 after the explicit IC 2560 formal-only
   decision; any future input carrying a pending-author status will still stop
   the release.
9. From the release root, render and check the provisional manuscript BHS table
   and publication figures from the audited machine products.  Set the manuscript directory
   explicitly, build into provisional locations, compare against the current
   manuscript copies, and copy only after reviewing the differences:
   ```bash
   AE8A42_MANUSCRIPT_DIR=/path/to/revision_package/working_clean

   python3 analysis/CountingBHS_and_Spectra/build_bhs_latex_table.py \
     --output /tmp/BHS_Params.provisional.tex
   cmp /tmp/BHS_Params.provisional.tex \
     "$AE8A42_MANUSCRIPT_DIR/BHS_Params.tex"

   python3 analysis/multi_freq_from_archive/build_continuum_latex_excerpt.py \
     --output /tmp/Continuum_from_Pipeline.provisional.tex
   cmp /tmp/Continuum_from_Pipeline.provisional.tex \
     "$AE8A42_MANUSCRIPT_DIR/Continuum_from_Pipeline.tex"

   python3 analysis/multi_freq_from_archive/build_bhs_flux_figure.py \
     --output-dir analysis/multi_freq_from_archive
   cmp analysis/multi_freq_from_archive/BHS_Flux_Graph.png \
     "$AE8A42_MANUSCRIPT_DIR/BHS_Flux_Graph.png"
   cmp analysis/multi_freq_from_archive/BHS_Flux_Graph.pdf \
     "$AE8A42_MANUSCRIPT_DIR/BHS_Flux_Graph.pdf"

   python3 analysis/multi_freq_from_archive/build_sed_figures.py \
     --output-dir /tmp/AE8A42_SpectraGraphs
   for source in NGC4258 NGC1194 NGC1068 NGC3079 NGC4945 Circinus; do
     cmp "/tmp/AE8A42_SpectraGraphs/${source}_Spectra.png" \
       "$AE8A42_MANUSCRIPT_DIR/${source}_Spectra.png"
   done

   python3 analysis/multi_freq_from_archive/build_contamination_tables.py
   cmp analysis/multi_freq_from_archive/machinetables/contamination/latex/contam_NGC4258.tex \
     "$AE8A42_MANUSCRIPT_DIR/contam_NGC4258.tex"
   cmp analysis/multi_freq_from_archive/machinetables/contamination/latex/contam_NGC1194.tex \
     "$AE8A42_MANUSCRIPT_DIR/contam_NGC1194.tex"
   cmp analysis/multi_freq_from_archive/machinetables/contamination/latex/contam_NGC1068.tex \
     "$AE8A42_MANUSCRIPT_DIR/contam_NGC1068_new.tex"
   cmp analysis/multi_freq_from_archive/machinetables/contamination/latex/contam_NGC3079.tex \
     "$AE8A42_MANUSCRIPT_DIR/contam_NGC3079_new.tex"
   cmp analysis/multi_freq_from_archive/machinetables/contamination/latex/contam_NGC4945.tex \
     "$AE8A42_MANUSCRIPT_DIR/contam_NGC4945.tex"
   cmp analysis/multi_freq_from_archive/machinetables/contamination/latex/contam_Circinus.tex \
     "$AE8A42_MANUSCRIPT_DIR/contam_Circinus_new.tex"
   ```
   After reviewing the generated files and any nonzero `cmp` results, update the
   manuscript copies explicitly:
   ```bash
   cp /tmp/BHS_Params.provisional.tex \
     "$AE8A42_MANUSCRIPT_DIR/BHS_Params.tex"
   cp /tmp/Continuum_from_Pipeline.provisional.tex \
     "$AE8A42_MANUSCRIPT_DIR/Continuum_from_Pipeline.tex"
   cp analysis/multi_freq_from_archive/BHS_Flux_Graph.png \
     analysis/multi_freq_from_archive/BHS_Flux_Graph.pdf \
     "$AE8A42_MANUSCRIPT_DIR/"
   cp /tmp/AE8A42_SpectraGraphs/*_Spectra.png \
     "$AE8A42_MANUSCRIPT_DIR/"
   cp analysis/multi_freq_from_archive/machinetables/contamination/latex/contam_NGC4258.tex \
     "$AE8A42_MANUSCRIPT_DIR/contam_NGC4258.tex"
   cp analysis/multi_freq_from_archive/machinetables/contamination/latex/contam_NGC1194.tex \
     "$AE8A42_MANUSCRIPT_DIR/contam_NGC1194.tex"
   cp analysis/multi_freq_from_archive/machinetables/contamination/latex/contam_NGC1068.tex \
     "$AE8A42_MANUSCRIPT_DIR/contam_NGC1068_new.tex"
   cp analysis/multi_freq_from_archive/machinetables/contamination/latex/contam_NGC3079.tex \
     "$AE8A42_MANUSCRIPT_DIR/contam_NGC3079_new.tex"
   cp analysis/multi_freq_from_archive/machinetables/contamination/latex/contam_NGC4945.tex \
     "$AE8A42_MANUSCRIPT_DIR/contam_NGC4945.tex"
   cp analysis/multi_freq_from_archive/machinetables/contamination/latex/contam_Circinus.tex \
     "$AE8A42_MANUSCRIPT_DIR/contam_Circinus_new.tex"
   ```
   A nonzero `cmp` result requires review: it can reflect changed content,
   rendering-library differences, or metadata, and does not by itself prove
   that either copy is stale. For the validated rc1 package, all six SED PNGs match
   the accepted files exactly; the BHS PNG differs in 80 pixels despite matching
   plotted data and reviewed appearance. No manuscript copy is changed by the
   rc1 validation. Use `scripts/render_publication.py` and the separate
   `RENDERING_ENVIRONMENT.txt` for the tested renderer.
   The IC 2560 decision is now canonical; run the final-release gate in step 8
   and require status 0 before preparing any release package.
   The NGC 4258 data-selection range appears in the generated Table 1 note,
   not in its adopted uncertainty cell.  The BHS-versus-flux figure omits
   horizontal BHS-uncertainty and vertical flux-uncertainty bars, as in the
   accepted manuscript. The uncertainties remain available in the tables;
   the non-probabilistic
   sidecar range is not used as a probability interval.
10. Generate the 230 GHz reference-flux comparison:
   ```bash
   cd analysis/multi_freq_from_archive
   python3 s230_reference.py \
     --source-text fitsumfiles/final/fitsummary_final_withextras.txt \
     --output machinetables/S230_extrapolation_comparison.csv
   ```
11. Generate the deterministic variability diagnostics:
   ```bash
   python3 build_variability_tables.py
   ```
12. Generate the deterministic contamination products:
   ```bash
   python3 build_contamination_tables.py
   ```

The contamination builder writes its CSV, JSON, and LaTeX products under
`machinetables/contamination/`.  For NGC 1068 it writes the canonical
`latex/contam_NGC1068.tex` table from the source-reported 21.9 GHz total HSA S1
model flux of `7 +/- 1 mJy` in Gallimore et al. (2024), together with a clearly
non-adopted `latex/contam_NGC1068_central_only.tex` diagnostic that preserves
the historical image-peak calculation for comparison.

The SMA continuum fitter reads
`analysis/sma_calibration/selected_clean_iterations.csv`. Paths beginning with
`analysis/` are resolved from the repository root; paths in a relocated custom
manifest may instead be relative to that manifest. If a manifest-selected clean
image is unavailable, fitting stops with an error rather than silently replacing
it with an `nitercal_0` image; a missing manifest also stops the run. An
`nitercal_0` image is used only for a source
without a selected clean image and is recorded as an upper limit.

Only `run_continuum_fits_clean.sh` enables the selected-clean manifest.  The
dirty-image runner continues to fit its explicitly supplied `nitercal_0` path
and cannot be redirected to a clean image by the manifest.

`Beam_Date_Flux_Graphs.ipynb`, `Variability_Search.ipynb`, `get_contams.ipynb`,
and `flux_size_graph.ipynb` are retained only as historical diagnostic records.
Their plotting cells write to notebook-specific diagnostic directories and
cannot delete or overwrite the canonical `SpectraGraphs` publication directory;
any other notebook outputs retain legacy, noncanonical names.  The deterministic
builders above, rather than the historical notebook cells, are the publication
table and figure paths.

Machine-readable tables retain full numerical precision.  Publication-facing
generators generally show uncertainty components to two significant figures,
match central values to the finest displayed uncertainty place, and show
derived values without uncertainties to three significant figures.

## Notes on Machine Tables

`fitsummary_machine.txt` gives the full continuum-fitting summary table.  The
companion `continuum_recovery_manifest.csv` is the authoritative portable
recovery crosswalk.  It covers all 501 table rows, with verified sizes and
hashes for 222 ALMA FITS products, 248 NVAS images, and 15 SMA CASA image
directories, plus an explicit literature-only disposition for 16 rows; the 485
local products total 6,040,874,990 payload bytes.  The early ALMA delivery
filenames without a modern `member.uid___...` prefix are mapped to reviewed
Member OUS and execution-block identities rather than inferred from filename
shape.  NGC 5495 contributes three additional catalog measurements and is
excluded from calculations for the stated 21-source sample. See
`DATA_ACCESS.md`.

The preapplication continuum correction snapshot contains 66 historical records
and changed 50 displayed flux fields and 54 displayed S/N fields without
changing a classification within that earlier fit-rerun scope.  The cumulative
manifest contains 84 rows after adding 18 newly admitted NGC 4945 products. The final source-excluded
scaled-MAD classifications and 89 admissions are instead controlled by the
separate exact reviewed-delta and addition ledgers.  Separately,
two failed NGC 4258 NVAS rows formerly presented as upper limits are removed
and replaced by explicitly cited Doi et al. (2013) detections.  A third,
artifact-dominated 2004-05-22 AH0847 row formerly reported as a 79.89 mJy
detection is removed without replacement.  In addition to
the formatter/frame/noise
corrections, all 45 NGC 4945 rows are now measured at source 18 rather than the
source-17 SSC/SNR-associated peak.  Its 349.40 GHz anchor changes from 15.8 to
13.1 mJy and its extrapolated `S_230` from 21.2 to 17.5 mJy.  NGC 4945
remains among the five selected sources above 10 mJy.

The diagnostic `variability_machine.txt` and `correlations_machine.txt`
products are included for reproducibility but are not rendered as paper
tables.  They are generated directly by `build_variability_tables.py`; the
builder uses `continuum_source_roles.csv` to restrict the variability
calculations to the 21-source analysis sample while retaining the NGC 5495
catalog rows.  The final 498 sample measurements produce 14,193 unordered
within-source pairs and a 105-row strict screen, including 103 different-day
and two same-day comparisons.  The legacy
`clean_variability_tables.py` repair utility is not part of the corrected
workflow.

Run the complete source regression suites in their documented environments:

```bash
python3 scripts/test_release.py \
  --fitting-python /path/to/new-wmagn-runtime/wm-agn-casa-python \
  --renderer-rootfs /path/to/ubuntu2404-render-rootfs \
  --output /path/to/new-test-output
```

This inventories both test directories, routes the two figure-artifact modules
to the locked renderer, and treats every nonzero test exit or source change
as a failure. It does not waive stale release evidence. These source tests
are separate from the actual multi-hour fitting run and scientific comparison;
neither substitutes for the other. See `INSTALL.md` and the original rc1
test report for results, intentional absent-private-document skips and limits.
