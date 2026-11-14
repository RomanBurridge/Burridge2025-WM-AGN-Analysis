# AE8A42 Software Correction Record

Working branch: `ae8a42-corrections`

Base release: Zenodo/Git tag `v1.0.8`, commit `9aa7f8b4645e53cce339fb308136b131139919d0`

Started: 2026-08-12

This file records the software corrections carried by Zenodo releases v1.0.9 and v1.0.10 accompanying the corrected paper. It distinguishes implemented and tested changes from external publication and delivery actions.

## Scope of the present software pass

- Deterministic black-hole mass, distance, shadow-size, and uncertainty
  calculations for the full 21-source sample.
- Deterministic generation of manuscript Table 1 and the BHS-size versus
  230 GHz flux figure.
- Continuum machine-table reproducibility and 230 GHz reference-flux metadata.
- Deterministic variability diagnostic tables.
- Deterministic 230 GHz thermal-dust and extended-jet contamination products.
- Tests and documentation for all workflows above.

## Authoritative Task 28 continuum closeout

- Audited the ALMA filename-selection rules and removed filename-only rejection
  of otherwise eligible `tt0`, `mfs`, `manual`, and `split.cal` continuum
  products.  The approved Task 28 inputs comprise 49 pre-cutoff omissions, 26
  post-cutoff target continuum datasets, and seven Stokes-I MFS products.
- Combined with the seven products admitted at the A8 checkpoint, the current
  addition policy contains 89 ALMA products relative to the immutable 409-row
  historical baseline.  The authoritative isolated run is
  `task28-authoritative-v2-230316-flag2-20260825`.
- The prepared inventory contains 544 inputs: 222 ALMA, 307 VLA, and 15 SMA.
  The reviewed raw publication table contains 495 rows (406 retained
  historical-common products plus 89 additions), and the complete diagnostic
  noise union contains 498 rows including three policy-bound boundary
  rejections.
- After the nine VLA bandwidth-smearing exclusions and one IC 2560 alias
  exclusion, 485 local products remain: 222 ALMA, 248 NVAS, and 15 SMA.
  Sixteen literature measurements bring the final continuum table and recovery
  manifest to 501 rows.  The 485 local payloads total 6,040,874,990 bytes.
- The 498 sample measurements generate 14,193 unordered within-source pairs
  and a 105-row strict screen, including 103 different-day and two same-day
  comparisons.  Three NGC 5495 rows remain additional catalog measurements.
- The final 19-source $S_{230}$ table contains 12 detections and seven upper
  limits selected from 11 ALMA, seven SMA, and one NOEMA measurement.  Five
  sources exceed 10 mJy: NGC 4258, NGC 1068, NGC 3079, NGC 4945, and Circinus.
- NGC 4945 fixed-nucleus photometry now covers 45 products while retaining the
  exact pre-Task-28 lineage of the 27 previously released products.  Seven
  legacy execution-block ALMA filenames are mapped to reviewed Member OUS and
  proposal identities in the recovery manifest.
- All continuum tables, $S_{230}$ values, variability products, contamination
  products, manuscript tables, and figures are regenerated from this state.
  The historical A8 and post-A8 checkpoints below remain audit history.

## Historical A8 continuum checkpoint (completed 2026-08-20 HST)

- An inclusive ALMA public-release census through 2025-01-27 found 112 Member
  OUS records absent from the historical input inventory.  Seven are eligible
  additions and 105 have documented exclusions.  The lightweight 112-row
  census, seven-row admission manifest, query/provenance record, and preserved
  pre-fix notebook are released under
  `analysis/multi_freq_from_archive/archive_inventory_audit/`; no archival FITS
  bytes are bundled there.
- Six eligible ESO products were missed because the legacy discovery notebook
  searched the raw product string for `SO`, which also matched the letters in
  `ESO`.  A seventh NGC 1068 product was independently omitted.  The shipped
  `get_urls.ipynb` now applies boundary-token chemical matching to the product
  basename.
- The final external four-pass inventory is 462 inputs: 140 ALMA, 307 VLA, and
  15 SMA. The exactly reviewed publication table contains 413 rows,
  comprising 406 retained historical-common rows and all seven publication-row
  additions. The complete generic diagnostic audit contains 416 rows: those
  413 products plus three exact, hash-bound, converged fits retained as noise
  diagnostics but excluded by the historical model-contour boundary gate.
  Downstream alias and literature assembly produced 425 measurements: 412
  retained local products (140 ALMA, 257 NVAS, and 15 SMA) plus 13
  literature-only measurements.  The final recovery payload is
  3,943,670,510 bytes.
- The final selected-anchor audit contains 19 sources: 12 detections and seven
  upper limits drawn from ten ALMA and nine SMA observations.  The admitted
  ESO 558-G009 X600 product is the selected reference measurement: `0.66 mJy` at
  `313.11 GHz` in a `0.71817 arcsec` beam, giving
  `S_230 = 0.8190724285773601 +/- 0.09125138464542781 mJy`.  Its Gao
  et al. VLA-B comparison coordinate is a nuclear proxy rather than an explicit
  disk center, and the J0648--1744 documentation lacks complete quantitative
  astrometric QA.
- The immutable
  `fitsumfiles/historical/fitsummary_fixtochosen_sorted_409_reference.txt`
  freezes the 409-row pre-scaled-MAD baseline and must remain distinct from the
  historical 413-row A8 table and current 495-row Task 28 table. Similarly, the
  run-specific 66-row preapplication correction snapshot is historical
  formatter/frame/RMS-era and NGC 4945 source-18 lineage, not the complete final
  scaled-MAD change set or the 84-row cumulative Task 28 lineage file.
  The run-specific `reviewed_current_to_generated_deltas.csv` and
  `reviewed_new_product_additions.csv`, bound by the exact review-freeze JSON,
  own the final delta and admission scope.
- The isolated schema-4 authoritative run is
  `cgcg-augmented7-local-scaled-mad-230316-flag2-authoritative8`.  All four CASA
  stages, the NGC 4945 source-18 correction, inventory and identity checks,
  no-clobber review rendering, schema-3 freeze, and schema-2 adjudication
  completed with exact reviewed matches.  The exact promoted hashes are recorded
  in the current-release fingerprint below.

## Historical post-A8 VLA bandwidth-smearing checkpoint (completed 2026-08-25 HST)

- Raw VLA FQ tables and the restoring-beam geometry were used to evaluate
  `x=(channel_width/center_frequency)*(source_offset/radial_beam_FWHM)` for
  every retained VLA measurement. The adopted limit is `x <= 0.5`.
- Nine rows exceed the limit and are excluded: four for NGC 3079, four for
  NGC 4258, and one for UGC 6093. The remaining 248 VLA catalog measurements
  pass. Exact inputs and decisions are frozen in
  `continuum_bandwidth_smearing_exclusion_policy.csv` and the generated
  `machinetables/vla_bandwidth_smearing_exclusions.csv` sidecar.
- The final continuum table now has 416 rows: 403 fitted products (140 ALMA,
  248 NVAS, and 15 SMA) plus 13 literature measurements. Its 413 sample rows
  generate 10,137 unordered pairs and a 95-row matched-pair screen.
- All 19 selected 200--400 GHz reference measurements and their `S_230`
  values are unchanged. The NGC 4258 extended-jet diagnostic shifts to
  `0.1402 +/- 0.0043 mJy`, still 1.1 percent of `S_230` after rounding. The
  NGC 3079 and NGC 4258 SEDs were regenerated; no source ranking or principal
  conclusion changes.
- NVAS processing logs confirm that primary-beam correction was already
  applied. No second primary-beam correction is introduced by this closeout.

## Changes implemented

Unless a paragraph explicitly describes the augmented architecture above, the
455/409/406/418/405 and 66-row counts in this section are historical
pre-augmentation checkpoints retained for audit continuity.

### Black-hole parameters, manuscript table, and figure

- Replaced the stateful black-hole notebook calculation with
  `build_bhs_params.py`, a deterministic 21-source builder driven by
  `bhs_source_inputs.csv`, `pesce_theta_g_quantiles.csv`, the attached
  `bhs_model_sensitivity_checks.csv`, the separately validated, unattached
  `bhs_audit_only_model_sensitivity_checks.csv`, and
  `bhs_nonadopted_ratio_fits.csv` provenance sidecars.
- Classified the sample into exactly six paired-posterior cases, thirteen
  shared-distance cancellations, zero independent mass--distance cases, and
  two sources without a robust primary mass or adopted shadow size.
- Used paired angular-gravitational-radius quantiles supplied by D. W. Pesce
  from the chains underlying the published results for NGC 4258,
  CGCG 074-064, UGC 3789, NGC 6264, NGC 6323, and NGC 5765b.  Their asymmetric
  statistical intervals and recommended 1 percent non-Keplerian allowances
  remain separate.
- Added a machine-readable NGC 4258 check from van der Boom & Kushnir (2026;
  DOI `10.1093/mnras/stag1327`).  The sidecar records that this is a
  component-tracking/data-selection reanalysis of a subset of the Argon et al.
  (2007) 1997--2000 VLBI epochs used by the primary result, identifies the 15
  accepted configurations and excluded configurations 2 and 17, preserves the
  published endpoint masses and distances, and derives the corresponding
  `0.5368--0.5489 microarcsec` shadow-size range with the same IAU-2015
  conversion used by the primary builder.
- Retained the adopted NGC 4258 paired-chain median and quantiles, its separate
  1 percent non-Keplerian allowance, and their existing quadrature display
  total without numerical change.  The new range is labeled
  non-probabilistic and `reported_separately_not_adopted`; it is attached only
  after the adopted row is calculated and never enters `size_unc_adopted_uas`,
  `size_total_unc_*`, or the adopted central BHS size.  The publication figure
  omits horizontal BHS-uncertainty bars and directs readers to Table 1 for the
  statistical and separate systematic/model terms.
- Added a second strictly typed model-sensitivity branch for Circinus from
  Kuo et al. (2020, Appendix A; DOI `10.1093/mnras/staa2260`).  Under its
  source-specific Keplerian assumptions (`r_in=0.11 pc`, `v_turb=30 km/s`,
  constant turbulent velocity, turbulent pressure dominating thermal pressure,
  and `Sigma proportional to r^-1`), the enclosed mass can be a one-sided
  3 percent underestimate.  This maps the unchanged `0.04152 microarcsec`
  center to a separate `0.0415--0.0428 microarcsec` check
  (`+0.0012 microarcsec`).  It is non-probabilistic, is not symmetrized,
  combined, or plotted, and does not replace the radius-only numerical
  component derived from Greenhill et al.'s angular-radius uncertainty.
- Added a third strictly typed model-sensitivity branch for NGC 3079 from
  Kuo et al. (2020, Appendix A; DOI `10.1093/mnras/staa2260`).  Under their
  wind-driven turbulence assumptions, their additive correction satisfies
  `Delta M_turb/M_enc=0.34`, or `M_BH=1.34 M_enc`, corresponding to
  `0.0119--0.0159 microarcsec` or a `+0.0040 microarcsec` shift.  This is a
  non-probabilistic model sensitivity: it is not symmetrized, combined with a
  formal error, or plotted as an error bar, and the adopted NGC 3079 center
  remains unchanged.  Kondratko et al.'s central value is explicitly typed as
  an approximate enclosed-mass black-hole proxy with no formal interval and an
  unquantified order-unity geometry/disk-mass caveat, so its complete shadow
  and physical-mass totals remain unavailable.  The adopted legacy
  Tully--Fisher uncertainty is labeled `reported`, not statistical.
- Added the corresponding fourth strictly typed NGC 4945 branch from the same Kuo
  et al. Appendix.  Its source-specific Keplerian model uses
  `r_in=0.16 pc`, `v_turb=60 km/s`, constant turbulent velocity, turbulent
  pressure dominating thermal pressure, and `Sigma proportional to r^-1` to
  give a one-sided 24 percent mass sensitivity.  This maps the unchanged
  `0.0388 microarcsec` center to a separate `0.0388--0.0481 microarcsec`
  check (`+0.0093 microarcsec`), which is likewise neither adopted,
  symmetrized, combined, nor plotted as an error bar.
- Added a fifth strictly typed branch for NGC 4388 from the same Kuo et al.
  Appendix.  Its source-specific Keplerian model uses `r_in=0.24 pc`,
  `v_turb=10 km/s`, constant turbulent velocity, turbulent pressure dominating
  thermal pressure, and `Sigma proportional to r^-1` to give a one-sided
  0.1 percent mass sensitivity.  This small but explicit check is retained
  as non-probabilistic and is not combined with or substituted for the adopted
  fitting-plus-position component.
- Added the missing NGC 3079 H2O line-VLBI provenance from Haschick et al.
  (1990) and Yamauchi et al. (2004), and clarified that Hagiwara et al. (2004)
  mapped OH absorption rather than the 22 GHz H2O maser line.  None supplies a
  formal replacement mass interval.
- For the thirteen fixed-distance mass normalizations, rescaled both the mass
  and its fixed-distance fitting interval to the adopted distance while
  canceling the shared distance factor exactly from `M/D`; the external
  distance uncertainty is not counted a second time in the shadow size.
- Corrected the Circinus uncertainty display inherited from R2.  R2 reported
  `0.042 +/- 0.013 microarcsec` after independently propagating the same
  distance scale carried by the mass.  After that shared distance cancels, the
  adopted center is `0.04152 microarcsec` and the surviving
  `+/-0.00078 microarcsec` component propagates only the published
  `5.3 +/- 0.1 mas` inner angular radius.  It is labeled radius-only rather
  than formal because Greenhill et al. quote no rotation-speed uncertainty or
  fitted inclination and describe the mass as a lower limit.  The logarithmic
  mass display now transforms Greenhill et al.'s printed
  `(1.7 +/- 0.3)e6 solar-mass` interval directly; the reconstructed
  `0.325e6 solar-mass` component total remains diagnostic only.
- Made two source-specific corrections that were obscured by the accepted R2
  table's single combined uncertainties and independent mass--distance
  propagation.  R2 reported UGC 6093 as `0.0184 +/- 0.0013 microarcsec` and
  ESO 558-G009 as `0.017 +/- 0.001 microarcsec`.  The revised UGC 6093 result
  separates the Zhao et al. (2018) formal interval from its 5 percent
  absolute-position systematic; the revised ESO 558-G009 result separates the
  Gao et al. (2017) 2.3 percent fitting-plus-position component from its
  6.0 percent combined midline/inclination model component.  The shared
  distance-dependent term cancels from `M/D` in both cases.
- Completed the ESO 558-G009 physical-mass and orientation audit.  The exact
  velocity-derived distance measurement term is now `27/73.9 = 0.3653586 Mpc`.
  The mass cell retains its fitting-plus-position and adopted-distance base
  interval, `7.225 (-0.025/+0.024)`, and separately displays the 6.003 percent
  allowance constructed from Gao et al.'s 5.2 percent midline and 3 percent
  inclination terms as `+/-0.026 dex (geom. sys.)`; the machine table retains
  the exact asymmetric logarithmic transforms and full quadrature total.  The
  shadow budget still combines that allowance with the distinct 2.3 percent
  fit-plus-position component, while the shared distance cancels.  Table 1 now
  gives Gao et al.'s estimated `98 +/- 1 degrees` inclination instead of
  `e.o.` and notes the weak disk-thickness constraint from the elongated
  north--south beam.
- Added Villaraos et al.'s ESO 558-G009 ratio as the tenth conditional
  non-adopted comparison.  It is 0.89 percent above the retained Gao ratio but
  reuses Gao's VLBI spots, fixes the inclination at 98 degrees, and selects a
  7 degree scatter angle; it is therefore labeled a same-data conditional
  check rather than independent confirmation or a superseding fit.
- Restored Zhao et al.'s source-reported, up-to-5 percent absolute-reference-
  position term to the UGC 6093 physical-mass budget.  The accepted working
  table had already propagated that term into the shadow-size uncertainty but
  omitted it from the mass fields.  The revised mass cell retains the formal
  fit plus adopted-distance interval, `7.388 (-0.026/+0.027)`, and displays the
  position term separately as `+/-0.022 dex (pos. sys.)`; the machine table
  also records the exact asymmetric log transforms and the full quadrature
  total.  Zhao et al.'s fitted recession-velocity interval is now retained
  asymmetrically rather than rounded to one symmetric distance term.  The
  shared distance and `H0` components still cancel only from `M/D`, so the
  adopted shadow size and BHS-figure point are unchanged.
- Documented that González-Juárez et al. (2024) report no usable UGC 6093 fit
  because relative per-spot position errors are comparable to the maser
  angular offsets.  This is distinct from Zhao et al.'s absolute reference-
  coordinate term and supplies neither a numerical alternate fit nor an
  additional systematic; no non-adopted-fit row is created.
- Restored the published NGC 3393 model sensitivity that was omitted from the
  accepted R2 display.  R2 reported only `0.063 +/- 0.005 microarcsec`; the
  revised table separates a `0.0636 +/- 0.0041 microarcsec` formal interval
  from the `-0.010/+0.0082 microarcsec` envelope produced by the alternate
  dynamical-center registration and rotation-law assumptions discussed by
  Kondratko et al. (2008).  This envelope is non-probabilistic and is neither
  symmetrized nor combined in quadrature with the formal term.
- Recorded separate mass, mass-uncertainty, and distance DOI provenance for
  every source, with `NA` used explicitly when no numerical mass uncertainty
  is adopted.  The calculation uses the IAU 2015 nominal solar mass parameter,
  the exact speed of light and astronomical unit, and the parsec definition.
- Kept paired statistical intervals, numerical systematics, reported
  undecomposed uncertainties, displayed model-sensitivity envelopes, separately
  reported non-probabilistic data-selection ranges, and unquantified model
  caveats semantically distinct.  No formal error is invented for NGC 4945 or
  NGC 3079; the Circinus angular-radius-only component is not presented as a
  complete formal uncertainty; and no primary shadow size is reported for
  Mrk 1029 or NGC 1320.
- Added a strictly validated non-adopted-fit sidecar for the conditional
  circular-disk $M/D$ fits reported for Mrk 1029 and NGC 1320 by
  González-Juárez et al. (2024; DOI `10.1051/0004-6361/202450098`).  The fits
  reuse the published Gao et al. VLBI positions rather than new observations;
  they are retained as literature checks but do not populate the adopted mass
  or shadow-size fields because their model does not resolve the complex maser
  geometry and dynamics that prevented robust primary masses.  The Mrk 1029
  row records the high-velocity-only
  `M/D=14400 +/- 1100 solar masses/Mpc` result, fixed `i=79 degrees`, fixed
  vertical center at zero, tuned 14-degree scattering angle, reduced
  chi-squared 1.385, and the reported 1.19 percent position-angle sensitivity.
- Added the González-Juárez et al. (2024) same-data NGC 1194 fit to that
  sidecar.  Its $M/D=(1.310\pm0.021)\times10^6~M_\odot\,\mathrm{Mpc}^{-1}$
  display value is 7.2 percent above the adopted Kuo normalization.  It remains
  a non-adopted literature comparison because it introduces no new VLBI
  measurements, fixes the inclination and vertical center, and gives a
  1-sigma recession-velocity interval that does not overlap the original
  value; the adopted Table 1 value and BHS-figure point do not change.
- Added the Villaraos et al. (2022; DOI `10.1093/mnras/stac2973`) same-data
  NGC 2273 fit to the non-adopted comparison sidecar.  Its symmetrized
  publication value, $M/D=(3.283\pm0.048)\times10^5~M_\odot\,
  \mathrm{Mpc}^{-1}$, is 12.5 percent above the adopted Kuo normalization;
  the exact asymmetric formal interval remains machine-readable.  The later
  result is conditional on a fixed vertical center and selected $15^\circ$
  scatter angle, so it is not combined with the adopted Gao-tabulated
  uncertainty and does not change the Table 1 row or BHS-figure point.
- Added the Villaraos et al. (2022) NGC 4388 fit to the same sidecar.  It reuses
  Kuo et al.'s five high-velocity spots and gives
  $M/D=4.233_{-0.268}^{+0.185}\times10^5~M_\odot\,\mathrm{Mpc}^{-1}$,
  or $0.04342_{-0.00275}^{+0.00190}$ microarcseconds, 4.3 percent below the
  adopted Kuo normalization.  Because the fit fixes the edge-on inclination
  and vertical center and does not resolve the absence of systemic masers or a
  demonstrated disk, it remains a non-adopted comparison and is not plotted.
- Added the Villaraos et al. (2022) NGC 5765b fit as the sixth reviewed
  non-adopted comparison.  It reuses the earlier Gao et al. VLBI data, fixes
  that fit's `94.5`-degree inclination, and selects a 10-degree scatter angle.
  Its $M/D=(3.727\pm0.013)\times10^5~M_\odot\,\mathrm{Mpc}^{-1}$ corresponds
  to $0.03823\pm0.00013$ microarcseconds, 0.8 percent above the adopted
  paired-chain median and within the separate 1 percent non-Keplerian
  allowance.  It is not adopted, combined, or plotted, and all primary
  NGC 5765b numerical fields remain unchanged.
- Added the Villaraos et al. (2022) CGCG 074-064 fit as the seventh reviewed
  non-adopted comparison.  Its Table 2 result,
  $M/D=2.749_{-0.015}^{+0.017}\times10^5~M_\odot\,\mathrm{Mpc}^{-1}$,
  corresponds to $0.028199_{-0.000154}^{+0.000174}$ microarcseconds, 0.5
  percent below the adopted paired-chain median.  It is not adopted, combined,
  or plotted.  Villaraos et al.'s Table 2 and Section 4.1.4 give inconsistent
  scatter-angle and reduced-$\chi^2$ values; the sidecar preserves that
  discrepancy as provenance without adopting either description.
- Added the Villaraos et al. (2022) NGC 6264 fit as the eighth reviewed
  non-adopted comparison.  It reuses the earlier Kuo et al. VLBI spot data and
  gives $M/D=(2.126\pm0.007)\times10^5~M_\odot\,\mathrm{Mpc}^{-1}$,
  corresponding to $0.021808\pm0.000072$ microarcseconds, 1.64 percent above
  the adopted paired-chain median.  It is conditional on a fixed vertical
  center and selected scatter angle and is therefore not adopted, combined,
  or plotted.  The sidecar explicitly records that its Table 2 and Section
  4.1.7 reduced-$\chi^2$ values are `1.304` and `1.301`, respectively.
- Added the Villaraos et al. (2022) UGC 3789 fit as the ninth reviewed
  non-adopted comparison.  It reuses the Reid et al. (2013) VLBI spot
  measurements but omits the systemic-maser acceleration information used by
  the primary geometric-distance fit, fixes the earlier Reid inclination, and
  selects a 10-degree scatter angle.  Its
  $M/D=(2.277\pm0.009)\times10^5~M_\odot\,\mathrm{Mpc}^{-1}$ corresponds to
  $0.023357\pm0.000092$ microarcseconds, 1.47 percent below the adopted
  paired-chain median.  It remains non-adopted, uncombined, and unplotted;
  Table 2 and Section 4.1.3 consistently report reduced $\chi^2=1.430$.
- Completed the NGC 2960/Mrk 1419 primary-error and model-check audit.  The
  adopted Kuo et al. (2018; DOI `10.3847/1538-4357/aabff1`) fiducial remains
  `0.017575 (-0.000085/+0.000099) microarcseconds`, with the numerical interval
  explicitly identified as a fixed-distance 68 percent marginalized formal
  term.  The fiducial and two density-constrained disk configurations define a
  separate, non-probabilistic `0.017504--0.017575 microarcsecond` envelope, a
  maximum 0.404 percent downward shift; configurations 2 and 4 are excluded
  because their inferred densities would quench masing.  Because the later
  paper does not document whether the legacy absolute-coordinate contribution
  was resolved, the residual astrometric systematic is marked unquantified, no
  complete shadow or physical-mass total is reported, and the earlier 3--5
  percent estimate is not transplanted as a new numerical error.  The displayed
  physical-mass interval combines only the available fit and adopted-distance
  components.
- Added the Villaraos et al. (2022) NGC 2960 result as the eleventh reviewed
  non-adopted comparison.  It reanalyzes the older Kuo et al. (2011) VLBI spots
  under fixed `i=89 degrees`, fixed vertical center, and a selected 15-degree
  scatter angle, giving
  `$M/D=1.738 (-0.016/+0.015)e5 solar masses/Mpc` or
  `0.017828 (-0.000164/+0.000154) microarcseconds`, 1.44 percent above the
  retained primary.  Its Table 2 and Section 4.1.6 give conflicting center and
  reduced-chi-squared descriptions, and varying the fixed vertical center
  changes `M/D` by up to 4 percent outside the formal interval.  These caveats
  are preserved in provenance; the comparison is not adopted, combined, or
  plotted, and exact spot-level overlap with the unpublished Kuo et al. (2018)
  input set is not claimed.
- Completed the NGC 6323 primary, orientation, and model-check audit.  The
  adopted Pesce et al. (2020b; DOI `10.3847/2041-8213/ab75f0`) paired-chain
  shadow remains
  `0.009549 (-0.000056/+0.000061) stat +/- 0.000095 systematic
  microarcseconds`.  Table 1 now uses that same fit's directed $r=0$
  inclination, `91.5 +/- 0.3 degrees`, instead of the older Kuo et al. value;
  no mass, distance, shadow, or plotted coordinate changes.
- Added NGC 6323 as the seventh strictly typed model-sensitivity record.  Kuo
  et al. (2018; DOI `10.3847/1538-4357/aabff1`) configurations 1, 3, and 5
  define the separate predecessor-fit, non-probabilistic
  `0.009460--0.009537 microarcsecond` envelope; the no-density configurations
  2 and 4 are excluded.  The range lies within the current primary's separate
  one-percent allowance and is not transferred, combined, or plotted.
- Added the Villaraos et al. (2022) NGC 6323 result as the twelfth reviewed
  non-adopted comparison.  It reanalyzes the older Kuo et al. (2011) VLBI
  spots, fixes the Kuo et al. (2015) inclination, selects an 11-degree scatter
  angle, and fits four parameters including the vertical center.  Its
  `$M/D=(0.916 +/- 0.005)e5 solar masses/Mpc` result corresponds to
  `0.009396 +/- 0.000051 microarcseconds`, 1.61 percent below the paired-chain
  median.  It is not adopted, combined, or plotted; the Table 2/Section 4.1.2
  reduced-chi-squared discrepancy is preserved in provenance.  The lower
  arXiv:2601.00989v2 value is an unrefereed audit-only check, not a sidecar row
  or superseding primary fit.
- Completed the J0437+2456 primary, distance-component, orientation, and
  comparison audit.  The Pesce et al. (2021; DOI
  `10.3847/1538-4357/abde3d`) nine-parameter high-velocity-only fit remains the
  primary, giving `0.004191 +/- 0.000293 microarcseconds` formally after the
  shared distance cancels.  The velocity-derived distance component now uses
  exact `40/73.9 Mpc` rather than a rounded decimal; this changes only
  distance and physical-mass diagnostics.  Because the fit's inclination
  posterior matches its broad prior, Table 1 now reports `unconstr.` with
  Pesce provenance instead of carrying forward Gao et al.'s approximate
  `81 +/- 1 degrees`.  The model contribution and complete shadow and
  physical-mass totals remain unquantified.
- Added the Villaraos et al. (2022) J0437+2456 result as the thirteenth
  reviewed non-adopted comparison.  It reuses the Gao et al. spots, fixes
  `i=81 degrees` and the vertical center, uses zero added scatter, and gives
  `$M/D=0.425 (-0.012/+0.010)e5 solar masses/Mpc`, or
  `0.004360 (-0.000123/+0.000103) microarcseconds`, 4.02 percent above the
  retained primary.  Its approximately 2 percent vertical-center sensitivity
  is preserved in provenance; the comparison is not adopted, combined, or
  plotted.
- At that historical checkpoint, added J0437+2456 as the eighth strictly typed
  model-sensitivity record.
  Deepshikha et al. (2024; DOI `10.1017/S1743921323002491`) report a
  one-dimensional, edge-on black-hole plus Mestel-disk decomposition with
  `M_BH=(1.85 +/- 0.23)e6 solar masses` at an imported 68.27 Mpc.  Holding the
  reported distance fixed gives the separate
  `0.002434--0.003125 microarcsecond` black-hole-only span.  The exact retained
  spot set and confidence semantics are unreported, so this range is
  non-probabilistic, is not a complete model envelope, and neither replaces
  nor augments the adopted formal interval.
- Completed the Mrk 1029 distance, geometry, and later-literature audit without
  creating an adopted primary.  Amid Gao et al.'s two published velocity
  systems, the canonical distance record selects the Huchra
  `9076 +/- 32 km/s` value used in the distance construction and now carries
  the exact `32/73.9 Mpc` measurement component.  Gao's
  `i=79 +/- 2 degrees` remains an indicative thin-disk geometry estimate, not a
  fitted inclination posterior.  No superseding full fit was found through
  2026-08-20; later WISDOM tables compile Gao's indicative values rather than
  refitting the source.
- Added Mrk 1029 as the ninth current strictly typed model-sensitivity record.
  Deepshikha et al. (2024; DOI `10.1017/S1743921323002491`) report a
  non-probabilistic, one-dimensional, edge-on black-hole plus Mestel-disk
  decomposition with `M_BH=(0.57 +/- 0.16)e6 solar masses` and
  `M_disk=(0.93 +/- 0.65)e6 solar masses` at an imported
  `126.787 +/- 0.059 Mpc`.  Holding that distance fixed gives the standalone
  `0.0003317155--0.0005906155 microarcsecond` black-hole-only span.  Exact spot
  selection, source-specific distance provenance, and confidence construction
  are unreported.  Because no robust primary exists, the span has no
  directional primary comparison; it is not adopted, combined, or plotted.
- Completed the NGC 1320 distance, geometry, and later-literature audit without
  creating an adopted primary.  Gao et al. (2017) report two highly variable
  redshifted/blueshifted maser groups, no systemic VLBI masers, and several
  spots away from the apparent major axis, allowing a strong warp or
  non-disk/outflow interpretation.  The canonical distance record selects
  Gao's `V_LSR=2649 +/- 16 km/s` value and now carries the exact `16/73.9 Mpc`
  measurement component.  Table 1 reports the inclination as `unconstr.`;
  `i=90 degrees` in the later restricted analyses is a fixed assumption, not a
  measured inclination.  No superseding full fit was found through 2026-08-20.
- Expanded the González-Juárez et al. (2024; DOI
  `10.1051/0004-6361/202450098`) NGC 1320 non-adopted ratio record.  It reuses
  the same Gao redshifted/blueshifted spots and fits `M/D`, `z_rec`, `x0`, and
  `y0`, while fixing `i=90 degrees` and selecting a 17-degree scattering
  angle.  Its `M/D=149700 (-6900/+7100) solar masses/Mpc` result corresponds
  to `0.015356 (+0.000728/-0.000708) microarcseconds`, with reduced
  chi-squared 1.127 and 2.64 percent position-angle sensitivity.  The
  restricted fit does not resolve Gao's missing-systemic, off-axis,
  warp/outflow limitations and remains non-adopted, uncombined, and unplotted.
- Added NGC 1320 as the tenth attached model-sensitivity record.  Deepshikha et
  al. (2024; DOI `10.1017/S1743921323002491`) report a non-probabilistic,
  one-dimensional, edge-on black-hole plus Mestel-disk diagnostic with
  `alpha=0.479 +/- 0.078`, `M_BH=(4.71 +/- 0.58)e6 solar masses`, and
  `M_disk=(0.10 +/- 1.14)e6 solar masses` at an imported
  `D=38.4 +/- 0.44 Mpc`.  Holding that distance fixed gives the standalone
  `0.01103254--0.01413127 microarcsecond` black-hole-only span.  With no robust
  primary, it has no directional comparison and is not adopted, combined, or
  plotted.
- Added the separately validated
  `bhs_audit_only_model_sensitivity_checks.csv` with exactly one record.  It
  preserves Kuo et al.'s (2020, Appendix A; DOI `10.1093/mnras/staa2260`)
  conditional analytic NGC 1320 turbulence result (`r_in=0.07 pc`,
  `v_turb=40 km/s`, `Sigma proportional to r^-1`, Keplerian rotation,
  `eta=2.5`, and ratio factors 1.00--1.01).  Because no robust primary exists,
  its attachment status is `audit_only_unattached_no_robust_primary`; it does
  not populate the generated row or become an adopted uncertainty, shadow
  size, or figure point.
- Audited the CGCG 074-064 disk-boundary calculation of Kuo et al. (2024;
  DOI `10.1093/mnras/stae1558`).  It is a physical-context model rather than a
  dynamical $M/D$ refit and implies a disk-to-black-hole mass ratio of only
  about 0.04 percent under its stated assumptions.  It is therefore recorded
  as qualitative context, not as a probability interval or an added shadow-size
  uncertainty.
- Updated the NGC 5765b Table 1 inclination from the earlier Gao et al. (2016)
  $r=0$ intercept, $94.50\pm0.25$ degrees, to the internally consistent
  Pesce et al. (2020b) joint-fit value, $72.4\pm0.5$ degrees.  The convention
  and `(0)` label are unchanged; the mapped warped disk remains near edge-on.
  This is an orientation-column correction only and does not change the mass,
  distance, shadow size, figure point, or conclusions.
- Updated the NGC 6264 Table 1 inclination from the older Kuo et al. (2013)
  value, $89.5\pm0.9$ degrees, to the internally consistent Pesce et al.
  (2020b) joint-fit value, $91.3\pm2.3$ degrees at $r=0$.  This is an
  orientation-only provenance correction: the adopted mass, distance, paired
  shadow-size quantiles, 1 percent allowance, and BHS-figure point are
  unchanged.
- Updated the UGC 3789 Table 1 inclination from Reid et al.'s value at
  $r=0.60$ mas, $90.6\pm0.4$ degrees, to the internally consistent Pesce et
  al. (2020b) joint-fit intercept, $84.9\pm0.6$ degrees at $r=0$.  Pesce's
  fitted inclination warp, $di/dr=7.7\pm1.0$ degrees per mas, gives a nominal
  value near $89.5$ degrees at 0.60 mas before covariance, so this is primarily
  a same-fit/reference-radius provenance correction rather than evidence for a
  physical 5.7-degree reorientation.  The adopted mass, distance, paired
  shadow-size quantiles, 1 percent allowance, and BHS-figure point are
  unchanged.
- Reviewed arXiv:2601.00989v2 as an unpublished restricted reanalysis of the
  same NGC 5765b and CGCG 074-064 data.  Its CGCG result,
  $M/D=(2.704\pm0.018)\times10^5~M_\odot\,\mathrm{Mpc}^{-1}$, corresponds to
  $0.02774\pm0.00018$ microarcseconds, 2.2 percent below the adopted median.
  These results are model- and prior-conditional and are neither treated as
  superseding published fits nor added to the adopted calculations.  The
  observing-log audit likewise found no later usable 22-GHz line-VLBI science
  dataset for either source.
- Corrected the NGC 4388 uncertainty semantics without changing its adopted
  central value: the `+/-0.2e6` solar-mass term is fixed-distance fitting plus
  source position, while the separate 11 percent Tully--Fisher calibration term
  is a shared systematic that cancels from $M/D$.  Table 1 now labels these
  `(fit+pos.)` and `(sys.)`, and reports the maser-disk inclination as
  unconstrained rather than treating the assumed edge-on geometry as a
  measurement.  Kuo et al.'s Table 4 value of 8.4 million solar masses is
  retained; the paper's 8.5-million narrative/Table 6 value is documented as a
  1.2 percent internal discrepancy within the available fixed-distance term.
- Corrected the MCP III Kuo et al. bibliography record from the online-first
  2010 date to the 2011 journal year used by ApJ volume 727 and the published
  citation; the internal BibTeX key remains unchanged for source stability.
- For NGC 1194 and NGC 2273, retained the conservative symmetric mass
  uncertainties tabulated by Gao et al. (2017), but replaced the unsupported
  `residual` label with `reported`.  Gao et al. do not decompose those two
  compilation values or specify a confidence construction; the central
  normalizations remain tied to Kuo et al. (2011), and no numerical value
  changes in this terminology correction.
- Added `build_bhs_latex_table.py` and `bhs_latex_metadata.csv` so manuscript
  Table 1 is rendered from `BHS_params.csv` without duplicating scientific
  numbers in presentation metadata.
- Added `build_bhs_flux_figure.py`, which reads only the audited BHS table and
  selected 230 GHz continuum table and deterministically writes
  `BHS_Flux_Graph.pdf` and `BHS_Flux_Graph.png`.  The figure distinguishes
  calibrated totals, formal-only intervals, model-sensitivity spans, and
  sources without formal size errors.  The separately reported NGC 4258
  data-selection range is not substituted for its unchanged adopted error bar.
- Converted `BH_Parameters.ipynb` into a readable diagnostic companion and
  prevented the historical size/flux notebook from overwriting the canonical
  table or publication figure.
- Implemented the formal IC 2560 fit term confirmed by N. Nakai.  No separately
  quantified model systematic is available, so IC 2560 is now
  explicitly classified as formal-only with an unquantified model contribution
  and an incomplete total interval; the shared distance still cancels exactly.
  The displayed physical-mass bounds likewise contain only the available fit
  and distance components, not a complete uncertainty total.  The later
  numerical IC 2560 analysis in arXiv:2409.00772v1 is retained only as audit
  history: v2 explicitly removed it because its input data were not officially
  published, and the peer-reviewed 2025 version contains no IC 2560 fit.  It is
  therefore neither adopted nor entered in the conditional-fit sidecar.

### Continuum fitting and 230 GHz anchors

- Added strict, portable resolution of `selected_clean_iterations.csv`.
- Made manifest use opt-in and enabled it only in the clean-image runner; the
  dirty runner cannot be redirected to a clean image.
- Missing selected images now stop the clean run instead of silently falling
  back to `nitercal_0`; all manifest-required iterations are documented.
- Added explicit conversion and validation of CASA component-flux and
  image-brightness units for table values, upper limits, and plotted images.
- Replaced the sub-1-mJy formatter that deleted zero digits with ordinary
  two-decimal formatting, routed fit coordinates from each image's celestial
  reference frame, and removed the non-independent CASA formal-fit term from
  the source-excluded background-RMS S/N denominator.
- Recorded the earlier affected output in the run-specific preapplication
  correction snapshot. That historical 66-row snapshot contains 50 changed flux fields and 54
  changed S/N fields, with overlaps; no detection/upper-limit classification
  changes occurred within that earlier fit-rerun scope.  The separate exact
  reviewed-delta and addition ledgers supersede it as the complete scaled-MAD
  change inventory.
- At the historical RMS-era checkpoint, accepted the regenerated NGC 2960
  85.88 GHz upper limit of `0.47*` mJy in
  place of `0.45*` mJy.  The immutable ICRS FITS input and fitted model are
  unchanged; coherent conversion of the J2000 NED position shifts the
  3-arcsec fitting square by 0.0291 arcsec, changes its sampled pixel mask from
  182 to 169 pixels, and changes three times the RMS from 0.452726 to
  0.467196 mJy.  `ngc2960_frame_correction_audit.json` preserves the input
  hash, coordinates, masks, arithmetic, decision, and completed downstream
  propagation record.  That old classification is now retained only as frame
  and RMS audit evidence: the final scaled-MAD table classifies the row as a
  detection, while its selected 225.53 GHz SMA anchor remains
  unchanged.
- Corrected the NGC 4945 physical anchor from Emig et al. (2020) source 17,
  an SSC/SNR-associated component, to source 18, the AGN/maser nucleus at ICRS
  `13:05:27.4798 -49:28:05.404`.  `rerun_ngc4945_nuclear.py` remeasures all
  27 NGC 4945 images with the centroid and restoring beam fixed.  The historical
  lineage used the source-excluded background RMS once; the final authoritative
  workflow applies the same source-excluded local scaled-MAD policy as other
  ordinary final fits.  Exact A8 adjudication retains all 27 as detections.
- Preserved the rejected 349.40 GHz source-17 result, the accepted source-18
  result, a free-centroid check, and 0.18/0.24 arcsec region/model-sensitivity
  checks in `ngc4945_component_audit.csv`.  These diagnostics are caveats, not
  additional independent formal-error terms.
- Added telescope, date, beam, portable reference, flux status, and stable
  observation ID to `S230_extrapolation_comparison.csv`.
- Added `s230_anchor_identity_audit.csv`, which assigns all 19 selected S230
  observations an explicit physical-component disposition and links NGC 4945
  to both the accepted source-18 and rejected source-17 audit records.
- Added the deterministic 19-row `s230_anchor_astrometry_audit.csv` and its
  builder, policy, WCS snapshot, schema, and tests.  The audit is locked to all
  current S230 statuses, observation IDs, product references, and hashes; it
  records fitted/fixed and nuclear/proxy coordinates in a common ICRS frame,
  signed offsets, geometric-beam separations, and explicit uncertainty,
  phase-calibrator, and astrometric-QA availability.  At the historical
  pre-augmentation checkpoint, its 11 detections had evidence-backed nuclear,
  ambiguous, or blended dispositions and its eight upper limits had no
  invented component.  The final augmented audit contains 12 detections and
  seven upper limits.  NGC 4945 retains source18 as canonical and source17
  as rejected evidence.
- Updated the NGC 5765b dynamical-center provenance to combine the Gao et al.
  (2016) phase-referenced absolute position with the fitted center offsets of
  the adopted Pesce et al. (2020b) model.  The fitted-center update changes the
  native declination by only 0.0005 mas.  Because Gao's phase calibrator comes
  from the VLBA Calibrator Survey, its coordinate is now correctly classified
  as ICRS rather than FK5 (Beasley et al. 2002; DOI `10.1086/339806`); removing the spurious frame-bias rotation changes
  the descriptive separation from 0.0620 to 0.0558 arcsec, or from 0.200 to
  0.180 geometric beams, without changing the identity disposition.  Liang et
  al. (2024; DOI
  `10.1093/mnras/stad3675`) independently identify the selected component as
  the compact nuclear source in the same ALMA program; their integrated flux
  is not substituted for the canonical beam-matched measurement.
- At the historical pre-augmentation checkpoint, corrected selected central
  `S_230` values for NGC 1320 (`0.780445 mJy`),
  NGC 3393 (`0.309386 mJy`), and NGC 5765B (`0.394812 mJy`).  Selected
  flux/SNR-derived measurement terms decrease to `0.018442`, `0.042617`, `2.030285`, and
  `1.814934 mJy` for Circinus, NGC 1068, NGC 3079, and NGC 4258,
  respectively.  The NGC 4945 source-18 correction changes its 349.40 GHz
  input from `15.84` to `13.08 mJy`, its `S_230` from `21.226174` to
  `17.527674 mJy`, and its approximate `S_230` measurement term from
  `0.347662` to `0.350673 mJy`; the selected observation and status were
  unchanged in that checkpoint.  The augmented table selects the newly admitted
  ESO 558-G009 313.11 GHz X600 product; its adjudicated fit is `0.66 mJy`
  in a `0.71817 arcsec` beam and its propagated S230 flux density is
  `0.8190724285773601 +/- 0.09125138464542781 mJy`.
- Documented that `flux/SNR`-derived measurement terms exclude calibration and
  spectral-index/model systematics and are not Gaussian-fit covariance errors.
- Identified the two IC 2560 BB091 image products as 2007 and 2009 NVAS
  reductions of the same 1998-01-06 observation.  The later `BB0091` reduction
  is the sole retained measurement; `continuum_product_alias_policy.csv` and
  the generated `continuum_product_aliases.csv` preserve both filenames,
  hashes, matching FITS-header evidence, fit values, and the exclusion reason.
- Removed the legacy special handling of three failed NGC 4258 NVAS images.
  The 1997-01-04 and 1996-12-31 free fits lie 1.679 arcsec (10.4 beam major
  axes) and 0.553 arcsec (6.5 beam major axes) from the maser nucleus; the
  latter image is strongly artifact dominated, while the former has only a
  marginal nuclear feature and substantially worse noise than the published
  reprocessing. The 2004-05-22 AH0847 image has bright diagonal stripe
  artifacts across the nucleus, so its legacy 79.89 mJy row is also
  unreliable. The active fitter now loads the hash-locked exclusion policy,
  verifies each artifact, and skips all three products.
- Added `continuum_product_exclusion_policy.csv` and the generated
  `continuum_product_exclusions.csv`, preserving all three filenames, immutable
  SHA-256 hashes, FITS-header evidence, legacy table rows, centroid/image
  diagnostics, decisions, and replacement dispositions. The final continuum
  table uses the same-observation Doi et al. (2013) reprocessed detections for
  the first two products:
  `3.5 +/- 0.2 mJy` at 15.365 GHz on 1997-01-04 and `3.1 +/- 0.4 mJy` at
  22.460 GHz on 1996-12-31.  Their scalar beams are the geometric means of the
  published Table 1 major/minor beams; neither value is attributed to NVAS.
  The AH0847 artifact has no replacement row.
- Added `continuum_source_roles.csv` and build-time inventory checks: exactly
  21 continuum sources have role `sample`, while NGC 5495 is retained with role
  `ancillary`.  The 230 GHz anchor and variability builders consume the same
  manifest so the ancillary rows cannot enter sample-level calculations.

### Historical pre-augmentation portable continuum-product recovery

- Added `build_continuum_recovery_manifest.py` and, at that historical
  checkpoint, a deterministic 418-row
  `machinetables/continuum_recovery_manifest.csv`, with an order-independent
  record ID and exact product-or-citation mapping for every canonical continuum
  measurement.
- Recorded archive-specific recovery keys and authoritative retrieval services:
  133 ALMA FITS rows have Member OUS UIDs, 257 NVAS image rows have exact
  program tokens and filenames, and 15 SMA CASA-image rows have project
  `2022B-H002`, date 2023-03-16, track, and selected image path.  The early
  ALMA execution-block-named product is explicitly mapped to project
  `2012.1.00912.S`, Member OUS `uid://A002/X5a9a13/X528`, and execution block
  `uid://A002/X5b1929/X57d`.
- Locked audit-time byte sizes and SHA-256 values for all 405 locally fitted
  products in `continuum_product_checksums.csv`.  Regular images use ordinary
  file SHA-256; CASA image directories use `sha256-relpath-tree-v1`, which
  hashes sorted relative names, entry types, link targets, file sizes, and file
  contents while excluding absolute paths and filesystem metadata.
- Marked the remaining 13 measurements explicitly as literature-only and
  linked each exact table citation to its DOI; no local-file checksum is
  invented for those rows.
- Added a lightweight `--check` that needs no large data and an optional
  `--archive-root` verification that resolves and recomputes all 405 product
  hashes.  Neither the checksum lock, generated manifest, nor summary contains
  an author-local path.

### Historical pre-augmentation variability matched-pair screen

- Replaced notebook table writes with `build_variability_tables.py`, which has
  no NED, distance, or black-hole dependency.
- Regenerated all 10,597 within-source pairs and the 98-row historical
  frequency/beam matched-pair screen from the 415 measurements belonging to
  the 21-source sample. All 418 canonical continuum records remain in the
  continuum product, including the three additional NGC 5495 catalog rows.
- The corrected numerical fields change 1,088 all-pairs records and six
  matched-screen records, while physical pair membership and both row counts
  remain unchanged.
- The later NGC 4945 fixed-nucleus rerun changes all 351 within-source pairs
  containing its 27 measurements and six matched-screen rows; physical pair
  membership among the retained rows remained fixed, and the selected screen
  remained at 98 rows at that historical checkpoint.
- Restored the omitted NGC 1320 pair and the true second-frequency values.
- Added stable record/pair IDs, complete input provenance, a header/schema,
  absolute symmetric contrasts, correct independent formal-error propagation,
  explicit upper-limit censoring, range-beam handling, and deterministic dates.
- Marked the two remaining same-day rows as non-temporal.  The spurious
  IC 2560 self-comparison is absent because the older reduction is excluded
  upstream and preserved in the explicit alias sidecar rather than silently
  deleted.
- The four geometry-matched censored pairs involving the two failed NGC 4258
  upper-limit rows disappear. The Doi replacements do not enter the strict
  matched screen because their published geometric-mean beams fail the
  2.6 percent beam criterion. Removing the artifact-dominated AH0847 row
  removes its 86 all-pairs comparisons but no selected-screen row, leaving
  10,597 all-pairs and 98 selected rows.
- Marked the old cleaner report and notebooks as legacy and prevented them from
  overwriting the corrected machine products.

### Contamination products

- Replaced manual notebook/table edits with `build_contamination_tables.py` and
  explicit literature-input and adoption-policy CSVs.
- Fixed the target at exactly 230 GHz, with strict dust `nu > 253 GHz` and jet
  `nu < 207 GHz` eligibility; reported measurement-scale terms use the same projection factor and
  input upper limits remain censored.
- At the historical pre-augmentation checkpoint, generated complete candidate, selected, continuum-anchor, JSON, LaTeX, and
  20--2000 K modified-blackbody products with physical-line or literature IDs
  and portable provenance. That table had 281 candidates and 14 selected policy
  results.  The final augmented build contains 283 candidates and 14 selected
  policy results.
- At the historical pre-augmentation checkpoint, corrected NGC 4258 dust to
  `2.0743 +/- 0.6486 mJy`, NGC 4258's VLA date, NGC 1194 jet to
  `<0.09337 mJy` without an invented sigma, NGC 3079 jet to
  `0.9050 +/- 0.1876 mJy`, and Circinus to the ATCA provenance.  The four
  Elmouttie et al. (1998) Circinus inputs are now identified correctly as
  central peak brightnesses rather than integrated flux densities.  The
  selected `7.11 +/- 0.71 mJy` projection is numerically unchanged, and its
  machine-readable input preserves the source table's blank integrated-flux
  field and ATCA provenance.
  The final scaled-MAD table classifies the NGC 1194 1.49 GHz row as a
  detection.  Its regenerated conservative attribution bound is
  `0.056341350943041435 +/- 0.007673842405753396 mJy`, or
  `3.364822225211587 percent` of S230.
- Resolved production query Q8 by replacing the ambiguous NGC 4258 dust
  asterisk with an explicit upper-bound symbol (`<=`) and requiring the
  dust-model caveat to be stated in words in the table note.
- Consistent 230-relative selection gives NGC 1068 dust `0.0002170 mJy` and
  NGC 4945 dust `0.0782213 mJy`; the rejected NGC 4945 jet remains diagnostic
  only and is not emitted as a table row.
- Resolved the NGC 1068 jet provenance directly from Gallimore et al. (2024):
  Section 2.3 states that the 21.9 GHz S1 continuum model contains six Gaussian
  components and two point sources and reports their total HSA flux as
  `7 +/- 1 mJy` on printed page 9.  The summed-components result is now adopted
  and produces the canonical table; the historical image-peak result remains a
  clearly non-adopted diagnostic.

### Marked scientific impact and cause

- The obsolete perfect-correlation and independent-quadrature prescriptions
  are replaced by paired posterior ratios or exact shared-distance
  cancellation.  The resulting shadow-size values and uncertainty components
  change manuscript Table 1 and the BHS-versus-flux figure.
- The corrected BHS calculation retains NGC 4258 as the only sample source
  angularly resolvable on an Earth--L2 baseline, but angular resolution and an
  arcsecond-scale total flux do not by themselves establish shadow-scale
  correlated-flux detectability.
- At the historical pre-augmentation checkpoint, three faint-source `S_230`
  central values changed by `+1.32%` (NGC 1320),
  `-11.4%` (NGC 3393), and `+2.56%` (NGC 5765B).  These changes do not alter
  the four sources above 10 mJy, the brightness ranking relevant to the main
  space-VLBI argument, or any selected 230 GHz observation at that checkpoint.
  The final augmented architecture changes the ESO 558-G009 selected
  observation to the admitted 313.11 GHz X600 product: `0.66 mJy` in a
  `0.71817 arcsec` beam, giving
  `S_230 = 0.8190724285773601 +/- 0.09125138464542781 mJy`.
- Correcting the NGC 4945 physical component lowers its selected `S_230` by
  17.4 percent, from `21.2262` to `17.5277 mJy`.  It remains the
  second-brightest selected source, the same four sources remain above
  10 mJy, and the paper's qualitative brightness ranking is unchanged; the
  continuum table, contamination table, spectral/variability diagnostics,
  and BHS-versus-flux figure nevertheless require the corrected numbers.
- The submitted NGC 1068 table used the 21.9 GHz image peak as though it were a
  single integrated component and gave `0.140 mJy`, whereas the submitted
  method said to sum all resolved components.  The source-reported total HSA
  input projects to `2.16 +/- 0.31 mJy` (about 24 percent of `S_230`), so the
  revised paper now describes the constraint as subdominant and retains it.
- NGC 4258's simple dust bound changes from `0.87 +/- 0.32 mJy` to
  `2.07 +/- 0.65 mJy`, and its modified-blackbody maximum changes from 13
  percent to 29.9 percent.  The conclusion that dust does not dominate remains,
  but the submitted percentages and any stronger wording do not.
- The two NGC 4258 centimeter-wave table arrows become correctly cited
  detections.  This changes the NGC 4258 SED and continuum/variability
  diagnostics but not its selected 225.51 GHz anchor, `S_230`, contamination
  adoption, black-hole size, source-above-10-mJy set, or headline feasibility
  conclusion.
- Removing the artifact-dominated AH0847 79.89 mJy row leaves 86 NGC 4258
  continuum measurements and removes one non-selected contamination candidate;
  the 14 selected contamination results and NGC 4258 adoption are unchanged.
- This happened because the Gallimore et al. (2024) notebook input contained
  only one row labeled `Center`; the generic summing function therefore summed
  one component while the independently edited manuscript described a full
  component sum.  Source-dependent target frequencies, hard-coded later cells,
  stale notebook state, and a faulty metadata formatter explain the other
  contamination discrepancies.  The new deterministic generator and tests
  replace that fragmented notebook-to-LaTeX path.

### Workflow safeguards

- Replaced the missing `get_chosenfreq.ipynb` handoff with
  `chosen_frequency_anchors.py`.  The replacement records frequency, beam,
  WCS frame, fitted coordinates, portable product identity, input SHA-256,
  CASA-model tree SHA-256, run context, and the reviewed component
  disposition; exact ties are deterministic.
- Made component identity a fail-closed part of chosen-frequency
  consolidation.  A generic minimum-beam winner must match the reviewed
  product, frame, coordinate tolerance, and handling in
  `chosen_frequency_component_policy.csv`, derived from
  `s230_anchor_identity_audit.csv`.  A synthetic test with a fainter reviewed
  nucleus and a brighter off-nuclear free-fit result in the smallest-beam
  image confirms that the off-nuclear component cannot silently win even when
  the product identity itself matches.
- Preserved the nine historical ALMA anchor identities in
  `chosen_frequency_anchor_reference.csv`, but explicitly marked the NGC 4945
  source-17 SSC/SNR-associated seed as rejected.  Its input and model-tree
  hashes were verified against all nine preserved hardware artifacts.  The
  existing isolated final-stage evidence remains 409 accepted rows and 46
  exclusions with exact scientific-field agreement.  Added
  `historical_final_stage_validation.json` to preserve its counts, CASA/code
  provenance, artifact hashes, source-evidence hashes, environment limitation,
  and explicit historical-only scope without bundling the generated model
  trees.
- Added `run_continuum_four_pass.py`, whose deterministic `prepare` command
  originally built and verified the historical 455-input isolated layout and
  now extends that baseline by the exact addition policy.  For final release it
  must verify 462 inputs and copy 447 ALMA/VLA FITS inputs with source/copy
  SHA-256 equality, together with all SMA CASA products, so
  the preserved archive cannot be mutated.  The runner refuses Git worktrees,
  directory-symlinked CASA inputs,
  changed input inventories, and pre-existing outputs; gives each legacy stage
  closed standard input; runs the four historical modes in order with
  notifications off; and preserves per-stage logs and hashes.  It then
  requires the guarded NGC 4945 source-18 fixed-nucleus rerun/check and an
  exact policy-derived `406 + P` comparison against the immutable historical
  409-row baseline after the three reviewed NGC 4258 exclusions.  With the
  seven final publication-row additions, `P = 7` and the generated table must
  contain 413 rows. The uncorrected generic
  fixed-to-chosen output is explicitly intermediate.
- Completed the historical pre-scaled-MAD corrected run on 2026-08-18 (HST)
  under run ID
  `current-corrected-230316-flag2-retry3`.  All four generic passes, the
  nine-anchor check, the mandatory NGC 4945 source-18 apply/check, and the
  post-run inventory guard passed.  The runner produced 406 rows that exactly
  match the frozen 409-row historical raw reference after the three policy exclusions,
  with no missing, extra, duplicate, or field-mismatched records.  All 406
  product paths are release-relative, including all 15 SMA paths.
- Refactored the fit search geometry through the tested
  `square_bmaj_region_width_arcsec()` helper.  The region is intentionally a
  square whose two full side widths come from `Bmaj`; the unused `Bmin` width
  variable was removed so the implementation no longer resembles an ellipse.
- Updated `post_processing_finalcontinuum.ipynb` to call the deterministic
  builders and removed every automatic legacy variability, contamination, and
  size-derived notebook execution.
- Added an explicit `rerun_ngc4945_nuclear.py --apply-recorded` correction
  stage after combining raw fitter summaries, so post-processing cannot
  silently restore the rejected source-17 measurements.  A live CASA rerun
  remains the authoritative way to regenerate the 27-row sidecar from the
  public FITS inputs, while `--check` verifies its integration without CASA.
- Strengthened that check to require all seven accepted/rejected/sensitivity
  component records, the immutable 27-product input SHA-256 mapping, and
  consistent unrounded flux, background RMS, rounded flux, S/N, and fit-box
  arithmetic.  Deleted diagnostics or tampered hashes/numerical fields now
  fail before `--apply-recorded` can alter the sorted table.
- Made all 15 SMA provenance entries in the sorted and final fit summaries
  repository-relative under `analysis/sma_calibration/`.  Both summary
  builders normalize recognized legacy paths and reject unknown absolute paths;
  the downstream machine table, S230, variability, contamination, and BHS
  numerical products are unchanged.
- Routed black-hole export through `build_bhs_params.py --require-final`, so a
  release export cannot silently bypass any unresolved author-response gate.
- Preserved the original hardware checkout and the v1.0.8 baseline; all
  black-hole corrections are isolated on the correction branch.

## v1.0.9 finalization record

- Reserved Zenodo version DOI `10.5281/zenodo.22051996` and inserted it into
  the manuscript bibliography, response documents, release description, and
  citation metadata.
- Updated `CITATION.cff` to version v1.0.9, release date 2026-08-26, and the
  exact reserved version DOI.  The release retains the MIT license and now
  includes the complete license text in `LICENSE`.
- Corrected the software attribution metadata: Roman N. Burridge is the sole
  software author, Creator, and copyright holder; Geoffrey C. Bower remains a
  coauthor of the associated paper and appears in Zenodo only as a Contributor
  with role `Supervisor`.
- Regenerated the clean manuscript, official redline/marked PDF, response
  portal text, delivery-package manifests, and ZIP archives from the final
  DOI-bearing source, then reran the regression, TeX, checksum, and unpacked-
  package audits.
- Zenodo publication and ApJ delivery are external actions and are not
  performed by this repository.  Each requires separate explicit author
  authorization and post-action verification.

The v1.0.9 record was subsequently published at
`10.5281/zenodo.22051996`.  That published archive remains unchanged.

## v1.0.10 finalization record (DOI reservation pending)

- Prepared v1.0.10 as a new Zenodo version because the public v1.0.9 archive
  predates the authoritative Task 28 continuum closeout and retains the earlier
  425-row continuum product.
- Packaged the exactly adjudicated Task 28 scientific state: a 501-row final
  continuum table (485 local image measurements and 16 literature records),
  the recovery and exclusion sidecars, and all regenerated downstream
  230 GHz, contamination, table, and figure products.
- Retained the black-hole parameter and uncertainty corrections introduced in
  v1.0.9.  No already-published v1.0.9 file is overwritten.
- Reserved v1.0.10 version DOI `10.5281/zenodo.22121490` and inserted it into
  the citation metadata and release description.  The manuscript bibliography,
  response documents, and final archive must use the same DOI before upload or
  delivery.

Every contamination adoption in `contamination_adoption_policy.csv` is now an
implemented decision rather than a pending choice.  NGC 1068 adopts the
summed-components result required by the manuscript method, with its scalar
input traced to Gallimore et al. (2024), Section 2.3, printed pages 8--9; the
central-only image-peak calculation is retained only as a non-adopted
diagnostic.

The continuum product decisions are also implemented rather than pending:
the later IC 2560 NVAS reduction is the sole measurement with the older alias
preserved in a provenance sidecar, and NGC 5495 remains a catalog-only source
outside the 21-source analysis sample.

## Publication precision and source-fidelity closeout (2026-08-19)

- Retained full numerical precision in machine-readable tables and calculations
  while giving the publication renderers one explicit display policy:
  uncertainty components normally use two significant figures, central values
  follow the finest displayed uncertainty place, and derived values without an
  uncertainty normally use three significant figures.
- Regenerated manuscript Table 1 with deterministic decimal rounding and
  source-faithful inclination entries.  NGC 1068 now uses the nearly edge-on
  Gallimore et al. (2024) geometry, NGC 5765b uses the updated Pesce et al.
  (2020b) $r=0$ warp intercept, and CGCG 074-064 retains both its statistical
  and systematic inclination terms and its separate $\pm0.009$ dex
  marginal-mass model-selection systematic.  That marginal mass term is not
  propagated independently into the paired shadow-size interval.  The
  published asymmetric intervals are retained for UGC 6093 and NGC 6323.
- Added `build_continuum_latex_excerpt.py`, which generates the ten printed
  Table 2 rows from `fitsummary_machine.txt`; this fixed the stale Circinus S/N
  entry and makes half-up display rounding reproducible.
- Regenerated all six contamination tables with the common display policy and
  three-significant-figure non-probabilistic spectral-index bracket endpoints;
  the Circinus table also carries the corrected central-peak ATCA provenance.
- Regenerated all six SED figures with distance-derived angular-to-linear
  insets and quadrature totals of the listed numerical distance components.
  Figure captions now define detections, upper limits, inset uncertainties, and
  the logarithmic geometric-mean beam color scale.
- Corrected the BHS figure axis label and synchronized its caption and prose
  with all four plotted total-flux thresholds (0.1, 1, 10, and 100 mJy per
  beam), the plotted uncertainty semantics, and the adopted baseline markers.
- Routed the plotting cells in all four historical notebooks to
  notebook-specific diagnostic directories.  Regression tests prevent them
  from deleting or overwriting the canonical `SpectraGraphs` publication
  products; other notebook outputs retain legacy, noncanonical names.

## Verification record

All continuum counts and hashes recorded below are dated historical
pre-augmentation checkpoints unless an entry explicitly identifies the final
462-input augmented execution.  They remain valuable audit evidence but cannot
serve as final-release acceptance values.

Command:

```bash
python3 -m unittest discover -s analysis/multi_freq_from_archive/tests -p 'test_*.py' -v
```

Result on 2026-08-12: 43 tests passed; `git diff --check` passed; all three
edited notebooks parsed as valid JSON; deterministic outputs contained no
author-local `/home/` paths.

Additional black-hole integration verification on 2026-08-17:

```bash
python3 -m unittest discover -s analysis/CountingBHS_and_Spectra/tests -p 'test_*.py' -v
python3 -m unittest discover -s analysis/multi_freq_from_archive/tests -p 'test_*.py' -v
python3 analysis/CountingBHS_and_Spectra/build_bhs_params.py --check
python3 analysis/CountingBHS_and_Spectra/build_bhs_params.py --check --require-final
```

- 33 black-hole parameter and manuscript-table tests passed.
- 64 integrated continuum, variability, contamination, BHS utility, figure,
  notebook, and workflow-guard tests passed.
- The ordinary deterministic `--check` passed for the 21-row BHS table.
- At this 2026-08-17 checkpoint, the `--require-final` command exited with
  status 2 as designed and named IC 2560 as the remaining author-decision gate;
  that historical checkpoint was resolved by the formal-only decision below.
- The manuscript Table 1 regression and checked-in PDF/PNG publication figure
  matched their deterministic builders byte for byte.
- `git diff --check` passed after the release-documentation update.

Continuum-correction integration verification on 2026-08-17 regenerated the
continuum, S230, variability, contamination, and repository BHS-figure products
twice with identical SHA-256 hashes.  After the NGC 4945 and anchor-identity
integration, strict sidecar/path guards, SED regeneration, and frame-routing
checks, all 143 integrated multi-frequency/continuum tests and all 33
black-hole tests passed;
the ordinary BHS deterministic check and `git diff --check` also passed.  The
IC 2560 gate was subsequently closed on 2026-08-19 by recording the
formal-only treatment and leaving the model contribution unquantified.

Final black-hole closeout verification on 2026-08-19 regenerated the canonical
21-row table and manuscript Table 1 without changing the IC 2560 central value
or formal interval.  Both `build_bhs_params.py --check --require-final` and the
guarded export exited with status 0; all 37 black-hole tests and all 175
integrated continuum tests passed.

The portable recovery audit was run against the preserved project tree on
2026-08-17 and regenerated for the accepted third exclusion. All 405 fitted
products resolved at their exact release-relative
locations and matched their recorded sizes and hashes; 13 literature-only rows
were separately identified, and the local payload total is 3,929,414,510
bytes. Ten focused manifest/hash tests and the
data-independent deterministic `--check` passed.

The chosen-frequency reference check independently matched all nine preserved
input hashes and all nine CASA-model tree hashes.  Focused tests cover
deterministic tie-breaking, component-policy enforcement, the synthetic
fainter-nucleus/brighter-off-nuclear case, square-`Bmaj` region construction,
workspace isolation, four-mode ordering, and the mandatory NGC 4945 source-18
post-pass guard.

The historical pre-scaled-MAD corrected full CASA execution completed on
2026-08-18 (HST).  Its
six stage-artifact snapshots and logs are preserved under
`end_to_end_validation/current_four_pass_20260818_retry3/release/analysis/multi_freq_from_archive/workflow_runs/current-corrected-230316-flag2-retry3/`.
Independent read-only verification found:

- 406 generated rows, SHA-256
  `da0c6aa568b70cbba0899926538a9e2ae16de26d44ec4b76823f628d5f5f6778`,
  exactly equal to the 409-row raw reference after the three hash-locked
  exclusions;
- an unchanged 455-input prepared inventory (133 ALMA, 307 VLA, and 15 SMA):
  the runner rehashed all 440 ALMA/VLA FITS contents and natively rechecked the
  SMA count and absence of directory symlinks; a separate post-run deterministic
  tree-hash comparison then matched all 15/15 SMA discovery trees and all 6/6
  selected-clean image/residual trees to the original archive;
- nine validated chosen-frequency anchors; all 406 product paths are relative,
  resolve inside the isolated release, and include 15 portable
  `analysis/sma_calibration/...` paths; and
- a passing NGC 4945 source-18 check covering 27 measurements and all seven
  accepted/rejected/sensitivity component records.

The exact-comparison JSON has SHA-256
`cc04f7c735f7ea1203d621a5e25a04a6de9d57bee554a236044cf24472b50e7a`;
the input-integrity JSON has SHA-256
`ea6928f95fe943cb9ca928b69c32a9ccf12e02c04f0a660da98da4911518edd4`.
The aggregate SHA-256 of the sorted 15-tree discovery-name/hash mapping is
`241b550ab23cbd9f5560248c0ae731fc37e81b339fbcbd5ee8b6aa7364c120c1`;
the corresponding six-tree selected-clean mapping SHA-256 is
`4cab219804e4d2b6f22f25a6dd5cbfc6d683e43fdd0f887ae3631bdb8ebfc4a8`.
These two SMA mapping hashes are independent post-run audit evidence; the
runner's prepare/integrity JSON does not yet preserve them automatically.
After this evidence was integrated, all 13 four-pass-orchestrator tests, all 9
chosen-anchor tests, and all 5 NGC 4945 guard tests passed; `git diff --check`
also passed.

### Authoritative A8 current-release fingerprint

The schema-4 A8 run completed and was exactly adjudicated on 2026-08-20 (HST).
Its 462 inputs produced 413 publication rows and 416 generic diagnostic-noise
rows, including three reviewed boundary rejections.  The downstream continuum
table contains 425 rows (412 image-based plus 13 literature), of which 422
belong to the analysis sample.  Those sample rows produce 10,787 unordered
within-source pairs and a 98-row strict screen (96 different-day plus two
same-day).  The 19 selected S230 rows comprise 12 detections and seven upper
limits from ten ALMA and nine SMA observations.  The recovery manifest records
140 ALMA, 257 NVAS, 15 SMA, and 13 literature rows and 3,943,670,510 local
payload bytes.

The final NGC 1194 1.49 GHz measurement is a `0.70 mJy` detection with local
scaled-MAD-derived `sigma = 0.09534186870062653 mJy`; its conservative projected
jet attribution bound is `0.056341350943041435 +/- 0.007673842405753396 mJy`, or
`3.364822225211587 percent` of S230.  NGC 4945 has
`S_230 = 17.527674085288034 +/- 2.485842303969371 mJy`; the uncertainty is the
propagated scaled-MAD measurement contribution and excludes calibration,
spectral-index/model, structured-background, aperture-choice, and variability
terms.

Authoritative A8 and promoted review artifacts:

- workflow plan: `00b0b3dff3c5cb00920934eaef87f503a4c99172e0fd7e2074d1eb0eff40c9eb`
- 413-row sorted fitter table: `2756c1b5abc1d548dea76fba826602a8ddbe04b18424aded4cdf7c68bdc839c9`
- 416-row final noise audit: `f4c4a1e2e3170e6b3a2db222305bfdd86f70a2166ddd411c5141849b57f6c93c`
- reviewed delta ledger: `0da4af0439d549c4a62759d20d2f80291e5eafd15c20288801406a97c810a2da`
- reviewed addition ledger: `479a34f7bfb9d61d2a835d96542f8409d44011158f9aa2ed32a4d96432b668fe`
- schema-3 review freeze: `5aa6b7c2fdefaf20b9fa751424a58ac1c5d8fe807b31c7f75fe3c56f0670322d`
- schema-2 adjudication: `f29c75068f4c1329156586f4bf6a68daa6e2bc331709a79e8e11464c62eb6234`

Regenerated downstream scientific artifacts:

- final 425-row fit summary: `b9c13bc6db28419bdb785595ef9f0729bfa0ef74873fdb8dede9842d2da4c1e5`
- `fitsummary_machine.txt`: `2dcf99465685cb807d451da607acd70a3a6e62760764a47ee03a0a71c4032965`
- `fitsummary_aas_mrt.txt`: `831a9e9f458cab5726281249ae7aef070639ebe3cea538ec7a8309296cf95fd6`
- recovery manifest: `8b5712f8a921b25d73218a74d60f8ada83929d89bcdb94267034ea44396267a4`
- S230 table: `12171f167d75a70ce30893178cebe8e7ba33e9a2ef54f1785d33b211c3e18db7`
- S230 identity and astrometry audits: `537dd92fe5254ba679b505a3cae0a3c608dd88b6fb5457759c5e4d6487ef6e73`, `0c15c4b2bd53cd318a9ba7a1e269cccb50aaa81439fa14a01827bc3916b08e95`
- all-pairs and strict variability tables: `494ba5f5f8e34c06e90ed66e4920b31b66965a26ec6c81faa1fe005fa4829832`, `2e6ada606755ea356eebc83038cac2a94467227dd20083ab0b743b1944fe2f37`
- contamination candidates, selections, and result audit: `89f1d7c7b9933e1e9d60686293965fc7d537ffe96a22566c6bcc047260a3a33b`, `cd760e92651b232ff23ba7e261ec416467cb2d0f8cd947c6b4aaef9a97c88286`, `a3008464d2aa844f2918091b61fb60a99007f3f8ab445ddba5dbdff1299ecb70`
- A8-packaged `BHS_params.csv` (superseded by the subsequent source-by-source
  black-hole audit; not a final packaging hash):
  `76427aff0c2018b52e7c3d1a5a1c6475722383c4811eb171cf72aabac8606fa4`
- BHS PDF/PNG: `ae38e406d5a7eecd92a01e2dfb3084fe4e95893fcb30eef4d3324ea3ed7a5825`, `52b57dd8a2ba4e41f574ba3ed5e3635c81a74782c954406a8b19c7b3dd72093e`
- six-SED manifest: `d7187a92158e1ed6f708714febcd04440e11255c1071a2cae7f0fde4a0cf284d`

Historical pre-augmentation regenerated file line counts include a header where
the format has one:

- `fitsummary_machine.txt`: 418 lines (418 canonical continuum measurements;
  the AAS version adds metadata header lines).
- `continuum_product_aliases.csv`: 3 lines (2 products in one alias group).
- `continuum_product_exclusions.csv`: 4 lines (3 failed NGC 4258 products, 2
  exact literature replacements, and 1 explicit no-replacement disposition).
- `continuum_product_checksums.csv`: 406 lines (405 fitted-product hashes).
- `continuum_recovery_manifest.csv`: 419 lines (418 canonical continuum rows).
- `continuum_recovery_manifest_summary.json`: 23 lines (coverage, type counts,
  input hashes, and directory-checksum definition).
- `continuum_source_roles.csv`: 23 lines (21 sample sources plus one ancillary
  source).
- `validated_continuum_corrections.csv`: 67 lines (66 corrected records).
- `ngc2960_frame_correction_audit.json`: 112 lines (one hash-locked frame/RMS
  correction decision and its downstream-regeneration gate).
- `ngc4945_nuclear_photometry.csv`: 28 lines (27 fixed-nucleus measurements).
- `ngc4945_component_audit.csv`: 8 lines (one rejected component, one
  canonical component, and five validation/sensitivity records).
- `s230_anchor_identity_audit.csv`: 20 lines (19 selected S230 observations).
- `s230_anchor_astrometry_audit.csv`: 20 lines (19 selected S230 observations;
  11 detections and 8 upper limits).
- `chosen_frequency_component_policy.csv`: 10 lines (9 reviewed generic-anchor component decisions).
- `chosen_frequency_anchor_reference.csv`: 10 lines (9 preserved historical anchor identities).
- `historical_final_stage_validation.json`: 42 lines (compact historical validation evidence).
- `S230_extrapolation_comparison.csv`: 20 lines (19 sources).
- `variability_machine.txt`: 99 lines (98 screening rows).
- `correlations_machine.txt`: 10,598 lines (10,597 pairs).
- `variability_likely_duplicates.txt`: 1 line (header only; zero surviving
  likely-duplicate groups).
- `contamination_candidates.csv`: 282 lines (281 candidates).
- `contamination_selected.csv`: 15 lines (14 policy results).
- `continuum_anchors.csv`: 7 lines (6 manuscript contamination sources).
- `modified_blackbody_20_2000K.csv`: 5 lines (4 dust sources).

The six SED angular-to-linear insets are now derived directly from the adopted
distances in `BHS_params.csv` rather than from a manually copied scale map.  A
single displayed inset uncertainty uses the quadrature total of the numerical
distance components.  This corrects NGC 1068 from the obsolete
`0.068 +/- 0.010 kpc` scale to `0.0520 +/- 0.0025 kpc` and makes the NGC 4258
and NGC 1194 inset uncertainties complete (`0.00054` and `0.019 kpc`,
respectively).

Superseded interim closeout SHA-256 records (2026-08-19; retained only as
historical audit evidence):

> **Do not use this block for packaging.** Subsequent source-by-source science
> and metadata audits changed the files below.  The authoritative SHA-256
> record will be regenerated only after the scientific freeze and all required
> deterministic rebuilds, as required by the upload checklist.

- `bhs_source_inputs.csv`: `d2a7a95706238c0530a44db20f545068c78efd5ba26646bfb57ee2406431f3f7`
- `pesce_theta_g_quantiles.csv`: `3a4088b8dc7c76c2cab3814cd4d61076f860ef476c3544360ceebe814290f668`
- `bhs_latex_metadata.csv`: `b7bd49bafc6882afa8898019ffd7b6407d7ea738d2c8eb137daf466b9230471a`
- `BHS_params.csv`: `7805aa348ca7a9c14d5726950fa5ca592966b598eeae432838a2abdadc9a8af2`
- `build_bhs_latex_table.py`: `b96df6bad18a45f1004a97dfce7f1f5d596c61c0bfedabfa9d2a567290471e59`
- `tests/test_build_bhs_latex_table.py`: `0a92076fed3c6661c50657a4f474ce70dd5f9c6f55e7ac6a3effcb9d31e99fcb`
- adjacent manuscript `BHS_Params.tex`: `07e7a4b7e1cd08f9d073509c8f6c0467e2d604c8289b7bd4aff0ec9de632365f`
- `build_bhs_flux_figure.py`: `54510cbea73b1572fc900e2bec3c2ff35f49878ff8cd2616dc64aba290c7d3c3`
- `tests/test_build_bhs_flux_figure.py`: `b90a99a069e3595c58f8f7ca353912f08cc9a216e519629b785bb935d0c23769`
- `BHS_Flux_Graph.pdf`: `29166efd897b61dc171683811ea1822ad92443d183be3bcd76088ab7ab20e5ce`
- `BHS_Flux_Graph.png`: `60c72f9ae4e0d258134abab4a88d521498a30f75bef71e3915772c465882554a`
- `build_continuum_latex_excerpt.py`: `4e4bf286924db5b33eb249426471921045e972b15461116074d7c3bfa2e53f82`
- `tests/test_build_continuum_latex_excerpt.py`: `6f7403ec8d2a1816525de6d0f84e468675d4fd5bba78ed8df9fbf41da545ba31`
- adjacent manuscript `Continuum_from_Pipeline.tex`: `315391288acb47ed750bc6c38d085c8c21cc158de7f2c9c70c8c8c3c1aec2a1f`
- `Beam_Date_Flux_Graphs.ipynb`: `417f71460747891f7041896c9301fc5049353483d46f282437eb2089f076489f`
- `Variability_Search.ipynb`: `74928274c439bb8f4191d0063f1c7d46dd9f898f51729327786c15994f364402`
- `flux_size_graph.ipynb`: `e2b460838da3f1b846bbf50bde99401b11900229011a034e66fce32e16ea0be1`
- `get_contams.ipynb`: `c152d36137bf395cdce8860dec93873f3aa92f97e34b83cda42ec41325722a83`

- `continuum_product_alias_policy.csv`: `bf1172d33090c77c9a4d96879080dc6a4f80626d02eeebfd49e7226b12a20084`
- `continuum_product_exclusion_policy.csv`: `6c8e51c45ee40752a951739f51ad377e02ff9e535af1f60cc69da77853d3c334`
- `continuum_source_roles.csv`: `713b02cf4c04683a6a731151dac0e8214f672aa62ceeb1bec598e02153be846f`
- `validated_continuum_corrections.csv`: `48e0bc4280d9ee23f0e3f8578649aaededd5bf14109aa6f3f19c0377d5bce448`
- `ngc2960_frame_correction_audit.json`: `bef4c577c5098270f9c54bd9bdce4ec58bf96b21b479180c817ba7aaf80f2084`
- `ngc4945_nuclear_photometry.csv`: `1afe285e91679a712c10b4a40f8cc8d2ed15694b1c1b228dca353d9f703d70b8`
- `ngc4945_component_audit.csv`: `cac748518a9b0f8f57fbea2e9bb10b378354f507ee4058547dfb1803064e7594`
- `s230_anchor_identity_audit.csv`: `39674065e1f02f516d9c9ddf8dc13b12dacb0eb788d33b058c1d66f6b9cf0e16`
- `s230_anchor_astrometry_audit.csv`: `c1c1fc989c83fe97e44a2e377431245fcfc9e876282f88c8cc1f9f10278e23b4`
- `s230_anchor_astrometry_schema.md`: `cf1df9f3400a63579104c4f440bb3748be941b935e589b5dac92ee2ac966fda4`
- `s230_astrometry_policy.csv`: `3f8da01946af07b6a6695dd234a9ea17cc5ec5214cd271b38ec7adf35abf6b26`
- `s230_wcs_evidence.csv`: `87f5d8002fe621c593eb96d5b5dcdbd8bcaec40f91952e949615907faa0d4a25`
- `extract_s230_wcs_evidence.py`: `6999e81626b78c88a84f0ffae446e1124431d8c7cc5e24db622c5be605fb6974`
- `build_s230_astrometry_audit.py`: `0c60f7f8d0e0cafd4ded2f1eb715d6ec7df0517725620f851ec61d0c91ed501d`
- `run_continuum_four_pass.py`: `0b05d9d4718c8a6db809427c763942ef000838d5f568440e023d2f22996706b6`
- `CONTINUUM_FOUR_PASS_WORKFLOW.md`: `7854dd235f31843925932095a9f99f3e3c3f7bd7a3f71312a7d2eaaf53b86f44`
- `chosen_frequency_component_policy.csv`: `c432dbd67dcd3ed419d4841bf802c2e4cd463b9d2696489cee1f655d5fc6f244`
- `chosen_frequency_anchor_reference.csv`: `bc05af0d56942813e5c20cb46c0eeff5433a20256aa52d0568b75c3f69169fb0`
- `historical_final_stage_validation.json`: `9d2acbf1ab136ca8192825276f2e793795c0cfcca6e521a372efc22a85327abe`
- `fitsummary_fixtochosen_sorted.txt`: `c1b3f8e5f3829fab4eecc8bdf8d0e537ee77d68ab89f711189f897f2f2bdf4f8`
- `fitsummary_final_withextras.txt`: `ba48dfb4e2502e179258f7dc364f49f433e0097c6c69a03040c6720745cc07a8`
- `fitsummary_machine.txt`: `1bd64fecc4b99d6f48f890693ca8c385a9e2f2d0a7fab6be2dc2c25654a3b330`
- `fitsummary_aas_mrt.txt`: `3559ebaa7b1c4e58f49ef417e676e20eb412ab1f0f8cf396880e9b604cb9eeea`
- `AAS72796R2_datafile2_corrected.txt`: `3559ebaa7b1c4e58f49ef417e676e20eb412ab1f0f8cf396880e9b604cb9eeea`
- `continuum_product_aliases.csv`: `99b445d9d043c36446663d0bbd8704ddab6ee35d7a4cf6f1135ec3941896bc5c`
- `continuum_product_exclusions.csv`: `3266ce5d90cd2db6f5b3e6695b3c3f077aa25760bf6ca9af821dce3b6d4715d1`
- `build_continuum_recovery_manifest.py`: `f5eddcd82083dca45c07fb59ad52ed74e80b55460197c6e7baab8a757eedcdcb`
- `continuum_product_checksums.csv`: `d1189aee3e7b0521d894ebc5f68c7c4bf6c484857fe74805224c5137606acded`
- `continuum_recovery_manifest.csv`: `e18b6ff5ee6439267cc9303bd6d5835ba3ed9c7299878300de03ae4062892b21`
- `continuum_recovery_manifest_summary.json`: `0d30445ff36d9e3a19824144e42d8b3dee5ffc704b524bce75bd078e346cc7e7`
- `S230_extrapolation_comparison.csv`: `8a82fd37129babf4c4d04ef26d5d417f680c74207d9a708c65f3ee688e33e186`
- `variability_machine.txt`: `7bf0ed29f9c1ad1b39e85a1569901f743a8e00ea870a9cf0ea0d8d34d0400861`
- `correlations_machine.txt`: `a9a0ac1b6f36f4a015224bb05242a6906aeb853e6a7da6147d978b8a4cb141b1`
- `variability_likely_duplicates.txt`: `495aca755b81520b7bbb2a5031753b4271c43ae587859069d4cf31e62d1edd21`
- `variability_tables_schema.txt`: `3ace665f95e9262274d90506d378a250f7660f2cfdca8772a14ff99967945fd6`
- `contamination_literature_inputs.csv`: `feeb712db9168f5e358817e2cf2524cd7f8cb0396c219ff3b4f07f2473d8ebf5`
- `contamination_adoption_policy.csv`: `5c4af72c2b8bdf56bf67837c2a42eb477f2852159af2233cb75b84e177371f5b`
- `contamination_candidates.csv`: `69ac4164be05dbf81251da2fa2459a15fb7dfedf30da77fc10dbf1de361577a5`
- `contamination_selected.csv`: `4bcb49acb4507ef7c153d4f6cd008e072ea13ec2519a1c69447c321cfae2be61`
- `contamination_results.json`: `f23277ba8ae36dcb6bbc75ffb49566ffffb101711528da1a0f5ba64fafe49df0`
- `continuum_anchors.csv`: `9b1eb4737b344cd2d51d3ddf9d8e7889405742dea184f00fcf6107e9600f5b88`
- `modified_blackbody_20_2000K.csv`: `506600556581c66a142d1112a92b16c437101f74e0b594886ac57025104c7f9e`
- `build_contamination_tables.py`: `4866a236b587c88c36f5b49f8e4c605c755f134cb3014a67aea27107c913cc78`
- `tests/test_contamination_tables.py`: `90f4fc1b5f2a0f25ded722e8bc1be622bd203b68724ca0c342c5feeeeb659fcd`
- `latex/contam_Circinus.tex`: `f5b971e0844e9b3aac387ee737fe26df9d3e755aa4d67b1c658a2272488f2aa8`
- `latex/contam_NGC1068.tex`: `98e34ce36f7e24fff53455eec8884a46f2dc180f6463ec05040b39c6fbf7acc0`
- `latex/contam_NGC1068_central_only.tex`: `83a92ef855504fa1b97df4245651075d2ba6b9b4b5b1a2b82d8152447daa326e`
- `latex/contam_NGC1194.tex`: `0ef824b45b301a3ecc55b4f7ca0bca9bf6d726da11f350d4eca6520bd2d68e95`
- `latex/contam_NGC3079.tex`: `79c5d86e8cb2fa14075538fa0955ab86c0c8b1efae739f6f42b55dd18a00b038`
- `latex/contam_NGC4258.tex`: `ac9127d62a2a1d76e6ed60b90252b74157066cbccb36993875380ebb6abada12`
- `latex/contam_NGC4945.tex`: `9c7a7aa162a4f41caac8a3b37ed7644805e8682df9d806240637cd771a25a9c5`

- `build_sed_figures.py`: `24ec75dfe0350363a248a6e9501677d569e4f7ec8077503d21f778c93c0437ea`
- `tests/test_build_sed_figures.py`: `e6336583034ae08e8560efc4f155ee90ff13338eba4d38c54a13d7db23383446`
- `sed_figure_manifest.json`: `3424a2d99e4a72e5475635dd81c7de6d920ea103d31412affce6407de8625a69`
- `Circinus_Spectra.png`: `fe5150b8811ad2e0fc79bb61eb8be6d4d5fc6c6a5ee2dc856d1483a959a47679`
- `NGC1068_Spectra.png`: `fe4c9342e934e6934b6ef113d746e568d006e43fcedbadf87daac0cc16a7d9de`
- `NGC1194_Spectra.png`: `263d265492f568f46df7a82bbf1cbd1cc8d337386dfcdb6bdaa07d3f2346f7f2`
- `NGC3079_Spectra.png`: `60c8fa260952117c7cc176dd7af3870c8196f5cdd1381184cb259794e4cf9dfe`
- `NGC4258_Spectra.png`: `9b66500ea4e36bc49633049ff541050ddedb95e61c8aa36f026c88173fa5db89`
- `NGC4945_Spectra.png`: `0ad849f378ce5d999a2d7790117c53077b46fe74887c66c5cffe170d1a80e078`

This was the now-superseded post-IC-2560 snapshot: it followed the IC 2560
formal-only decision, deterministic publication rounding, regeneration of
Tables 1 and 2, contamination and SED table/figure products, the BHS figure,
and quarantine of the four historical notebook plotting cells.  It is not the
current release fingerprint.

### Authoritative Task 28 final-release fingerprint (2026-08-26)

This block supersedes every historical fingerprint above for packaging and
publication validation.  The isolated Task 28 CASA rerun, exact reviewed-delta
adjudication, all downstream deterministic builders, the manuscript tables and
figures, and the independent numerical audit completed successfully.  The
final catalog contains 501 measurements (485 image-based and 16 literature),
the selected-$S_{230}$ table contains 19 sources (12 detections and seven upper
limits), and the recovery manifest binds 485 local payloads totaling
6,040,874,990 bytes.

- `fitsummary_fixtochosen_sorted.txt`: `1bb7db5304ccdfa154a07f691cfa8f90ea40755f05c7c60fb7dca13c39cbf75b`
- `230316_2_final_publication_noise_audit.csv`: `3fcf7f23b7647b3e84f7e35bba7dbd7fa7708f9aeca4191aefcf793172b7d6a2`
- `reviewed_local_scaled_mad_freeze.json`: `a5df0f1ea29e0a9e714126a3b7545a4db116c57f1da2eae071ce47385b7be2e4`
- `reviewed_current_to_generated_deltas.csv`: `f3093fe1437438ade2dd96cc6aa2177fd9a37487d0d62194139f628b65953909`
- `reviewed_new_product_additions.csv`: `ba8fbe138d8236a7c1acee444a552b5c4af90646b8340bc39b8759fecca635fa`
- `reviewed_delta_adjudication.json`: `64fe1ec2799274ab951d02cf70e0867b2eb3e480e27703f46d31a073b53ea752`
- `fitsummary_final_withextras.txt`: `fe14f16400a3e880176d34202033e037630b3b8911a159e56cb9a5f65c055a3a`
- `fitsummary_machine.txt`: `7b2cc7ffcedf35bf0e791946213288a75c78b520d98f873895535b29c9b5f194`
- `fitsummary_aas_mrt.txt`: `1498cba1c6a42acea20be17ccf9c4c97e92239ac43d2a679c32af2478a38a688`
- `AAS72796R2_datafile2_corrected.txt`: `1498cba1c6a42acea20be17ccf9c4c97e92239ac43d2a679c32af2478a38a688`
- `S230_extrapolation_comparison.csv`: `010bc028f70ae8db8bb37cfc8c82e30fe2d86b2ebae2b03fbd8a08acd679c1bc`
- `variability_machine.txt`: `e2a01b7f1069472fd7c068cdcc19e12610580f50c68b35e68e79768511119314`
- `correlations_machine.txt`: `9e1d98aa45ea1bb544607c1ce6c6db229c21c7ab06fa827f6bdc70f04ed4ed59`
- `contamination_selected.csv`: `2cb0929078009c51472c759839b77f0c2ab877ea686f0049709a38388908f12d`
- `contamination_results.json`: `fce093602e6dc6806c2e425ff6e19c126ae239cf1eb1e19aa5c8ce26756cec19`
- `continuum_recovery_manifest.csv`: `150e42aa8bea8538a161469bc9cf88e47d64ae1f8e54259403916100529e6892`
- `continuum_recovery_manifest_summary.json`: `251e6cc4d23e1814f49ed50bdcfaec4cf5fa97f8f7f2675f4d8469ad63d12b3a`
- `s230_wcs_evidence.csv`: `abc586e5b280d8d5b15723fc346cbbd049e9f911e4484779c0d0bdb801a68795`
- `s230_anchor_astrometry_audit.csv`: `3e47ee06865b999b29827c1bebfbfd23f15699d919c9ed8e53f8e02395d7baa2`
- `ngc4945_nuclear_photometry.csv`: `365878f9099c075a8702bf158523782b62c265d50592e3a2e49df773057c3650`
- `ngc4945_component_audit.csv`: `034da862c91e4e4691054b01d0b15f2e1dab5080b1e1c0eb0b875269e25ba173`
- `validated_continuum_corrections.csv`: `93a6095ed704fc94a421c6bda6f81ac5b5f52c16c8aecf33ca1f8864931e3f11`
- `BHS_params.csv`: `39163060c5024f057eaec3c55c421ec1316b18f95bfd701231ebbf451fca3fc7`
- `BHS_Flux_Graph.pdf`: `52f919229f59c64124e0f2743e7924cec4ec79c5d3ad42cb2abd6375d43a19a3`
- `BHS_Flux_Graph.png`: `c829b24d64c517028576c2ebd1412442ac95e1a3306046e2934d81e1c99e53ee`
- `sed_figure_manifest.json`: `8edef9a73f737b0877eddbcbc6ec032ca545fc50895f26c111acc51b005a5528`
- `NGC4258_Spectra.png`: `bbf6e4231c4945afca6b2b6ae538775a0ebc07442fb216364ed0eabae9098f1d`
- `NGC1194_Spectra.png`: `410c05bf8d0b61a487442827385dd618ac5c96cd332407093f04672d8d66141a`
- `NGC1068_Spectra.png`: `347906bf2cf19fe7c5d980d0cc6937efd52d8008f4ded3103d3f2a92493c4e8f`
- `NGC3079_Spectra.png`: `73efddf740ce211f3c0c2cf0e38e1e75203258cefef350da147948efc027a80b`
- `NGC4945_Spectra.png`: `9c45f9a5edcdd3dacf23e09fc0e6fd48364583ef23cd935e2565d77022daaef2`
- `Circinus_Spectra.png`: `3c797c93bdf996134174853c9964f50d1c6bc13648b8586d1cc4204eb5fdad6e`
