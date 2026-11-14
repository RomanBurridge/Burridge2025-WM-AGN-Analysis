# WM-AGN v1.0.11-rc1: unpublished repair candidate

This candidate is for the accepted private ApJ manuscript AAS72796R3 /
AE8A42. It is not a new scientific analysis, an update to the older arXiv
paper, or a published Zenodo release. The original v1.0.10 archive and
observations are preserved. No new version DOI is asserted.

## Repairs

- Added the missing `reproject` dependency and tested Python/CASA version
  pins, including a compatible CPython 3.12 Linux CASA wheel.
- Added a pinned analysisUtils installer and an isolated CASA launcher with
  explicit data, configuration and cache locations. No private workstation
  installation or user CASA configuration is required.
- Initialized NGC 4945 runs from authenticated historical 27-row photometry
  and 66-row correction assets. These are historical metadata, not a replay
  of newly fitted values; actual fitting, correction and final guards remain
  mandatory. New preparation also works after legitimate final diagnostics
  have been promoted.
- Froze the exact 544 discovery inputs and six selected-clean supports.
  Preparation verifies their hashes and does not silently add new query
  results. Public recovery includes the five ALMA Stokes-I transformations.
- Made those five transformations deterministic at the FITS serialization
  layer while retaining the original pixel, coordinate, beam and full-file
  integrity checks. No admission threshold or scientific policy was relaxed.
- Corrected fitter output paths to follow the actual input source rather
  than inherited loop state, including the NGC 5765b naming case.
- Added no-clobber recovery and MIR-conversion behavior. Verified downloads
  and derived files are installed without replacing a destination created
  concurrently; incomplete user downloads are not deleted.
- Added a read-only SMA spectral-window/channel preflight before calibration.
  Calibration assumptions, averaging, line exclusions and adopted fluxes
  have not been changed.
- Restored the four accepted Table 1 uncertainty markers and associated
  note; omitted both horizontal and vertical uncertainty bars in Figure 1
  without removing uncertainty data; aligned the sensitivity notebook's
  benchmark with the accepted value.
- Aligned the printed continuum excerpt with the accepted seven-column
  Table 2, retaining provenance in the machine-readable/recovery products.
  The ten displayed scientific rows and numerical values are unchanged.
- Supplied a separately specified publication-rendering environment and a
  read-only rendering wrapper. Plot appearance is not assumed independent
  of Matplotlib, fonts or FreeType versions.

## Scientific scope and limits

The adopted mass/distance inputs, uncertainty treatment, spectral slopes,
source/component decisions, and contamination prescriptions are unchanged.
Exact fitting from preserved historical products and fresh reduction of a
new archive delivery are different tests. The independent NGC 4258 SMA
discrepancy remains separate and has not replaced the paper's flux.

Large telescope products and preinstalled environments are not in the ZIP.
Fresh MIR-to-MS conversion and exact recovery of all historical SMA products
are not claimed from installation or unit tests. File/hash checks prove input
identity, not the completeness of an astrophysical model; the paper's stated
conditional dust/jet and disk-model limitations remain.

See `CANDIDATE_TEST_REPORT.txt` for actual integration results and any
differences, `INSTALL.md` for installation, and `DATA_ACCESS.md` for inputs.
The historical `CORRECTION_CHANGELOG.md` is retained as the v1.0.10 scientific
change record; it is not evidence that this candidate's new execution passed.
