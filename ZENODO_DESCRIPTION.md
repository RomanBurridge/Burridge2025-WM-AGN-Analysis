# Zenodo description — v1.0.11

Zenodo reserved the version-specific DOI `10.5281/zenodo.23759601` for
v1.0.11; the record becomes public when published. The previous
published v1.0.10 remains at
[10.5281/zenodo.22121490](https://doi.org/10.5281/zenodo.22121490); that DOI
identifies the older version. The all-versions DOI is
[10.5281/zenodo.18135780](https://doi.org/10.5281/zenodo.18135780).

## Description

WM-AGN Analysis Code contains the analysis software and lightweight
machine-readable products for *When the Shadow Meets Its Measure: Assessing
the Feasibility of Submillimeter Black Hole Shadow Imaging in Megamaser Disk
AGN*, accepted by ApJ as AAS72796R3 / AE8A42.

This release repairs clean installation, authenticated NGC 4945 workflow
initialization, exact archive-input selection, deterministic ALMA Stokes-I
serialization, and file/path safety. It retains the paper's adopted numerical
inputs, measurements, uncertainty treatments and conditional contamination
prescriptions. Its table renderer restores the accepted uncertainty markers;
the BHS-versus-flux figure omits error bars, as in the accepted manuscript.
Fitting and exact publication rendering use separately specified environments.

The frozen fitter inventory contains 544 discovery inputs: 222 ALMA files,
the exact 307 historical NVAS products, and 15 SMA discovery-image trees.
Six selected-clean SMA image/residual trees are additional support products,
not additional discovery inputs. The 89 admitted ALMA products remain
hash-bound. The included publication products comprise the 501-row continuum
catalog, the 21-source black-hole parameter table, the 19-source common-frequency
reference table, contamination and variability products, and publication
tables and figures. Large observations and installed environments are not
bundled; `DATA_ACCESS.md` describes their recovery.

`CANDIDATE_TEST_REPORT.txt` records the actual rc1 verification and its
remaining limitations, inherited by v1.0.11. The stable promotion changes only
version metadata and documentation; `STABLE_RELEASE_PROMOTION.txt` records
the separate packaging checks. Scientific source and products are unchanged. Exact historical fitting uses the
hash-matching historical image products. Fresh SMA archive processing remains
a separate limitation, including an unresolved NGC 4258 flux difference; this
release does not adopt a replacement measurement or claim exact recovery
of every historical SMA product from the present archive. Optional fresh MIR
conversion and reconstruction of original maser posterior chains must not be
described as verified unless separately demonstrated.

## Existing metadata to preserve unless a verified correction is required

- Title: WM-AGN Analysis Code
- Creator: Roman N. Burridge; ORCID 0009-0000-3265-7726
- License: MIT
- Version: v1.0.11
- Publication date: 2026-10-28 (Pacific/Honolulu)
- Version DOI: 10.5281/zenodo.23759601 (reserved by Zenodo)
- Keywords: black holes; megamasers; water masers; active galactic nuclei;
  VLBI; submillimeter astronomy; millimeter astronomy; SMA; ALMA; VLA; CASA


Before publication verify the actual draft, file checksums and preserved
creators/contributors, MIT license, AAS community, repository and arXiv
relationship. Cite the new version using its verified public Zenodo record.
