# Suggested Zenodo Description

This record contains the analysis code and lightweight machine-readable tables associated with the manuscript *When the Shadow Meets Its Measure: Assessing the Feasibility of Submillimeter Black Hole Shadow Imaging in Megamaser Disk AGN*.

The code supports the paper's analysis of water-megamaser disk AGN as possible targets for future submillimeter/millimeter space-VLBI black-hole-shadow imaging. It includes scripts and notebooks for SMA calibration and imaging, archival ALMA/VLA continuum-product fitting, continuum flux compilation, selection of 200--400 GHz anchor measurements, conversion of those measurements to a common 230 GHz reference frequency, black-hole-shadow angular-size calculations, VLBI baseline/sensitivity calculations, and thermal-dust/extended-jet contamination estimates.

The archive includes the full machine-readable continuum-fitting table, the 230 GHz reference-flux table used by the revised manuscript, and cleaned diagnostic variability/correlation machine tables. It does not include raw SMA data, ALMA or VLA archival FITS files, CASA image products, CASA logs, or measurement sets. The included `DATA_ACCESS.md` file describes how to retrieve the external data from the ALMA Science Archive, NRAO/NVAS archive services, and the SMA public archive.

Primary dependencies include Python, NumPy, SciPy, pandas, Matplotlib, Astropy, Astroquery, Requests, PyUVData, Jupyter, and CASA 6 with `casatasks`, `casatools`, and `analysisUtils`.

This version is intended to correspond to the final revised analysis for the AAS Journals manuscript.

This v1.0.8 release corrects the machine-readable continuum table and the software used to generate it. The corrected table restores two omitted rows, removes duplicate NGC4258 literature rows, restores the full Doi et al. reference text, uses `SMA_this_work` for the SMA rows from this work, preserves the fit-results numeric precision in the machine-readable table, and updates the AAS/CDS byte-by-byte header line for the FWHM `*` entry.

The manuscript table `tab2.txt`, the AAS production file `datafile2.txt`, and `analysis/multi_freq_from_archive/machinetables/fitsummary_machine.txt` refer to the same continuum table in different contexts. These corrections do not change the printed figures, printed science tables, or scientific conclusions of the accepted manuscript.

## Suggested Keywords

black holes; megamasers; water masers; active galactic nuclei; VLBI; submillimeter astronomy; millimeter astronomy; SMA; ALMA; VLA; CASA; Astroquery

## Suggested Community

Submit the record to the `AAS Journals` Zenodo community.
