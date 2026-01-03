## To Calibrate SMA Images

Upload SMA Observing Script with the following naming scheme:
- analysis/sma_calibration/txt_220131.txt
- analysis/sma_calibration/txt_230316.txt

Convert SMA MIR files to .ms format: 
- analysis/sma_calibration/sma2casa.py

Change SMA observations to the following naming scheme:
- analysis/sma_calibration/220131_128_master.ms
- analysis/sma_calibration/230316_128_master.ms

Molecular transitions to flag: analysis/sma_calibration/molecular_transitions.txt
	Martín, S., Mangum, J. G., Harada, N., et al. 2021,
	Astronomy & Astrophysics, 656, A46
	doi: 10.1051/0004-6361/202141567

Calibrate SMA data: analysis/sma_calibration/calibrate_sma.ipynb	
- set "date" and "binsize" to the corresponding SMA observation
- set "flagnum" to flag edges, or molecular lines
- set calnum to correspond with chosen calibrators


## Collect Archival Data from ALMA and VLA
## Produce SED Table and Plots, Variability Table and Plots, BHS Size vs Flux Plot, 
## Estimate Contamination Due to Thermal-Dust, Extended-Jets, and Intrinsic Variability,
## Estimate spectral index from SMA lower and upper side bands

Text file of coordinates to use for archival collection: analysis/multi_freq_from_archive/namesandcoords.txt

Obtain archival data from namesandcoords text file:
- get_urls.ipynb

Produce continuum peak extraction for point sources using intelligent Gaussian fitting: 
- analysis/multi_freq_from_archive/get_continuum_params.ipynb
	- Based on: https://lweb.cfa.harvard.edu/rtdc/SMAdata/process/tutorials/sma_in_casa_tutorial.html
	- Run this script 4 times in the following order, setting each successive parameter to 1 and the rest to 0: prelimchecker, fit_sizenoloc, getchoosercoords, fixtochosenfreq
		- 'prelimchecker' looks for bad images
		- 'fit_sizenoloc' fits the gaussian without setting a bounding region
		- 'getchoosercoords' looks for an image with the best resolution in a given frequency range
		- 'fixchosenfreq' fits all observations to the chosen image from getchoosercoords

Add extra points and produce BHS size vs Flux plot, SED plots, Variability Plots, and estimate contaminations due to thermal-dust and extended-jets: 
- analysis/multi_freq_from_archive/post_processing_finalcontinuum.ipynb
	- executes the following: Beam_Date_Flux_Graphs.ipynb, Variability_Search.ipynb, flux_size_graph.ipynb, and get_contams.ipynb under the same parent directory

Upon obtaining Variability Plots, to produce Variability Estimates: 
- analysis/multi_freq_from_archive/Contam_variability.ipynb

Spectral Index Extraction using SMA Upper and Lower Side-band Images: 
- analysis/multi_freq_from_archive/get_alpha.ipynb

## Calculations of BHS Parameters and VLBI_sensitivity

BHS Parameters: 
- BH_Parameters.ipynb

VLBI_sensitivity: 
- VLBI_sensitivity.ipynb