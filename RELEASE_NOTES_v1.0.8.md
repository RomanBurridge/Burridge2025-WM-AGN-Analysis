# v1.0.8 Correction Release

This release corrects the machine-readable continuum table and the software used to generate it.

The manuscript table `tab2.txt`, the AAS production file `datafile2.txt`, and `analysis/multi_freq_from_archive/machinetables/fitsummary_machine.txt` refer to the same continuum table in different contexts.

Corrected files include:

- `analysis/multi_freq_from_archive/build_fitsummary_machine.py`
- `analysis/multi_freq_from_archive/fitsumfiles/sorted/fitsummary_fixtochosen_sorted.txt`
- `analysis/multi_freq_from_archive/fitsumfiles/final/fitsummary_final_withextras.txt`
- `analysis/multi_freq_from_archive/machinetables/fitsummary_machine.txt`
- `analysis/multi_freq_from_archive/machinetables/fitsummary_aas_mrt.txt`
- `analysis/multi_freq_from_archive/machinetables/AAS72796R2_datafile2_corrected.txt`
- `analysis/multi_freq_from_archive/machinetables/fitsummary_aas_mrt_changes.txt`
- `analysis/multi_freq_from_archive/machinetables/README.md`

The corrections do not change the printed figures, printed science tables, or scientific conclusions of the accepted manuscript.
