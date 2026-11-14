# ALMA and VLA archive-file selection audit

Status: automated provenance and selection screening complete; every scientific disposition implemented in the authoritative Task 28 rollout  
Audit date: 2026-08-25  
Manuscript, catalog, figures, and release staging changed: **complete 501-row continuum release synchronized**

## Bottom line

The cryptic filenames did not cause a wrong source, frequency, date, or unit to survive silently in the current catalog. Every local FITS header is readable, and the downstream fitter obtains frequency, date, beam, units, and sky coordinates from FITS metadata rather than from the filename.

The important defect is omission. The historical ALMA acquisition notebook treated filename tokens as scientific rejection rules. In particular, it rejected `tt0`, `mfs`, and `manual`, even though ALMA uses `tt0` for the ordinary intensity image in a multi-term product, `mfs` for a multifrequency-synthesis image, and `manual` to describe how the image was produced. It also rejected otherwise valid continuum products containing `split.cal`. The earlier fail-closed archive census inherited those assumptions and therefore was not complete.

**Disposition implemented 2026-08-26:** filenames no longer determine scientific eligibility. All 49 pre-cutoff filename-only omissions, all 26 straightforward post-cutoff target continuum datasets, and one Stokes-I measurement from each of the seven MFS/polarization deliveries were admitted. The regenerated catalog retains one full-band measurement per source/MOUS, preferring an aggregate continuum image, then a full-band `tt0` intensity image, and then a suitable target Stokes-I MFS image; alternate reductions and per-spectral-window products are not counted separately. For five IQUV MFS deliveries, only the Stokes-I plane was extracted into a provenance-linked two-dimensional FITS image; polarization-amplitude products remain outside the continuum catalog. Every admitted product passed the standard fitting and visual-quality checks.

## Files already downloaded

### Pre-Task-28 working-release files

- 447 local FITS products were inspected: 140 ALMA and 307 VLA/NVAS.
- All 447 headers were readable and all use a supported Jy-per-beam brightness unit.
- All 140 ALMA files are joined to an official observation identifier. The one legacy NGC 4945 filename containing an execution-block UID is mapped through the recovery manifest to `uid://A002/X5a9a13/X528`.
- Every ALMA file contains the adopted source coordinate. All 140 selected products are aggregate continuum images. Of these, 138 have a standard PB-corrected filename; the two older naming exceptions have valid headers and established provenance.
- A fresh, **serialized** 2-arcsec NVAS coordinate query reproduces all 307 local VLA filenames exactly. The legacy NVAS CGI must not be queried concurrently because simultaneous POST requests can return crossed result sets.
- All 248 retained VLA images contain the intended coordinate. The eight images whose WCS does not contain the intended coordinate are among the 59 inputs excluded before the final VLA catalog.
- The eight outside-WCS files form two byte-identical NVAS products copied into several source directories. None survives in the publication catalog.
- One VLA filename/header frequency difference greater than 2% occurs in an already excluded NGC 3079 image. No retained row depends on that filename frequency.

### Thirteen Task 28 ALMA downloads

All 13 renamed FITS files in `new_alma_fits/` have a unique match to an official ALMA DataLink product using source, observing time, exact byte size, and FITS metadata. All contain the intended coordinate and have ALMA/Jy-per-beam headers. Their mapping is frozen in `task28_downloaded_alma_provenance_audit.csv`.

## Resolved current-file review flags

These cases were reviewed explicitly after the VLA bandwidth-smearing decision:

- The Circinus product `member.uid___A001_X133d_X114.SN1996cr_sci.spw4_6_16_18.cont.I.pbcor.fits` is retained as an explicitly off-axis, PB-corrected measurement: the image targets SN 1996cr, 23.8 arcsec from the nucleus, but contains the adopted Circinus coordinate.
- Twelve ALMA images have a nominal phase/image center more than 2.1 arcsec from the adopted source coordinate. Eleven are source-named mosaics or offset fields whose PB-corrected image contains the nucleus; the SN 1996cr product carries the separate off-axis disposition above.

## Implemented VLA off-axis decision

The retained NVAS products had already passed through primary-beam correction. The per-run AIPS processing logs record the `PBCOR`--`RADCOR`--`NVAS08` sequence; the exported FITS headers do not preserve that full processing history, so applying another primary-beam correction would double-correct the images.

Bandwidth smearing was audited separately with
`x=(channel_width/center_frequency)*(source_offset/radial_beam_FWHM)`. Raw VLA FQ tables establish the relevant channel widths, and the adopted catalog rule is `x <= 0.5`. Nine measurements exceed the threshold: four for NGC 3079, four for NGC 4258, and one for UGC 6093. They were removed from the catalog and every downstream pair table and SED. All 248 remaining VLA catalog measurements satisfy the threshold under the audited channel widths or the conservative 50 MHz maximum used where appropriate. The exact inputs and results are preserved in `VLA_BANDWIDTH_SMEARING_AUDIT.md`, `vla_bandwidth_smearing_exclusion_policy.csv`, and `vla_bandwidth_smearing_exclusions.csv`.

## Historical ALMA logic

The notebook first required the literal substring `cont`, then removed products containing category or chemical substrings. It chose among survivors using source-number fragments and the literal preferences `pbcor`, `ari_l`, and `regcal`. If more than one product remained, the workflow printed the ambiguity but later used the first item.

The post-download coordinate safeguard was not functional: it compared the wrong nested source field, reused stale variables in the VLA loop, and was guarded by inconsistent `downloadyn` values. The independent FITS/WCS audit above replaces that safeguard.

Official ALMA documentation describes the relevant product terms:

- [ALMA delivered calibration and imaging products](https://help.almascience.org/kb/articles/pdf/what-calibration-and-imaging-products-will-be-delivered-to-me)
- [ALMA QA2 data products](https://almascience.nrao.edu/portal/documents-and-tools/cycle7/alma-qa2-data-products-for-cycle-7)

Official NVAS documentation defines the otherwise cryptic VLA position-image filename fields and states that the directory coordinate is the observed field/phase center:

- [NVAS readme guide](https://www.vla.nrao.edu/astro/nvas/read.shtml)
- [VLA primary-beam field of view](https://science.nrao.edu/facilities/vla/docs/manuals/oss/performance/fov)

## ALMA products absent from the local release

The live TAP/DataLink audit queried 287 public science source/MOUS intersections and classified 18,714 FITS products. Product-level counts are large because one observation can deliver many spectral-window, calibration, polarization, and alternate-reduction images. Collapsing to one review row per source/MOUS leaves 88 source/MOUS pairs with a continuum-intensity product but no local file:

| Screening class | MOUS pairs | Meaning |
|---|---:|---|
| Pre-cutoff possible omissions | 49 | Released by 2025-01-27 and hidden by the old filename rules. |
| Post-cutoff straightforward candidates | 26 | A target aggregate-continuum or `tt0` image is available. Thirteen are the already downloaded Task 28 images. |
| Manual MFS/polarization review | 7 | A target MFS product exists, but there is no ordinary one-image Stokes-I aggregate product. |
| Off-target coordinate intersections | 6 | Supernova/other-target observations whose image footprint intersects the adopted galaxy coordinate. |

The 49 pre-cutoff possible omissions consist of:

- 40 sub-200 GHz target `tt0` intensity images;
- five NGC 4388 aggregate-continuum images near 100, 221, 263, 300, and 349 GHz hidden by `split.cal`;
- the 223.7 GHz high-resolution NGC 1068 manual continuum image `uid://A001/X879/X2d7`;
- three aggregate continuum images above 400 GHz: one Circinus and two NGC 1068 datasets.

The 26 straightforward post-cutoff candidates consist of the 13 already downloaded Task 28 images plus 13 newly recovered `tt0` intensity datasets. Two of the latter overlap 200--400 GHz, both for Circinus; eleven are below 200 GHz.

The seven MFS/polarization cases contribute one reviewed Stokes-I measurement per source/MOUS. They include two polarization observations previously summarized too simply as calibrator-only products: Circinus `uid://A001/X2d20/X3580` and NGC 4945 `uid://A001/X3621/X965`. Both DataLink records contain target IQUV MFS products as well as calibrator products; the rollout extracts only the target Stokes-I plane.

Complete one-row-per-MOUS evidence is in `alma_archive_mous_candidate_review.csv`; every product and historical rule outcome is in `alma_archive_fits_product_audit.csv`.

## Replacement selection rules

Future acquisition should use the following explicit policy:

1. Query public science observations by source coordinate and retain the source/MOUS identity from TAP.
2. Expand the complete DataLink product list and identify the science target using archive target metadata plus product target labels. A coordinate intersection alone is insufficient.
3. Select a PB-corrected Stokes-I continuum-intensity product. Treat aggregate `cont` and intensity `tt0` as eligible product roles. Treat `manual` as provenance, not as a rejection.
4. Do not automatically turn every per-spectral-window `mfs` image into a catalog row. When no aggregate continuum/`tt0` image exists, review whether a science-target aggregate MFS Stokes-I product is available and comparable.
5. Do not prefer QA2, ARI-L, regular-calibration, or self-calibration products solely from a filename. Record the reduction origin and make any preference explicit for that MOUS, preserving astrometry when localization matters.
6. Read frequency, observing date, beam, units, celestial frame, and source coverage from the FITS header. Preserve the official URL, observation identifier, exact byte size, and SHA-256 digest.
7. Query NVAS serially, retain every returned image for screening, reject files whose WCS excludes the source, and evaluate phase-center offset, primary-beam response, and smearing before admitting an off-axis measurement.

## Completed review order

1. The retained Circinus/SN 1996cr ALMA row and all off-axis VLA rows were adjudicated.
2. The six pre-cutoff aggregate-continuum omissions in the 200--400 GHz band were reviewed.
3. The 40 pre-cutoff sub-200 GHz `tt0` datasets and three pre-cutoff products above 400 GHz were reviewed by source.
4. The 13 newly discovered post-cutoff `tt0` datasets and seven MFS/polarization-only cases were reviewed.
5. The final approved inventory was processed by the isolated four-pass fitter and propagated through the catalog, selected $S_{230}$ values, tables, figures, manuscript, and release metadata.

The resulting release contains 501 continuum measurements: 485 image-based rows and 16 literature rows. The 498 measurements belonging to the 21-source sample generate 14,193 within-source pairs and a 105-pair strict diagnostic screen. The 19-source $S_{230}$ comparison contains 12 detections and seven upper limits from 11 ALMA, seven SMA, and one NOEMA measurement.
