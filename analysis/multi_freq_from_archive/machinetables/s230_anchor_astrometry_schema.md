# S230 anchor astrometry and component-identity audit

## Purpose and row grain

`s230_anchor_astrometry_audit.csv` is the machine-readable astrometry record for the exact 19 observations currently selected by `S230_extrapolation_comparison.csv`.

Each row represents one selected S230 anchor, not one literature source or one possible image component.  A row therefore carries the selected observation ID and reference verbatim so that an old astrometry decision cannot silently attach to a replacement S230 product.

The table contains 12 detections and 7 upper limits.  Upper-limit rows deliberately use `NA` for component, reference-coordinate, transform, and offset fields because no component was adopted.

## Deterministic build

From the repository root:

```bash
python3 analysis/multi_freq_from_archive/build_s230_astrometry_audit.py
python3 analysis/multi_freq_from_archive/build_s230_astrometry_audit.py --check
```

The builder joins four release records:

1. `machinetables/S230_extrapolation_comparison.csv` — current observation, status, beam, and product selection.
2. `s230_anchor_identity_audit.csv` — accepted identity disposition and required handling.
3. `s230_wcs_evidence.csv` — hash-verified WCS extracted from each actual selected image.
4. `s230_astrometry_policy.csv` — component/reference coordinates, literature provenance, availability flags, and disposition evidence.

The builder fails if the four source sets differ from the exact 19 current S230 sources; any observation ID, selected reference, or status is stale; a detection lacks coordinate/disposition evidence; an upper limit contains a component coordinate; or the NGC 4945 source18/source17 records drift.

`s230_wcs_evidence.csv` is refreshed separately with CASA Python because CASA image tables are among the selected products:

```bash
<CASA_PYTHON> analysis/multi_freq_from_archive/extract_s230_wcs_evidence.py --archive-root <ARCHIVE_ROOT>
<CASA_PYTHON> analysis/multi_freq_from_archive/extract_s230_wcs_evidence.py --archive-root <ARCHIVE_ROOT> --check
```

The WCS extractor verifies every selected product against `continuum_recovery_manifest.csv` and records the recovery ID and SHA-256 in the snapshot.  No workstation path is written to a release table.

## Coordinate frames and transforms

The common comparison frame is ICRS.

- `ICRS` native coordinates are copied without rotation.
- `FK5_J2000` native coordinates are transformed with the USNO Circular 179 frame-bias angles used by Astropy for FK5 at equinox J2000: eta0 = -19.9 mas, xi0 = +9.1 mas, and delta-alpha0 = -22.9 mas.
- The pure-Python builder implements and tests the same three-dimensional rotation, so rebuilding the final audit does not require Astropy.
- `coordinate_transform_method` records the component and reference operation separately.

The source papers generally label equatorial positions “J2000” without a separate `RADESYS` keyword.  The policy encodes those labels as `FK5_J2000` to make the transformation explicit and reproducible.  Radio positions whose published provenance explicitly ties them to an ICRF/ICRS realization remain `ICRS`; in particular, the NGC 5765B phase-reference zero point comes from the VLBA Calibrator Survey.  If original visibility metadata later demonstrate that another quoted number is directly ICRS/ICRF despite its J2000 label, change that policy row's native frame and rebuild; do not silently reinterpret the committed number.

The selected-image WCS is not inferred from telescope name.  `wcs_direction_reference_code` comes from the actual image coordinate system.  Raw FITS `RADESYS`, `EQUINOX`, and `CTYPE` are retained when a FITS header exists; CASA image-table inputs use the explicit `NA_CASA_IMAGE_TABLE` sentinel for FITS-only keywords.

## Offset convention

Offsets always mean component minus the row's reference coordinate:

- `signed_east_offset_arcsec` is positive when the selected component is east of the reference.
- `signed_north_offset_arcsec` is positive when it is north of the reference.
- `total_separation_arcsec` is the great-circle separation.
- `separation_geometric_beams` is total separation divided by the exact geometric beam recorded in the current S230 table.

The signed components use the local tangent-plane approximation for east and north; the total uses the haversine formula.  At the subarcsecond separations in this table their difference is below the reported precision.

Beam separation is descriptive, not a statistical significance.  Every freely fitted selected component lacks a preserved positional covariance, and selected-image phase-reference QA is not packaged.  This matters especially for Circinus, whose old absolute maser position uncertainty is larger than the 0.02063 arcsec ALMA beam, and for the resolved S1 region of NGC 1068.  `positional_significance_status` therefore fails closed rather than dividing an offset by an invented error.

## Field groups

| Field group | Meaning |
|---|---|
| `source` through `selected_geometric_beam_arcsec` | Exact current S230 selection keys and observing metadata. |
| `identity_disposition` through `disposition_evidence` | Adopted nuclear, ambiguous, blended, upper-limit, or conservative-upper-limit decision and its evidence. |
| `component_*` | Selected fitted/fixed component coordinate, native frame, ICRS coordinate, provenance, and positional-uncertainty availability. |
| `reference_*` and `authoritative_nuclear_coordinate_status` | Nuclear/dynamical-center coordinate or explicitly labeled proxy, its citation, transform, and uncertainty availability. |
| `signed_*`, `total_separation_arcsec`, `separation_geometric_beams` | Recomputed component-minus-reference astrometric comparison. |
| `selected_image_phase_calibrator_*` | Whether a phase-calibrator identity was recoverable for the selected image. |
| `selected_image_astrometric_qa_*` | Availability of quantitative source-level astrometric QA; no unavailable quantity is imputed. |
| `recovery_record_id` through `wcs_transform_convention` | Hash-verified selected-image WCS metadata and its pixel-to-world convention. |
| `accepted_component_record` and `rejected_component_*` | NGC 4945 source18 acceptance and preserved source17 rejection; `NA` elsewhere. |

All unavailable values use the literal `NA`, except FITS-only keywords on CASA image tables, which use `NA_CASA_IMAGE_TABLE` so “not applicable” is distinguishable from “not recovered.”

## Evidence boundaries by source

- **Direct or derived nuclear/dynamical centers:** Circinus (Greenhill et al. 2003), NGC 1068 S1 (Gallimore et al. 2004), NGC 3393 (Kondratko et al. 2008), NGC 4945 source18 (Emig et al. 2020), and NGC 5765B (Gao et al. 2016).
- **Published maser nuclear proxies rather than explicit disk centers:** ESO 558-G009 (the Gao et al. 2017 VLA-B maser correlation position), NGC 1194, and NGC 4388 (positions compiled by Kamali et al. 2017 for the latter two).
- **Ambiguous:** IC 2560 uses a published maser-position proxy whose RA uncertainty is comparable to the measured separation; NGC 1320 has only a NED query coordinate in the audited evidence, which is not independent registration proof.
- **Blended:** NGC 3079 has an absolute reference-maser position but not an absolute disk-center position, and its 5.21 arcsec beam blends the nuclear region; NGC 4258 has an accurate relative disk fit but no audited absolute reference-feature coordinate, and its 5.63 arcsec beam is likewise blended.
- **Upper limits:** Mrk 1029, NGC 2273, NGC 2960, NGC 6264, NGC 6323, UGC 3789, and UGC 6093 have no adopted component coordinate.

For NGC 4945, `accepted-source18-fixed-canonical` is the fixed AGN/maser-nucleus coordinate used by the released S230 flux.  `rejected-source17-archived-anchor` remains in the same row as explicit evidence that the SSC/SNR-associated historical anchor was not reused.

## Calibration and QA coverage

The 9 SMA rows retain their phase-calibrator names from `analysis/sma_calibration/txt_230316.txt`.  That establishes which calibrator was used, but the release does not contain a quantitative per-source astrometric QA metric, so the QA status remains unavailable.

The 10 selected ALMA anchors remain 10 selected ALMA FITS products after X600 replaces X604 for ESO 558-G009.  The released X600 member documentation identifies J0648-1744 as its phase calibrator, making ESO 558-G009 the sole ALMA calibrator-identity exception in this audit.  The other nine selected ALMA images do not retain phase-calibrator identities in their image headers, and no selected ALMA row has a complete quantitative source-level astrometric QA record; those unavailable fields remain explicit rather than imputed.
