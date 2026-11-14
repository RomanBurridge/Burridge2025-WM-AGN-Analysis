# ALMA continuum-inventory gap audit

Snapshot date: 2025-01-27 (inclusive). Audit performed 2026-08-20 against the public ALMA Science Archive TAP and DataLink services.

## Scope and query

The audit used the 22 coordinates in the preserved `namesandcoords.txt` (SHA-256 `71640aff436969cb4bc1d3a99f5bb29dc585608e6fbb0c677ab9b1d731bb5687`) and a 2 arcsec cone at each coordinate. The TAP endpoint was `https://almascience.eso.org/tap/sync`. The ADQL template was:

```sql
SELECT DISTINCT proposal_id, member_ous_uid, target_name, obs_release_date,
       schedblock_name, spatial_resolution, t_min, t_max, frequency,
       band_list, antenna_arrays
FROM ivoa.obscore
WHERE science_observation='T'
  AND data_rights='Public'
  AND INTERSECTS(
        CIRCLE('ICRS', <source_ra_deg>, <source_dec_deg>, 0.000555556),
        s_region
      )=1
```

Rows with `obs_release_date` later than 2025-01-27 were excluded. Distinct member OUS identifiers were compared with the 133 ALMA products in the preserved 455-input inventory. DataLink manifests for every unmatched MOUS were expanded through the ESO service, with the NRAO service as a failover. The executable audit logic is `archive_query_audit_executed.py`. The source used for the live execution had SHA-256 `bd738d0bd9ac965126bf9d09f8ab698dd802aacfbbe7079b8e3682f4b7fa2495`; the bundled copy changes only its absolute input/output defaults to release-relative paths so it remains portable.

## Fail-closed census

There are 112 distinct pre-cutoff unmatched source/MOUS pairs. Every DataLink expansion completed: 7 have an eligible non-line, non-Taylor, primary-beam-corrected continuum intensity FITS and 105 do not. The complete disposition is `pre20250127_unmatched_mous_disposition.csv`; it contains no `unknown` classifications. Exclusion counts are: 46 only Taylor-term or other non-intensity continuum FITS, 36 no continuum-named FITS, 13 no FITS products, 9 excluded by the existing filename-category rules, and 1 no eligible pbcor continuum intensity FITS.

The seven admitted primary products are frozen in `admitted_primary_products.csv`. The FITS payloads are deliberately not bundled in this lightweight evidence directory; the manifest preserves each official URL, release-relative destination, byte size, and SHA-256 digest needed to recover and verify it. The recorded product hashes therefore describe external archive inputs, not files covered by this directory's `EVIDENCE_SHA256SUMS.txt`.

Primary archive evidence is the TAP metadata above, the fully expanded DataLink manifests represented in the 112-row census, and the direct official product URLs in the seven-row manifest. For X600, the official NRAO delivery bundle is `https://almascience.nrao.edu/dataPortal/2013.1.00824.S_uid___A001_X122_X600_001_of_001.tar` (SHA-256 `ded936ad689664ac123b1780304d600d26a321b3dba825846a431424b1a978a2`) and its official README is `https://almascience.nrao.edu/dataPortal/member.uid___A001_X122_X600.README.txt`.

## Root cause and product identity

Six eligible products have `ESO` in the target/product URL (two ESO558-G009 and four Circinus/ESO97 products). The historical notebook applied a naive molecular-line exclusion using substring membership over the full URL and included the token `SO`. Therefore the `SO` inside `ESO` falsely excluded all six.

The bundled pre-remediation implementation is `get_urls_pre_boundary_fix.ipynb` (SHA-256 `bc898ea9d487d042b01513797d729cb49c96e52bc86a80472b4a67768ad51c22`). Its `chemicals` list contains `SO`, and its selection loop evaluates `if chem in url` before accepting URLs containing `cont`. This snapshot is historical evidence, not the current discovery implementation.

The current corrected discovery notebook is the sibling `../get_urls.ipynb` (SHA-256 `34f0da16ffcff43e1d4c2e20cd929f5875c9a1c2d8fd6de6b74d87fb0fe1127e` at packaging). It applies chemical screening only to `os.path.basename(url)` and treats the short `SO` token with a filename-boundary regular expression, so the `SO` letters embedded in `ESO` no longer reject a product.

NGC1068 `uid://A001/X122/X4d8` passes the literal historical filename filter. Its absence is an independent historical selection/download omission; no scientific exclusion is supported by the preserved rules.

ESO558-G009 X600 and X604 are distinct observations, not duplicate reductions: X600 is scheduling block `ESO_558-_a_07_TE`, ASDM `uid://A002/X96bfab/X12a6`, MJD 57004.284625--57004.298670 and the approximately 311.8--314.4 GHz tuning; X604 is scheduling block `ESO_558-_b_07_TE`, ASDM `uid://A002/X96bfab/X149a`, MJD 57004.322492--57004.334983 and a distinct approximately 315--319 GHz tuning.

For X600, the original delivered QA2 product is the sole admitted primary and the reviewed chosen-frequency reference. Its historical getchooser reconstruction has S/N 7.85306575176 and geometric beam 0.7181698841 arcsec, beating the pre-augmentation X604 anchor beam of 0.87561 arcsec. The same-MOUS ARI-L product is preserved only as a non-row alternate reduction and must not be inserted into the publication inventory.

## Minimal remediation

Insert exactly the seven primary FITS at the release-relative paths in `admitted_primary_products.csv`, then run the complete augmented four-pass pipeline and all downstream regeneration/validation. A direct insertion into final tables is not scientifically sufficient because X600 changes the deterministic chosen-frequency anchor for ESO558-G009, which changes the fixed-position/fixed-beam context and downstream derived artifacts. Re-downloading or re-screening all 455 original inputs is unnecessary: the 112-MOUS fail-closed census freezes the missing archive surface, while the augmented run should operate on 462 inputs (140 ALMA, 307 VLA, 15 SMA).

## Self-check

From the software-release root, verify the complete lightweight bundle with:

```bash
(cd analysis/multi_freq_from_archive/archive_inventory_audit && sha256sum --check EVIDENCE_SHA256SUMS.txt)
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest analysis.multi_freq_from_archive/tests/test_archive_inventory_audit.py -v
```

## Limitations

- The cutoff is an inclusive archive-publication cutoff (`obs_release_date <= 2025-01-27`), chosen to reproduce the preserved local download snapshot. It is not stated in the manuscript.
- Archive services are live. This audit freezes returned metadata and product dispositions, but the provider may later revise DataLink manifests or reprocess files.
- ESO DataLink completed for all 112 rows; the NRAO endpoint was configured only as failover. This is mirror-complete for the executed audit, not a claim that future mirror contents will remain byte-identical.
- The seven-row manifest records the archive-audit disposition at the time of the census. Final measurement classifications and anchor adoption are governed by the separately frozen augmented four-pass evidence.
