# Release verification checklist — v1.0.11

The previous v1.0.10 record at `10.5281/zenodo.22121490` remains preserved.
v1.0.11 has its own version-specific DOI assigned by Zenodo. Keep the archive,
internal version metadata, upload label and citation consistent.

The actual scientific validation below was performed for v1.0.11-rc1.
`CANDIDATE_TEST_REPORT.txt` retains that report verbatim. v1.0.11 changes only
version metadata and documentation; `STABLE_RELEASE_PROMOTION.txt` records
its separate version, integrity and scientific-payload byte-identity checks.
Do not describe the inherited tests as a fresh stable-version scientific run.

## Inherited scientific validation requirements

1. **Clean installation.** Follow `INSTALL.md` in a new CPython 3.12/Linux
   environment. Record the exact constraints, official CASA wheel, dependency
   check, isolated runtime setup and actual fitter/Measures import probe.
   Keep the exact-publication renderer separate. No inherited research
   environment or private cache may substitute for the documented setup.

2. **Exact data membership.** Validate the frozen 550-record input manifest:
   544 discovery inputs (222 ALMA, 307 historical NVAS, 15 SMA trees) and
   six selected-clean SMA image/residual support trees. Admit no later query
   extras silently. Retain all 89 ALMA admission identities and exact file
   hashes. Verify recovery, actual public-parent Stokes-I regeneration, and
   rejection of changed pixels, WCS, beams or parent identity. A 501-row
   final catalog is not a substitute for the full discovery manifest.

3. **Preparation and actual fitting.** Use a new isolated workspace and
   independent copies of the observations. Production preparation must
   initialize authentic NGC 4945 historical 27/66 records. Actually run the
   four fitting passes and the 45-image fixed-nucleus correction, then prove
   the 66-to-84 ledger transition and every integrity, identity, noise and
   classification gate. Preserve the scientific review/adjudication stage.
   Do not replay shipped final fitted values or weaken safeguards to claim
   a successful fresh fit. Preserve failed-run evidence.

4. **Numerical regression.** Compare all 495 fit rows, the 501-row final
   catalog, all 19 reference choices (12 detections and seven upper limits),
   adopted BHS values and uncertainty semantics, contamination choices and
   downstream products with the accepted private manuscript. Record every
   numerical, classification, membership or NA-status difference. Separate
   tiny non-adopted diagnostics from adopted measurements explicitly; do
   not apply a blanket tolerance that hides a scientific difference.

5. **Presentation.** Regenerate and inspect the seven PNGs and BHS PDF in
   the specified renderer. Match the accepted figure's omission of error
   bars and Table 1's four uncertainty daggers and explanatory note. Verify
   notebook examples against the accepted 0.075-microarcsecond benchmark.
   Explain intentional changes from the old public artifacts rather than
   blindly replacing golden files. Keep numeric scientific inputs unchanged.

6. **Safety and limitations.** Test refusal of pre-existing outputs, normal
   source/path handling, partial failures, concurrent recovery destinations,
   and SMA channelization before calibration. Retain the historical
   calibration/mask behavior unless a separately reviewed scientific change
   is approved. Disclose the unresolved fresh-SMA NGC 4258 difference. Do not
   replace its adopted flux, add an invented systematic, tune to the old
   answer, or claim unexecuted MIR/posterior-chain reconstruction succeeded.

7. **Self-contained release.** Include all required lightweight source,
   manifests, historical initialization assets, tests, installation/data
   instructions, changelog, actual test report and payload checksums. Exclude
   raw observations, FITS/MS/CASA trees, runtime installations, caches,
   private manuscripts/emails, credentials and validation workspaces.
   Inspect the ZIP and independently extract/check it. Record its filename,
   byte count and SHA256. `CANDIDATE_TEST_REPORT.txt` must distinguish completed
   tests, intentional skips, failures and unavailable evidence; a proposed
   command or green saved-output check is not a completed raw fitting run.

## Publication checks

Verify the final stable ZIP, checksum, concise change summary, original test
report and inherited limitations. Confirm the archive root and `CITATION.cff`
identify v1.0.11, the internal checksum inventory passes after a fresh
extraction, and scientific source/products remain byte-identical to rc1.

Use a new version of the existing Zenodo publication. Preserve v1.0.10 and
its files. Upload only the stable ZIP and its matching checksum. Preserve
creators/contributors, ORCIDs, MIT license, AAS community, languages,
development status, repository and arXiv relationship. Use the actual
publication date and new version DOI; do not invent a journal DOI.

Publish only with the author's authorization. After publication verify the
public record and downloaded archive checksum. Update the paper's software
citation using the actually published version-specific DOI and date.
