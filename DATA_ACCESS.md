# External Data Access

This Zenodo release intentionally excludes archival data files and generated image products. The external data can be retrieved from the original observatory archives.

## ALMA

ALMA continuum products can be retrieved from the ALMA Science Archive:

- Web archive: https://almascience.org/aq/
- Regional archive mirrors: https://almascience.eso.org/aq/ , https://almascience.nao.ac.jp/aq , and https://almascience.nrao.edu/aq/
- Programmatic access: `astroquery.alma`

The machine-readable continuum table lists full ALMA product filenames. Rows beginning with `member.uid___...` include the ALMA Member OUS identifier and product name. These identifiers can be used in the ALMA Science Archive query interface or through `astroquery.alma` to locate the corresponding products. ALMA also provides a Request Handler download script for command-line retrieval after selecting products in the archive.

Minimal programmatic sketch:

```python
from astroquery.alma import Alma

alma = Alma()
result = alma.query_object("NGC 1068", public=True)
print(result.colnames)
```

## VLA / NVAS

VLA data and products can be retrieved through NRAO:

- NRAO Data Archive / Archive Access Tool: https://data.nrao.edu
- NRAO VLA archive documentation: https://science.nrao.edu/facilities/vla/archive/index
- NVAS historical VLA image products: https://www.vla.nrao.edu/astro/nvas/
- Programmatic NVAS image-list access: `astroquery.nvas`

The machine-readable continuum table lists NVAS image filenames for historical VLA products. These filenames can be matched through the NVAS service. Raw and modern VLA data can be located through the NRAO Archive Access Tool; NRAO also provides TAP-based scripted metadata queries, while downloads are completed through the archive web/request workflow.

Minimal NVAS sketch:

```python
from astroquery.nvas import Nvas

urls = Nvas.get_image_list("12h18m57.5s +47d18m14s", radius="0d0m30s")
for url in urls[:5]:
    print(url)
```

## SMA

SMA data can be retrieved from the SMA public archive:

- SMA Data Archive: https://lweb.cfa.harvard.edu/rtdc/SMAdata/archives/
- SMA archive instructions: https://lweb.cfa.harvard.edu/rtdc/SMAdata/instructions/archive_use.html

The SMA observations used here are associated with programs `2022B-H002` and `2021B-H004`. Search the SMA archive by project code and observing date. The archive provides public data after the proprietary period and can provide raw MIR data or CASA measurement-set formats when available.

This release includes the observing-script text files used by the calibration scripts:

- `analysis/sma_calibration/txt_230316.txt`
- `analysis/sma_calibration/txt_220131.txt`

Raw SMA interferometry directories contain the science targets and calibration sources observed during the night; keep that structure intact when placing data under `analysis/sma_calibration/sma_data/`.
