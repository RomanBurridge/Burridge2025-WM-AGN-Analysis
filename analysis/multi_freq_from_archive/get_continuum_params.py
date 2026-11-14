#!/usr/bin/env python
# coding: utf-8

import csv
import hashlib
import os
import tempfile
from pathlib import Path

from selected_clean_manifest import (
    load_selected_clean_manifest,
    normalize_source_name,
    selected_image_path,
)
from continuum_output_paths import fit_output_directory
from continuum_units import (
    brightness_to_mjy_per_beam,
    brightness_to_mjy_per_beam_factor,
    flux_density_to_mjy,
)
from continuum_noise_inventory import (
    FORCED_FIT_POLICY_CONVERGED,
    FORCED_FIT_POLICY_NONCONVERGED,
    GENERIC_FINAL_NOISE_AUDIT_FIELDS,
)
from continuum_boundary_rejections import (
    load_boundary_rejection_policy,
    validate_boundary_rejection_record_identity,
)
from continuum_measurements import (
    FINAL_NOISE_APERTURE_BMAJ_MULTIPLIER,
    FINAL_NOISE_MIN_INDEPENDENT_BEAMS,
    SCALED_MAD_NORMAL_FACTOR,
    THREE_SIGMA_FROM_FWHM,
    angular_value_to_arcsec,
    background_rms_snr,
    calendar_date_sort_key,
    classify_ordinary_fit,
    effective_independent_noise_samples,
    final_noise_review_flags,
    forced_fit_diagnostic,
    format_two_decimal,
    gaussian_beam_area_arcsec2,
    independent_beam_areas,
    normalize_right_ascension_degrees,
    portable_release_product,
    require_fixed_restoring_beam_shape,
    require_minimum_independent_beams,
    square_aperture_coverage,
    square_bmaj_region_width_arcsec,
)
from celestial_frames import (
    ICRS,
    J2000,
    casa_sexagesimal_declination,
    colon_sexagesimal_declination,
    coordinate_pair_for_reference_code,
    coordinate_candidate_attempt_limit,
    discard_failed_coordinate_candidate_or_raise,
    normalize_celestial_reference_code,
)
from chosen_frequency_anchors import (
    SELECTION_RULE as CHOSEN_FREQUENCY_SELECTION_RULE,
    consolidate as consolidate_chosen_frequency_anchors,
    path_sha256,
    sexagesimal_dec_to_deg,
    sexagesimal_ra_to_deg,
    tree_sha256,
    write_manifest as write_chosen_frequency_manifest,
)


def load_reviewed_exclusion_policy(policy_path):
    """Return the exact product-to-hash map from the reviewed CSV policy.

    Exclusions are scientific decisions, so a missing, malformed, duplicated,
    or incompletely documented policy must stop the fitter rather than silently
    changing which archive products are measured.
    """

    required_fields = {"exclusion_id", "product", "sha256", "decision_reason"}
    policy_path = Path(policy_path)
    if not policy_path.is_file():
        raise FileNotFoundError(
            f"Reviewed continuum-product exclusion policy is missing: {policy_path}"
        )
    with policy_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing_fields = required_fields - set(reader.fieldnames or ())
        if missing_fields:
            raise ValueError(
                f"{policy_path}: missing required exclusion fields "
                f"{sorted(missing_fields)}"
            )
        exclusions = {}
        seen_ids = set()
        for lineno, row in enumerate(reader, 2):
            exclusion_id = (row.get("exclusion_id") or "").strip()
            product = portable_release_product((row.get("product") or "").strip())
            expected_sha256 = (row.get("sha256") or "").strip().lower()
            decision_reason = (row.get("decision_reason") or "").strip()
            if not exclusion_id or not product or not decision_reason:
                raise ValueError(f"{policy_path}:{lineno}: incomplete exclusion record")
            if exclusion_id in seen_ids:
                raise ValueError(
                    f"{policy_path}:{lineno}: duplicate exclusion_id {exclusion_id!r}"
                )
            if product in exclusions:
                raise ValueError(
                    f"{policy_path}:{lineno}: duplicate excluded product {product!r}"
                )
            if len(expected_sha256) != 64 or any(
                character not in "0123456789abcdef" for character in expected_sha256
            ):
                raise ValueError(
                    f"{policy_path}:{lineno}: invalid SHA-256 for {product!r}"
                )
            seen_ids.add(exclusion_id)
            exclusions[product] = expected_sha256
    if not exclusions:
        raise ValueError(f"{policy_path}: reviewed exclusion policy is empty")
    return exclusions


def verify_reviewed_exclusion_hash(product_path, expected_sha256):
    """Fail closed unless an excluded input is the reviewed immutable artifact."""

    product_path = Path(product_path)
    if not product_path.is_file():
        raise FileNotFoundError(product_path)
    digest = hashlib.sha256()
    with product_path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    actual_sha256 = digest.hexdigest()
    if actual_sha256 != expected_sha256:
        raise ValueError(
            f"Reviewed exclusion hash mismatch for {product_path}: "
            f"expected {expected_sha256}, found {actual_sha256}"
        )

fundir1 = str(Path(__file__).resolve().parent)
os.chdir(fundir1)

# In[ ]:
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--flagnum", type=int, required=True)
parser.add_argument("--calnum", type=int, required=True)
parser.add_argument("--mode", type=str, required=True, choices=["prelimchecker", "fit_sizenoloc", "getchoosercoords", "fixtochosenfreq"])
parser.add_argument("--chooser", type=int, default=0)
parser.add_argument("--choose_science", type=str, default="NGC4258")
parser.add_argument("--printimages", type=int, default=0)
parser.add_argument("--writesummary", type=int, default=1)
parser.add_argument("--nofixmultval", type=float, default=10)
parser.add_argument("--printbadfix", type=int, default=0)
parser.add_argument("--telechooser", type=int, default=0)
parser.add_argument("--telechoice", type=str, default="SMA")
parser.add_argument("--smapath", nargs="+", required=True)
parser.add_argument(
    "--selected_clean_manifest",
    type=str,
    default=str(Path(__file__).resolve().parents[1] / "sma_calibration" / "selected_clean_iterations.csv"),
)
parser.add_argument(
    "--use_selected_clean_manifest",
    action="store_true",
    help=(
        "Use the manifest-selected clean SMA image for matching tracks. "
        "Enable this only for the clean-image fitting pass."
    ),
)
parser.add_argument("--moveniter", type=int, default=0)
parser.add_argument("--niternum", type=int, default=0)
parser.add_argument("--datestr", type=str, default="")
parser.add_argument(
    "--enable_notifications",
    action="store_true",
    help="Opt in to external progress notifications (disabled by default).",
)

args = parser.parse_args()
moveniter = bool(args.moveniter)
niternum = args.niternum
datestr = args.datestr.strip()

smapaths = args.smapath
selected_clean_manifest = args.selected_clean_manifest.strip()

flagnum = args.flagnum
calnum = args.calnum

summary_prefix = f"{datestr}_{flagnum}" if datestr else f"{flagnum}"

nofixmultval = args.nofixmultval
printbadfix = bool(args.printbadfix)

prelimchecker = 1 if args.mode == "prelimchecker" else 0
fit_sizenoloc = 1 if args.mode == "fit_sizenoloc" else 0
getchoosercoords = 1 if args.mode == "getchoosercoords" else 0
fixtochosenfreq = 1 if args.mode == "fixtochosenfreq" else 0

chooser = bool(args.chooser)
choose_science = args.choose_science
printimages = bool(args.printimages)
writesummary = bool(args.writesummary)

telechooser = bool(args.telechooser)
telechoice = args.telechoice










# In[ ]:


#Start to use CASA

import os

import os
import numpy as np

import shutil
import os
import math
import numpy as np
from pathlib import Path
import sys
import shutil


from casatasks import listobs
from casatasks import imhead
from casatasks import imstat
from casatasks import exportfits
from casatasks import imval
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

import analysisUtils as au
from casatools import measures
au.metool = measures

from casatools import quanta
au.qatool = quanta


#libraries
from astroquery.alma import Alma
alma = Alma()
from astropy import units as u
from astropy import coordinates
import numpy as np
from astroquery.nvas import Nvas
import requests
import shutil
import os

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel

import warnings
from astropy.utils.exceptions import AstropyWarning

from collections import Counter

warnings.simplefilter('ignore', category=AstropyWarning)

#phone text Function
import requests
def textme(words,choice):
    if choice=='go':
        bot_token = "enter bot token"
        chat_id = "enter chat id"

    url=f'https://api.telegram.org/bot{bot_token}/sendMessage?chat_id={chat_id}&text={words}'

    requests.get(url)

#function: getimview image
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from casatasks import exportfits
import os

import math
import numpy as np

if args.use_selected_clean_manifest:
    selected_clean_image_lookup, selected_clean_original_sources = load_selected_clean_manifest(
        selected_clean_manifest,
        datestr,
        flagnum,
    )
else:
    selected_clean_image_lookup = {}
    selected_clean_original_sources = []

def get_selected_sma_image_path(source):
    return selected_image_path(selected_clean_image_lookup, source)

def get_sma_image_files_for_source(source):
    selected_image_path = get_selected_sma_image_path(source)
    if selected_image_path is not None:
        return [selected_image_path]

    matched_files = []
    for smapath in smapaths:
        if not os.path.isdir(smapath):
            continue

        for ifile in os.listdir(smapath):
            if source.lower() in ifile.lower() and 'image' in ifile:
                matched_files.append(f'{smapath}/{ifile}')

    return list(dict.fromkeys(matched_files))

def is_forced_sma_niter0_upper_limit(image_path, tele, date_code=""):
    if tele != "SMA" or image_path is None:
        return False

    if str(date_code).strip() not in {"220131", "230316"}:
        return False

    normalized_path = str(image_path).replace("\\", "/")
    return "/nitercal_0/" in normalized_path

def get_sma_source_names():
    normalized_sources = []
    original_sources = []
    seen_original_sources = set()

    for original_source in selected_clean_original_sources:
        if original_source in seen_original_sources:
            continue
        original_sources.append(original_source)
        normalized_sources.append(normalize_source_name(original_source))
        seen_original_sources.add(original_source)

    for smapath in smapaths:
        if not os.path.isdir(smapath):
            continue

        for filename in os.listdir(smapath):
            if 'image' not in filename:
                continue

            original_source = filename.split('_test.image')[0]
            if original_source in seen_original_sources:
                continue

            original_sources.append(original_source)
            normalized_sources.append(normalize_source_name(original_source))
            seen_original_sources.add(original_source)

    return normalized_sources, original_sources


def remove_unselected_chosen_frequency_source_folders(chosenfreq_dir):
    """Remove each unselected source folder at most once after anchor selection."""

    if not os.path.isdir(chosenfreq_dir):
        return []

    removed = []
    for entry in os.listdir(chosenfreq_dir):
        source_dir = os.path.join(chosenfreq_dir, entry)
        if not os.path.isdir(source_dir):
            continue
        files = os.listdir(source_dir)
        if any(filename.lower().endswith('.pdf') for filename in files):
            shutil.rmtree(source_dir)
            removed.append(source_dir)
    return removed


def geom_mean_radius_arcsec(major_arcsec, minor_arcsec, is_limit=False):
    if is_limit:
        return np.nan
    if major_arcsec is None or minor_arcsec is None:
        return np.nan
    if not np.isfinite(major_arcsec) or not np.isfinite(minor_arcsec):
        return np.nan
    if major_arcsec <= 0 or minor_arcsec <= 0:
        return np.nan
    return math.sqrt(major_arcsec * minor_arcsec)

def append_separation_and_radii(
            outfile,
            source,
            instrument,
            image_path,
            separation_deg,
            chosen_major_arcsec,
            chosen_minor_arcsec,
            interest_major_arcsec,
            interest_minor_arcsec,
            chosen_is_limit=False,
            interest_is_limit=False
      ):

            chosen_geom_radius_arcsec = geom_mean_radius_arcsec(
                        chosen_major_arcsec,
                        chosen_minor_arcsec,
                        is_limit=chosen_is_limit
            )

            interest_geom_radius_arcsec = geom_mean_radius_arcsec(
                        interest_major_arcsec,
                        interest_minor_arcsec,
                        is_limit=interest_is_limit
            )

            outdir = os.path.dirname(outfile)
            if outdir and not os.path.exists(outdir):
                        os.makedirs(outdir, exist_ok=True)

            with open(outfile, "a") as f:
                        f.write(
                                    f"{source}\t{instrument}\t{image_path}\t"
                                    f"{separation_deg:.16g}\t"
                                    f"{chosen_geom_radius_arcsec:.16g}\t"
                                    f"{interest_geom_radius_arcsec:.16g}\n"
                        )

def getimview(imviewimage):
    output_file = 'test.fits'
    if os.path.exists(output_file):
        os.remove(output_file)

    exportfits(imagename=imviewimage, fitsimage=output_file, overwrite=True)

    hdul = fits.open(output_file)
    wcs = WCS(hdul[0].header, naxis=2)
    data = hdul[0].data
    data = data[0, 0, :, :] if data.ndim > 2 else data

    plt.figure(figsize=(8, 6))
    ax = plt.subplot(111, projection=wcs)
    im = ax.imshow(data, origin='lower', cmap='viridis', interpolation='none')
    plt.colorbar(im, ax=ax, label='Intensity')
    ax.set_title(f'{imviewimage}')
    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')
    ax.coords.grid(True, color='white', ls='dotted')
    ax.coords[0].set_format_unit('hour')
    ax.coords[1].set_format_unit('deg')

    hdul.close()


# In[ ]:


# In[ ]:


#function: get all files for a source
def is_archive_image_file(path):
    return os.path.isfile(path) and path.lower().endswith(('.fits', '.imfits'))


def get_archive_image_files(tele, source):
    tele_root = tele.upper()
    tele_paths = os.listdir(tele_root) if os.path.isdir(tele_root) else []
    image_files = []

    for tele_source in tele_paths:
        if source not in tele_source:
            continue
        source_dir = f'{tele_root}/{tele_source}'
        for item in os.listdir(source_dir):
            item_path = f'{source_dir}/{item}'
            if is_archive_image_file(item_path):
                image_files.append(item_path)

    return image_files


def allfiles(source):
    finfiles=[]
    finfiles.append(get_archive_image_files('VLA', source))
    finfiles.append(get_archive_image_files('ALMA', source))

    for image_path in get_sma_image_files_for_source(source):
        finfiles.append([image_path])

    return(finfiles)


# In[ ]:


#function: get multitude of coordinates
#since we are finding the pixel in a fit file from the same image we dont need to think about coordinate conversion
def get_multcoords(imfitloc_gauss2, tele,debugger):
    pixcoord1=[]
    pixcoord2=[]

    ia.close()
    ia.open(imfitloc_gauss2)
    stats=ia.statistics()

    # Extract the maximum pixel value
    max_val = stats['max'][0]

    # Get the pixel coordinates of the maximum
    max_pos = stats['maxpos']

    pix1 = stats['maxpos'][0]
    pix2 = stats['maxpos'][1]



    return(pix1,pix2,max_val)



'''
for imultgauss in multgauss:
    if source in imultgauss[0]:
        j=-1
        for iimultgauss in imultgauss:
            j=j+1
            if j==0:
                continue

            #offset_divisor=iimultgauss[0]
            iimultgauss=iimultgauss[0]
            icoords=[iimultgauss[1], iimultgauss[2]]
            icoord1=icoords[0].split(':')
            selectcoord1=f'{icoord1[0]}h{icoord1[1]}m{icoord1[2]}s'
            selectcoord2=icoords[1].replace(':','.')

            #aicoords=[]
            #aicoords=au.J2000ToICRS(f'{icoords[0]} {icoords[1]}',verbose=False).split(', ')
            #aicoord1=aicoords[0].split(':')
            #aselectcoord1=f'{aicoord1[0]}h{aicoord1[1]}m{aicoord1[2]}s'                
            #aselectcoord2=aicoords[1].replace(':','.')

            ia.close()
            ia.open(file)
            csys33 = ia.coordsys()

            #for VLA and SMA using original imagefile
            ia.close()
            ia.open(file)
            csys33 = ia.coordsys()
            pixcoords33=csys33.topixel([selectcoord1,selectcoord2])
            pixcoord1.append(pixcoords33['numeric'][0])
            pixcoord2.append(pixcoords33['numeric'][1])  


return(pixcoord1,pixcoord2)
'''


# In[ ]:


#function: get highest resolution file
#allfile=[source,allfiles(source)]
#allfiles(source) is a 3D array, each element for a different telescope

#getbestrez(allfile,icoords,6)])
def sortfiles(getfiles,alltele,tele):
    i=-1
    j=-1
    gotit=0
    for telescopes in getfiles:
        i=i+1
        #iterate past the sourcename
        if i==0:
            continue
        k=-1
        for telescope in telescopes:
            k=k+1
            for file in telescope:
                if file==None:
                    continue
                if len(file)==0:
                    continue
                if tele!='NA':
                    bestrezfile=file
                    besttele=tele
                    return(bestrezfile,besttele)
                if 'fitfiles' in file:
                    if tele=='NA':
                        continue   
                breaker=0
                #mask files need to be filtered back... need to fix earlier in process.
                if 'mask' in file:
                    os.remove(file)
                    breaker=1
                if breaker==1:
                    continue
                j=j+1
                majunit=imhead(imagename=file)['restoringbeam']['major']['unit']
                majval=imhead(imagename=file)['restoringbeam']['major']['value']
                minunit=imhead(imagename=file)['restoringbeam']['minor']['unit']
                minval=imhead(imagename=file)['restoringbeam']['minor']['value']

                if majunit!='arcsec' or minunit!='arcsec':
                    if majunit=='deg':
                        majval=majval*60*60
                        majunit='arcsec'
                    if minunit=='deg':
                        minval=minval*60*60
                        minunit='arcsec'
                    else:
                        print(f'Beam Axis does not use arcsec or deg for units.  check file:\n{file}')

                #imhead appear to give diameter values (FWHM)
                beamrad=(majval)/2
                gotit=1
                if j==0:
                    bestrezfile=file
                    bestrez=beamrad
                    besttele=alltele[k]
                if j!=0:
                    if beamrad<bestrez:
                        bestrezfile=file
                        bestrez=beamrad  
                        besttele=alltele[k]
                i=i+1  
    if gotit==0:
        return(None,None)
    if gotit==1:
        return(bestrezfile,besttele)
    #send back coordinates based on the proper coordinate system


def image_celestial_frame(image_path):
    """Read the direction reference code without changing the shared image tool."""

    frame_image = image()
    frame_coordsys = None
    try:
        if not frame_image.open(str(image_path)):
            raise RuntimeError(f'Could not open image to read its celestial frame: {image_path}')
        frame_coordsys = frame_image.coordsys()
        reference_codes = frame_coordsys.referencecode(type='direction')
        return normalize_celestial_reference_code(reference_codes)
    finally:
        if frame_coordsys is not None:
            try:
                frame_coordsys.done()
            except Exception:
                pass
        try:
            frame_image.close()
        except Exception:
            pass


def _casa_ra(raw_ra):
    value = str(raw_ra).strip()
    if 'h' in value:
        return value
    parts = value.split(':')
    if len(parts) != 3:
        raise ValueError(f'Cannot format right ascension {raw_ra!r}')
    return f'{parts[0]}h{parts[1]}m{parts[2]}s'


def coherent_coordinate_pairs(native_ra, native_dec, native_frame):
    """Return coherent ``(ICRS, J2000)`` CASA coordinate pairs."""

    frame = normalize_celestial_reference_code(native_frame)
    native_ra = str(native_ra).strip()
    native_dec_colon = colon_sexagesimal_declination(native_dec)
    if frame == ICRS:
        icrs_raw = (native_ra, native_dec_colon)
        j2000_raw = tuple(
            au.ICRSToJ2000(f'{native_ra} {native_dec_colon}', verbose=False).split(', ')
        )
    else:
        j2000_raw = (native_ra, native_dec_colon)
        icrs_raw = tuple(
            au.J2000ToICRS(f'{native_ra} {native_dec_colon}', verbose=False).split(', ')
        )
    if len(icrs_raw) != 2 or len(j2000_raw) != 2:
        raise ValueError(
            f'Coordinate conversion did not return RA and Dec for {native_ra} {native_dec}'
        )
    icrs = [_casa_ra(icrs_raw[0]), casa_sexagesimal_declination(icrs_raw[1])]
    j2000 = [_casa_ra(j2000_raw[0]), casa_sexagesimal_declination(j2000_raw[1])]
    return icrs, j2000


def coordinates_for_image(image_path, icrs_coords, j2000_coords):
    """Select a coherent coordinate pair from the image's own WCS metadata."""

    frame = image_celestial_frame(image_path)
    return list(
        coordinate_pair_for_reference_code(
            frame,
            icrs=icrs_coords,
            j2000=j2000_coords,
        )
    )


def getbestrez(getfiles,icoords,tele,alltele,choosefile,choicefile):
    attempt_limit = coordinate_candidate_attempt_limit(
        getfiles,
        choosefile=bool(choosefile),
    )
    for _attempt in range(attempt_limit):
        if choosefile==0:
            bestrezfile,besttele=sortfiles(getfiles,alltele,tele)
        else:
            bestrezfile=choicefile
            besttele=tele
        if bestrezfile==None:
            return(None)
        if ':' in icoords[0]:
            icoord1=icoords[0].split(':')
            selectcoord1=f'{icoord1[0]}h{icoord1[1]}m{icoord1[2]}s' 
        else:
            selectcoord1=icoords[0]  
        if ':' in icoords[1]:
            selectcoord2=icoords[1].replace(':','.')
        else:
            selectcoord2=icoords[1]
        aicoords=[]
        aicoords=au.J2000ToICRS(f'{icoords[0]} {icoords[1]}',verbose=False).split(', ')
        aicoord1=aicoords[0].split(':')
        aselectcoord1=f'{aicoord1[0]}h{aicoord1[1]}m{aicoord1[2]}s'
        aselectcoord2=aicoords[1].replace(':','.')
        region=f'centerbox[[{selectcoord1},{selectcoord2}],3arcsec,3arcsec]]'
        aregion=f'centerbox[[{aselectcoord1},{aselectcoord2}],3arcsec,3arcsec]]'
        try:
            best_frame = image_celestial_frame(bestrezfile)
            search_region = (
                aregion
                if best_frame == ICRS
                else region
            )
            if tele == 'NA':
                maxposf = imstat(bestrezfile, region=search_region)['maxposf']
            else:
                maxposf = imstat(bestrezfile)['maxposf']
            native_parts = [part.strip() for part in str(maxposf).split(',')]
            if len(native_parts) < 2:
                raise ValueError(f'Unexpected maxposf for {bestrezfile}: {maxposf!r}')
            icrs_pair, j2000_pair = coherent_coordinate_pairs(
                native_parts[0],
                native_parts[1],
                best_frame,
            )
            aselectcoord1, aselectcoord2 = icrs_pair
            selectcoord1, selectcoord2 = j2000_pair
            # Keep each returned coordinate pair in one frame.  The
            # historical second pair mixed J2000 right ascension with ICRS
            # declination.
            return (
                [aselectcoord1, aselectcoord2],
                [selectcoord1, selectcoord2],
                bestrezfile,
            )
        except Exception as e:
            discard_failed_coordinate_candidate_or_raise(
                getfiles,
                bestrezfile,
                choosefile=bool(choosefile),
                error=e,
            )

    return None




# In[ ]:


#function: make nicknames for the sources
def nickname(ivsource):
    if len(ivsource.split('MRK'))>1:
        ivsource=ivsource.split('MRK')[1]
    if len(ivsource.split('NGC'))>1:
        ivsource=ivsource.split('NGC')[1]
    if len(ivsource.split('J0'))>1:
        ivsource=ivsource.split('J0')[1]
    if len(ivsource.split('UGC'))>1:
        ivsource=ivsource.split('UGC')[1]
    if len(ivsource.split('ESO'))>1:
        ivsource=ivsource.split('ESO')[1]
    if len(ivsource.split('IC'))>1:
        ivsource=ivsource.split('IC')[1]
    if len(ivsource.split('B'))>1:
        ivsource=ivsource.split('B')[0]
    if len(ivsource.split('b'))>1:
        ivsource=ivsource.split('b')[0]
    if len(ivsource.split('-G009'))>1:
        ivsource=ivsource.split('-G009')[0]
    return(ivsource)


# In[ ]:


#function: print the first line
def firstline():
    sta='target'
    stele='Telescope'
    sfl='flux(mJy)'
    ssma='snr'
    sfr='frequency(GHz)'
    sda='date (index)'
    sar='beam size (")'
    sfilename='file name'
    print(f"{sta.ljust(15)}{stele.ljust(15)}{sfr.ljust(19)}{sfl.ljust(16)}{sda.ljust(15)}{sar.ljust(19)}{ssma.ljust(15)}{sfilename}")


fluxerrfile_handle = None
final_noise_audit_path = None
final_noise_audit_records = {}
observed_boundary_rejection_products = set()

FINAL_NOISE_AUDIT_FIELDS = GENERIC_FINAL_NOISE_AUDIT_FIELDS


def fluxerr_firstline(doublegaussfit):
    sta='target'
    scomponent='component'
    stele='Telescope'
    sfr='frequency(GHz)'
    sda='date (index)'
    serr='flux err(mJy)'
    sfilename='file name'
    if doublegaussfit:
        print(f"{sta.ljust(15)}{scomponent.ljust(10)}{stele.ljust(15)}{sfr.ljust(19)}{sda.ljust(15)}{serr.ljust(16)}{sfilename}", file=fluxerrfile_handle)
    else:
        print(f"{sta.ljust(15)}{stele.ljust(15)}{sfr.ljust(19)}{sda.ljust(15)}{serr.ljust(16)}{sfilename}", file=fluxerrfile_handle)


def write_fluxerr_row(source, tele, freq, date, fluxerr_mjy, filename, component=None):
    if fluxerrfile_handle is None:
        return

    sfreq=f'{freq:.2f}'
    sdate=str(date)
    if fluxerr_mjy is None:
        sfluxerr='NA'
    else:
        try:
            if math.isnan(fluxerr_mjy):
                sfluxerr='NA'
            else:
                sfluxerr=format_number(fluxerr_mjy)
        except TypeError:
            sfluxerr=format_number(fluxerr_mjy)

    if component is None:
        print(f"{source.ljust(15)}{tele.ljust(15)}{sfreq.ljust(19)}{sdate.ljust(15)}{sfluxerr.ljust(16)}{filename}", file=fluxerrfile_handle)
    else:
        scomponent=str(component)
        print(f"{source.ljust(15)}{scomponent.ljust(10)}{tele.ljust(15)}{sfreq.ljust(19)}{sdate.ljust(15)}{sfluxerr.ljust(16)}{filename}", file=fluxerrfile_handle)


# In[ ]:


#function: write image estimates
def write_region(iregion,pixcoord1,pixcoord2,apixcoord1,apixcoord2,tele,file,region,aregion,majval,minval,rotval,rotunit,write,ofixer,doublegaussfit,second_gauss_pix, debugger, fixtochosenfreq, imagefound, celestial_frame):
    worked1=0
    with open(iregion, 'w') as file2:
        if write==1:
            original_stdout2 = sys.stdout
            sys.stdout = file2

        if celestial_frame == J2000:
            try:
                if fixtochosenfreq==0 or imagefound==0:
                    imstatter= imstat(imagename=file,region=region)
                    pixcoord1=imstatter['maxpos'][0]
                    pixcoord2=imstatter['maxpos'][1]
                else:
                    pixcoord1=pixcoord1[0]
                    pixcoord2=pixcoord2[0]
            except Exception as e:
                if write==1:
                    sys.stdout = original_stdout2
                return('bad',file)
            if ofixer==1:
                #fixer=',xyabp'
                fixer=',xyabp'
            if ofixer==0:
                fixer=',abp'
            if ofixer=='nofix':
                fixer=''

            if ofixer==0:
                try:
                    print(f'{imstat(imagename=file,region=region)["max"][0]},{pixcoord1},{pixcoord2},{majval}arcsec,{minval}arcsec,{rotval}{rotunit}{fixer}') 
                except Exception as e:
                    if write==1:
                        sys.stdout = original_stdout2
                    return('bad',file)

            if ofixer==1:
                try:
                    x = round(pixcoord1)
                    y = round(pixcoord2)
                    pix = f"{x},{y}"  # properly formatted box string
                    val = imval(imagename=file, box=pix)
                    val=val['data'][0]
                    print(f'{val},{pixcoord1},{pixcoord2},{majval}arcsec,{minval}arcsec,{rotval}{rotunit}{fixer}') 

                except Exception as e:
                    if write==1:
                        sys.stdout = original_stdout2
                    return('bad',file)

            if ofixer=='nofix':
                try:
                    print(f'{imstat(imagename=file,region=region)["max"][0]},{pixcoord1},{pixcoord2},{majval}arcsec,{minval}arcsec,{rotval}{rotunit}{fixer}') 
                except Exception as e:
                    if write==1:
                        sys.stdout = original_stdout2
                    return('bad',file)

        if celestial_frame == ICRS:
            try:
                if fixtochosenfreq==0 or imagefound==0:
                    imstatter= imstat(imagename=file,region=aregion)
                    apixcoord1=imstatter['maxpos'][0]
                    apixcoord2=imstatter['maxpos'][1]
                else:
                    apixcoord1=apixcoord1[0]
                    apixcoord2=apixcoord2[0]
            except Exception as e:
                if write==1:
                    sys.stdout = original_stdout2
                return('bad',file)
            if ofixer==1:
                #fixer=',xyabp'
                fixer=',xyabp'
            if ofixer==0:
                fixer=',abp'
            elif ofixer=='nofix':
                fixer=''

            if ofixer==0:
                try:
                    print(f'{imstat(imagename=file,region=aregion)["max"][0]},{apixcoord1},{apixcoord2},{majval}arcsec,{minval}arcsec,{rotval}{rotunit}{fixer}')
                except Exception as e:
                    if write==1:
                        sys.stdout = original_stdout2
                    return('bad',file)

            if ofixer==1:
                try:
                    x = round(apixcoord1)
                    y = round(apixcoord2)
                    pix = f"{x},{y}"  # properly formatted box string
                    val = imval(imagename=file, box=pix)
                    val=val['data'][0]
                    print(f'{val},{apixcoord1},{apixcoord2},{majval}arcsec,{minval}arcsec,{rotval}{rotunit}{fixer}')
                except Exception as e:
                    if write==1:
                        sys.stdout = original_stdout2
                    return('bad',file)

            if ofixer=='nofix':
                try:
                    print(f'{imstat(imagename=file,region=aregion)["max"][0]},{apixcoord1},{apixcoord2},{majval}arcsec,{minval}arcsec,{rotval}{rotunit}{fixer}') 
                except Exception as e:
                    if write==1:
                        sys.stdout = original_stdout2
                    return('bad',file)

        if doublegaussfit:
            doublefixer=''
            try:
                print(f'{second_gauss_pix[2]},{second_gauss_pix[0]},{second_gauss_pix[1]},{majval}arcsec,{minval}arcsec,{rotval}{rotunit}{doublefixer}') 
                worked1=1
            except Exception as e:
                if write==1:
                    sys.stdout = original_stdout2
                return('bad',file)
            if debugger:
                print('3.1.2')

        if write==1:
            sys.stdout = original_stdout2
    #worked1 shows that the file has no imstat under the proper region
    #worked2 shows that the regrid file has no imstat under the proper region
    return([worked1])


# In[ ]:


#function: makegauss image
from casatools import componentlist
cl = componentlist()
from casatools import quanta
qa = quanta()

def makegauss(coords,freq,frequn,file,imfitter,tele,majorax,minorax,posang):

    ppas1=imfitter['pixelsperarcsec'][0]
    ppas2=imfitter['pixelsperarcsec'][1]

    aspp1=1/ppas1
    aspp2=1/ppas2

    #divide the grid into smaller pieces
    if aspp1<aspp2:
        aspp=aspp1
    else:
        aspp=aspp2

    fulpix1=imstat(file)['trc'][0]
    fulpix2=imstat(file)['trc'][1]

    if os.path.exists('testg.fits'):
        os.remove('testg.fits')
    if os.path.exists('testg.im'):
        shutil.rmtree('testg.im')
        os.mkdir('testg.im')

    ra=coords[0]
    dec=coords[1]
    celestial_frame = image_celestial_frame(file)
    direction = f"{celestial_frame} {ra} {dec}"
    cl.done()

    cl.addcomponent(dir=direction, flux=1.0, fluxunit='Jy', freq=f'{freq}{frequn}', shape="Gaussian", majoraxis=majorax, minoraxis=minorax, positionangle=posang)

    aspp=f'{str(aspp)}arcsec'
    if os.path.exists("testg.im"):
        shutil.rmtree("testg.im")
    ia.fromshape("testg.im",[fulpix1,fulpix2],overwrite=True)
    cs=ia.coordsys()
    cs.setunits(['rad','rad'])
    cs.setreferencecode(celestial_frame, type='direction', adjust=False)
    cell_rad=qa.convert(qa.quantity(aspp),"rad")['value']

    cs.setincrement([-cell_rad,cell_rad],'direction')
    cs.setreferencevalue([qa.convert(ra,'rad')['value'],qa.convert(dec,'rad')['value']],type="direction")
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.modify(cl.torecord(),subtract=False)
    exportfits(imagename='testg.im',fitsimage='testg.fits',overwrite=True)


# In[ ]:


#function: makedot image
from casatools import componentlist
cl = componentlist()
from casatools import quanta
qa = quanta()

def makedot(coords,freq,frequn,file,imfitter,tele,type):

    ppas1=imfitter['pixelsperarcsec'][0]
    ppas2=imfitter['pixelsperarcsec'][1]

    aspp1=1/ppas1
    aspp2=1/ppas2

    #divide the grid into smaller pieces
    if aspp1<aspp2:
        aspp=aspp1
    else:
        aspp=aspp2

    fulpix1=imstat(file)['trc'][0]
    fulpix2=imstat(file)['trc'][1]

    if type=='orig':
        pointname='testp'
    if type=='fit':
        pointname='fitp'
    if os.path.exists(f'{pointname}.fits'):
        os.remove(f'{pointname}.fits')
    if os.path.exists(f'{pointname}.im'):
        shutil.rmtree(f'{pointname}.im')
        os.mkdir(f'{pointname}.im')
    ra=coords[0]
    dec=coords[1]
    celestial_frame = image_celestial_frame(file)
    direction = f"{celestial_frame} {ra} {dec}"
    cl.done()

    cl.addcomponent(dir=direction, flux=1.0, fluxunit='Jy', freq=f'{freq}{frequn}', shape="Point")
    aspp=f'{str(aspp)}arcsec'
    if os.path.exists(f"{pointname}.im"):
        shutil.rmtree(f"{pointname}.im")
        os.mkdir(f"{pointname}.im")
    ia.fromshape(f"{pointname}.im",[fulpix1,fulpix2],overwrite=True)
    cs=ia.coordsys()
    cs.setunits(['rad','rad'])
    cs.setreferencecode(celestial_frame, type='direction', adjust=False)
    cell_rad=qa.convert(qa.quantity(aspp),"rad")['value']

    cs.setincrement([-cell_rad,cell_rad],'direction')
    cs.setreferencevalue([qa.convert(ra,'rad')['value'],qa.convert(dec,'rad')['value']],type="direction")
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.modify(cl.torecord(),subtract=False)
    exportfits(imagename=f'{pointname}.im',fitsimage=f'{pointname}.fits',overwrite=True)


# In[ ]:


#funcion: make dot with contour image
import os
from casatools import image
from casatasks import exportfits
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np
from reproject import reproject_interp
from matplotlib.lines import Line2D

def makedotcontour(file, imfitloc, peak_value, beamsize, newfilename, printimages,color,pixlen,oicoords,freq2,frequnit,imfitter,alltele,tele, getchoosercoords, fixtochosenfreq, imagefound, prelimchecker, fit_sizenoloc, doublegaussfit):
    sepp=newfilename.split('.')
    source=sepp[0]
    freq=sepp[1]
    year=sepp[2]
    month=sepp[3]
    daynum=sepp[4].split('_')
    day=daynum[0]
    if isinstance(color, str):
        colors='one'
        bothcolorteller=False
    elif isinstance(color, list):
        bothcolorteller=True
        bothb=1
        bothr=1
        #first one is the selected fit... it may be red
        firstcolor=color[0]
        for icolor in color:
            if icolor=='r':
                bothb=0
            if icolor=='b':
                bothr=0
        if bothb==1:
            color='b'
        elif bothr==1:
            input('both are red')
            return('error')
        elif firstcolor=='1':
            color='b'
            if firstcolor=='b':
                input('the fixed location is a detection and the freely varying is a nondetection')
                return('error')

    # Define the paths
    wcs_main = None  # Initialize wcs_main to ensure it's accessible in the except block

    imagename = file
    contour_fits = "testp.fits"

    nicename=f'{source} {freq} {year}/{month}/{day}'

    # Zoom region size (in pixels, e.g., a square region of size `zoom_width x zoom_width`)
    zoom_width = pixlen*3

    # Temporary output FITS files
    output_file = "test.fits"
    imfit_fits = "imfit.fits"

    # Ensure previous files are cleared
    for f in [output_file, imfit_fits]:
        if os.path.exists(f):
            os.remove(f)

    try:
        # Export the main image and imfitloc as FITS files
        exportfits(imagename=imagename, fitsimage=output_file, overwrite=True)
        exportfits(imagename=imfitloc, fitsimage=imfit_fits, overwrite=True)          

        # Load main FITS file
        hdul_main = fits.open(output_file)
        wcs_main = WCS(hdul_main[0].header, naxis=2)
        data_main = hdul_main[0].data
        data_main = data_main[0, 0, :, :] if data_main.ndim > 2 else data_main
        data_main = data_main * brightness_to_mjy_per_beam_factor(
            hdul_main[0].header.get('BUNIT', '')
        )
        img_shape = data_main.shape  # Shape of the main image (pixels)

        # Load imfitloc FITS file
        hdul_imfit = fits.open(imfit_fits)
        wcs_imfit = WCS(hdul_imfit[0].header, naxis=2)
        data_imfit = hdul_imfit[0].data
        data_imfit = data_imfit[0, 0, :, :] if data_imfit.ndim > 2 else data_imfit


        # Replace NaNs in data_imfit
        data_imfit = np.nan_to_num(data_imfit, nan=0.0)


        # Reproject imfit data onto main image WCS
        reprojected_data, _ = reproject_interp(
            (data_imfit, wcs_imfit),
            wcs_main,
            shape_out=img_shape
        )          

        #Load contour FITS file if it exists
        if os.path.exists(contour_fits):
            hdul_contour = fits.open(contour_fits)
            wcs_contour = WCS(hdul_contour[0].header, naxis=2)
            dot_ra = wcs_contour.wcs.crval[0]
            dot_dec = wcs_contour.wcs.crval[1]
            dot_pixel = wcs_main.world_to_pixel_values(dot_ra, dot_dec)
        else:
            dot_pixel = None

        dot_pixel2 = None 
        if fixtochosenfreq==1 and imagefound==1:
            choosefile=1
            choicefile=imfitloc
            gausloc=getbestrez([source,[[imfitloc],[],[]]],oicoords,tele,alltele,choosefile,choicefile)
            #the gaussian location should be that of the imfitloc
            gauslocfinder=gausloc[2].replace('/fitfiles','')
            gauslocfinder=gauslocfinder.replace('.fits.fit','.fits')
            gauslocfinder_file=file.replace('.imfits','.fits')
            if gauslocfinder_file!=gauslocfinder:
                pass
                #textme(f"ERROR!",'go')
                #textme(f"{gauslocfinder}",'go')
                #textme(f"{gauslocfinder_file}",'go')
                #input('check')
            if gausloc=='None':
                return('None')
            if gausloc==None:
                return(None)
            type='fit'
            acoords=gausloc[0]
            coords=gausloc[1]
            plot_coords = coordinates_for_image(file, acoords, coords)
            makedot(plot_coords,freq2,frequnit,file,imfitter,tele,type)

            contour_fits2="fitp.fits" 
            if os.path.exists(contour_fits2):
                hdul_contour2 = fits.open(contour_fits2)
                wcs_contour2 = WCS(hdul_contour2[0].header, naxis=2)
                dot_ra2 = wcs_contour2.wcs.crval[0]
                dot_dec2 = wcs_contour2.wcs.crval[1]
                dot_pixel2 = wcs_main.world_to_pixel_values(dot_ra2, dot_dec2)
            else:
                dot_pixel2 = None
        if color=='b':
            choosefile=1
            choicefile=imfitloc
            gausloc=getbestrez([source,[[imfitloc],[],[]]],oicoords,tele,alltele,choosefile,choicefile)
            #the gaussian location should be that of the imfitloc
            gauslocfinder=gausloc[2].replace('/fitfiles','')
            gauslocfinder=gauslocfinder.replace('.fits.fit','.fits')
            gauslocfinder_file=file.replace('.imfits','.fits')
            if gauslocfinder_file!=gauslocfinder:
                pass
                #textme(f"ERROR!",'go')
                #textme(f"{gauslocfinder}",'go')
                #textme(f"{gauslocfinder_file}",'go')
                #input('check')
            if gausloc=='None':
                return('None')
            if gausloc==None:
                return(None)
            type='fit'
            acoords=gausloc[0]
            coords=gausloc[1]
            plot_coords = coordinates_for_image(file, acoords, coords)
            makedot(plot_coords,freq2,frequnit,file,imfitter,tele,type)

            contour_fits2="fitp.fits" 
            if os.path.exists(contour_fits2):
                hdul_contour2 = fits.open(contour_fits2)
                wcs_contour2 = WCS(hdul_contour2[0].header, naxis=2)
                dot_ra2 = wcs_contour2.wcs.crval[0]
                dot_dec2 = wcs_contour2.wcs.crval[1]
                dot_pixel2 = wcs_main.world_to_pixel_values(dot_ra2, dot_dec2)
            else:
                dot_pixel2 = None

        # Plot the main image
        plt.figure(figsize=(8, 6))
        ax = plt.subplot(111, projection=wcs_main)
        im = ax.imshow(data_main, origin='lower', cmap='viridis', interpolation='none')
        # Adjust colorbar
        cbar = plt.colorbar(im, ax=ax, label='mJy/beam')

        # Define contour levels dynamically based on reprojected data
        vmin = np.nanmin(reprojected_data)
        vmax = np.nanmax(reprojected_data)

        # quick sanity: any finite pixels?
        finite_count = np.isfinite(reprojected_data).sum()

        if not np.isfinite(vmax) and not np.isfinite(vmin):
            print("No finite pixels at all. Skipping contour.")
            levels_imfit = []
        else:
            if vmax > 0:
                peak = vmax
            else:
                peak = vmin

            levels_imfit = [0.5 * peak]  # FWHM (always positive)

        #fwhm_level2 = 0.5 * np.max(reprojected_data)  # Full Width at Half Maximum
        #levels_imfit = [fwhm_level2]

        #fwhm_level2_1 = 0.01*np.nanmax(reprojected_data)  # Use nanmax to ignore NaNs
        #fwhm_level2_2 = 0.5*np.nanmax(reprojected_data)  # Use nanmax to ignore NaNs
        #fwhm_level2_3 = 0.9*np.nanmax(reprojected_data)  # Use nanmax to ignore NaNs
        #levels_imfit = [fwhm_level2_1,fwhm_level2_2,fwhm_level2_3] 

        # Extract and analyze the contours
        contours = ax.contour(reprojected_data, levels=levels_imfit, colors=color, origin="lower")

        #Define the image boundaries explicitly
        x_min, x_max = 0, img_shape[1] - 1
        y_min, y_max = 0, img_shape[0] - 1

        # Debug: Print boundary details

        # Threshold for boundary classification
        threshold = 90  # Percentage threshold


        # Check and classify vertices along the image edges
        for path in contours.get_paths():

            vertices = path.vertices  # Contour vertices (x, y)
            x_coords, y_coords = vertices[:, 0], vertices[:, 1]


            percentages=[]
            for x, y in zip(x_coords, y_coords):
                # Calculate proximity to the nearest boundary
                proximity_x = min(abs(x - x_min), abs(x - x_max)) / (x_max - x_min)
                proximity_y = min(abs(y - y_min), abs(y - y_max)) / (y_max - y_min)
                proximity = min(proximity_x, proximity_y)
                percentage = (1 - proximity) * 100
                percentages.append(percentage)
                # Classify as "on the boundary" if percentage >= threshold
            highest=0
            for percentage in percentages:
                if percentage>highest:
                    highest=percentage
            if highest > threshold:

                #print(f"x={x:.2f}, y={y:.2f} -> Proximity to boundary: {percentage:.2f}%")
                # Close FITS files
                hdul_main.close()
                hdul_imfit.close()
                if os.path.exists(contour_fits):
                    hdul_contour.close()
                if color=='b':
                    if os.path.exists(contour_fits2):
                        hdul_contour2.close()   
                return('ONBOUNDRY')
        # Create a custom blue line for the legend
        if bothcolorteller==False:
            if '*' not in peak_value:   
                blue_line = Line2D([0], [0], color=color, lw=1.5, label=f'0.5 contour of {peak_value} mJy Gaussian Peak')
            elif peak>=0:
                blue_line = Line2D([0], [0], color=color, lw=1.5, label=f'LOW SNR 0.5 Gaussian Contour: 3X RMS: {peak_value} mJy')
            elif peak<0:
                blue_line = Line2D([0], [0], color=color, lw=1.5, linestyle="--", label=f'LOW SNR 0.5 Gaussian Contour: 3X RMS: {peak_value} mJy')

        else:
            blue_line = Line2D([0], [0], color='b', lw=1.5, label=f'0.01, 0.5, 0.9 contour mJy 2-Component Gaussian Peak \n {peak_value[0]} (free) , {peak_value[1]} (fixed)')
        handles = [blue_line]
        # Overlay red dot if available
        if dot_pixel:
            ax.plot(dot_pixel[0], dot_pixel[1], 'ro', markersize=2, label="Reference Dot")
            if fixtochosenfreq==0:
                locationofdot='Location from NED'
            if fixtochosenfreq==1 and imagefound==1:
                locationofdot='Location from chosen freq'
            if fixtochosenfreq==1 and imagefound==0:
                locationofdot='Location from NED (no chosen freq)'
            if doublegaussfit==1:
                locationofdot='Location from chosen freq'
            red_dot = Line2D([0], [0], color='red', marker='o', markersize=6, label=locationofdot, linestyle='')
            handles.append(red_dot)  
        if dot_pixel2:
            ax.plot(dot_pixel2[0], dot_pixel2[1], 'bo', markersize=2, label="Reference Dot2")
            blue_dot = Line2D([0], [0], color='blue', marker='o', markersize=6, label='Fitted Location for this image', linestyle='')
            handles.append(blue_dot)  

        black_line = Line2D([0], [0], color='black', lw=1.5, label=f'{beamsize}" FWHM beam')
        handles.append(black_line)
        # Add legend
        ax.legend(handles=handles, loc='upper right')
        # Apply zoom
        if dot_pixel:
            center_x, center_y = dot_pixel
            x_min = max(center_x - zoom_width / 2, 0)
            x_max = min(center_x + zoom_width / 2, img_shape[1])
            y_min = max(center_y - zoom_width / 2, 0)
            y_max = min(center_y + zoom_width / 2, img_shape[0])
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y_min, y_max)
        # Finalize plot
        if tele=='SMA':
            file='New Observation'
        ax.set_title(f'{file}\n{nicename}')
        ax.set_xlabel("Right Ascension")
        ax.set_ylabel("Declination")
        ax.coords.grid(True, color='white', ls='dotted')
        ax.coords[0].set_format_unit('hour')
        ax.coords[1].set_format_unit('deg')

        if prelimchecker:
            output_plot_filename=f'images/{flagnum}/prelim/{source}/{newfilename}.pdf'
            if not os.path.exists(f'images/{flagnum}/prelim/{source}'):
                os.makedirs(f'images/{flagnum}/prelim/{source}')
        elif fit_sizenoloc:
            output_plot_filename=f'images/{flagnum}/noloc/{source}/{newfilename}.pdf'
            if not os.path.exists(f'images/{flagnum}/noloc/{source}'):
                os.makedirs(f'images/{flagnum}/noloc/{source}')
        elif getchoosercoords:
            output_plot_filename=f'images/{flagnum}/chosenfreq/{source}/{newfilename}.pdf'
            if not os.path.exists(f'images/{flagnum}/chosenfreq/{source}'):
                os.makedirs(f'images/{flagnum}/chosenfreq/{source}')
        elif fixtochosenfreq==1:
            output_plot_filename=f'images/{flagnum}/fixedtochosen/{source}/{newfilename}.pdf'
            if not os.path.exists(f'images/{flagnum}/fixedtochosen/{source}'):
                os.makedirs(f'images/{flagnum}/fixedtochosen/{source}')
        elif doublegaussfit==1:
            output_plot_filename=f'images/{flagnum}/doublegauss/{source}/{newfilename}.pdf'
            if not os.path.exists(f'images/{flagnum}/doublegauss/{source}'):
                os.makedirs(f'images/{flagnum}/doublegauss/{source}')

        breaker=0
        i=0
        while breaker==0:
            i=i+1
            newfilename=output_plot_filename.replace('.pdf',f'_{i}.pdf')
            newdate=f'{year}/{month}/{day}_{i}'
            if not os.path.exists(newfilename):
                breaker = 1

        output_plot_filename=newfilename

        plt.savefig(output_plot_filename, dpi=100, bbox_inches='tight')
        if printimages==True:
            plt.show()
        plt.close()
        # Close FITS files
        hdul_main.close()
        hdul_imfit.close()
        if os.path.exists(contour_fits):
            hdul_contour.close()
        if color=='b':
            if os.path.exists(contour_fits2):
                hdul_contour2.close()       

    except Exception as e:
        print(f"Error1: {e} for {file}")

    return(newdate)


# In[ ]:


#function: make only dot
#funcion: make dot with contour image
import os
from casatools import image
from casatasks import exportfits
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np
from reproject import reproject_interp
from matplotlib.lines import Line2D

def makeonlydotcontour(file, imfitloc, noise_value, beamsize, newfilename, printimages,color,pixlen,oicoords,freq2,frequnit,imfitter,alltele,tele, getchoosercoords, fixtochosenfreq, imagefound, prelimchecker, fit_sizenoloc, doublegaussfit):
    sepp=newfilename.split('.')
    source=sepp[0]
    freq=sepp[1]
    year=sepp[2]
    month=sepp[3]
    daynum=sepp[4].split('_')
    day=daynum[0]
    if isinstance(color, str):
        colors='one'
        bothcolorteller=False
    elif isinstance(color, list):
        bothcolorteller=True
        bothb=1
        bothr=1
        #first one is the selected fit... it may be red
        firstcolor=color[0]
        for icolor in color:
            if icolor=='r':
                bothb=0
            if icolor=='b':
                bothr=0
        if bothb==1:
            color='b'
        elif bothr==1:
            input('both are red')
            return('error')
        elif firstcolor=='1':
            color='b'
            if firstcolor=='b':
                input('the fixed location is a detection and the freely varying is a nondetection')
                return('error')

    # Define the paths
    wcs_main = None  # Initialize wcs_main to ensure it's accessible in the except block

    imagename = file
    contour_fits = "testp.fits"

    nicename=f'{source} {freq} {year}/{month}/{day}'

    # Zoom region size (in pixels, e.g., a square region of size `zoom_width x zoom_width`)
    zoom_width = pixlen

    # Temporary output FITS files
    output_file = "test.fits"
    imfit_fits = "imfit.fits"

    # Ensure previous files are cleared
    for f in [output_file, imfit_fits]:
        if os.path.exists(f):
            os.remove(f)

    try:
        # Export the main image and imfitloc as FITS files
        exportfits(imagename=imagename, fitsimage=output_file, overwrite=True)
        exportfits(imagename=imfitloc, fitsimage=imfit_fits, overwrite=True)          

        # Load main FITS file
        hdul_main = fits.open(output_file)
        wcs_main = WCS(hdul_main[0].header, naxis=2)
        data_main = hdul_main[0].data
        data_main = data_main[0, 0, :, :] if data_main.ndim > 2 else data_main
        data_main = data_main * brightness_to_mjy_per_beam_factor(
            hdul_main[0].header.get('BUNIT', '')
        )
        img_shape = data_main.shape  # Shape of the main image (pixels)

        # Load imfitloc FITS file
        hdul_imfit = fits.open(imfit_fits)
        wcs_imfit = WCS(hdul_imfit[0].header, naxis=2)
        data_imfit = hdul_imfit[0].data
        data_imfit = data_imfit[0, 0, :, :] if data_imfit.ndim > 2 else data_imfit


        # Replace NaNs in data_imfit
        data_imfit = np.nan_to_num(data_imfit, nan=0.0)


        # Reproject imfit data onto main image WCS
        reprojected_data, _ = reproject_interp(
            (data_imfit, wcs_imfit),
            wcs_main,
            shape_out=img_shape
        )          

        #Load contour FITS file if it exists
        if os.path.exists(contour_fits):
            hdul_contour = fits.open(contour_fits)
            wcs_contour = WCS(hdul_contour[0].header, naxis=2)
            dot_ra = wcs_contour.wcs.crval[0]
            dot_dec = wcs_contour.wcs.crval[1]
            dot_pixel = wcs_main.world_to_pixel_values(dot_ra, dot_dec)
        else:
            dot_pixel = None


        # Plot the main image
        plt.figure(figsize=(8, 6))
        ax = plt.subplot(111, projection=wcs_main)
        im = ax.imshow(data_main, origin='lower', cmap='viridis', interpolation='none')
        # Adjust colorbar
        cbar = plt.colorbar(im, ax=ax, label='mJy/beam')

        # Define contour levels dynamically based on reprojected data
        vmin = np.nanmin(reprojected_data)
        vmax = np.nanmax(reprojected_data)

        # quick sanity: any finite pixels?
        finite_count = np.isfinite(reprojected_data).sum()

        #Define the image boundaries explicitly
        x_min, x_max = 0, img_shape[1] - 1
        y_min, y_max = 0, img_shape[0] - 1

        # Debug: Print boundary details

        # Threshold for boundary classification
        threshold = 90  # Percentage threshold


        # Overlay red dot if available
        if dot_pixel:
            ax.plot(dot_pixel[0], dot_pixel[1], 'ro', markersize=2, label="Reference Dot")
            if fixtochosenfreq==0:
                locationofdot='Location from NED'
            if fixtochosenfreq==1 and imagefound==1:
                locationofdot='Location from chosen freq'
            if fixtochosenfreq==1 and imagefound==0:
                locationofdot='Location from NED (no chosen freq)'
            if doublegaussfit==1:
                locationofdot='Location from chosen freq'
            red_dot = Line2D([0], [0], color='red', marker='o', markersize=6, label=locationofdot, linestyle='')
            handles = [red_dot] 
            red_line = Line2D([0], [0], color=color, lw=1.5, label=f'LOW SNR 0.5 Gaussian Contour: 3X RMS: {noise_value} mJy')
            handles.append(red_line)


        black_line = Line2D([0], [0], color='black', lw=1.5, label=f'{beamsize}" FWHM beam')
        handles.append(black_line)
        # Add legend
        ax.legend(handles=handles, loc='upper right')
        # Apply zoom
        if dot_pixel:
            center_x, center_y = dot_pixel
            x_min = max(center_x - zoom_width / 2, 0)
            x_max = min(center_x + zoom_width / 2, img_shape[1])
            y_min = max(center_y - zoom_width / 2, 0)
            y_max = min(center_y + zoom_width / 2, img_shape[0])
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y_min, y_max)
        # Finalize plot
        if tele=='SMA':
            file='New Observation'
        ax.set_title(f'{file}\n{nicename}')
        ax.set_xlabel("Right Ascension")
        ax.set_ylabel("Declination")
        ax.coords.grid(True, color='white', ls='dotted')
        ax.coords[0].set_format_unit('hour')
        ax.coords[1].set_format_unit('deg')

        if prelimchecker:
            output_plot_filename=f'images/{flagnum}/prelim/{source}/{newfilename}.pdf'
            if not os.path.exists(f'images/{flagnum}/prelim/{source}'):
                os.makedirs(f'images/{flagnum}/prelim/{source}')
        elif fit_sizenoloc:
            output_plot_filename=f'images/{flagnum}/noloc/{source}/{newfilename}.pdf'
            if not os.path.exists(f'images/{flagnum}/noloc/{source}'):
                os.makedirs(f'images/{flagnum}/noloc/{source}')
        elif getchoosercoords:
            output_plot_filename=f'images/{flagnum}/chosenfreq/{source}/{newfilename}.pdf'
            if not os.path.exists(f'images/{flagnum}/chosenfreq/{source}'):
                os.makedirs(f'images/{flagnum}/chosenfreq/{source}')
        elif fixtochosenfreq==1:
            output_plot_filename=f'images/{flagnum}/fixedtochosen/{source}/{newfilename}.pdf'
            if not os.path.exists(f'images/{flagnum}/fixedtochosen/{source}'):
                os.makedirs(f'images/{flagnum}/fixedtochosen/{source}')
        elif doublegaussfit==1:
            output_plot_filename=f'images/{flagnum}/doublegauss/{source}/{newfilename}.pdf'
            if not os.path.exists(f'images/{flagnum}/doublegauss/{source}'):
                os.makedirs(f'images/{flagnum}/doublegauss/{source}')

        breaker=0
        i=0
        while breaker==0:
            i=i+1
            newfilename=output_plot_filename.replace('.pdf',f'_{i}.pdf')
            newdate=f'{year}/{month}/{day}_{i}'
            if not os.path.exists(newfilename):
                breaker = 1

        output_plot_filename=newfilename

        plt.savefig(output_plot_filename, dpi=100, bbox_inches='tight')
        if printimages==True:
            plt.show()
        plt.close()
        # Close FITS files
        hdul_main.close()
        hdul_imfit.close()
        if os.path.exists(contour_fits):
            hdul_contour.close()  

    except Exception as e:
        print(f"Error1: {e} for {file}")

    return(newdate)


# In[ ]:


#format beamsize number
def format_number(number):
        return format_two_decimal(number)


# In[ ]:


# Convert image RMS values only after validating the CASA image brightness unit.
def image_rms_mjy_per_beam(imagename, region):
        brightness_unit = imhead(imagename=imagename, mode="get", hdkey="bunit")
        if isinstance(brightness_unit, dict):
                brightness_unit = brightness_unit.get("value", brightness_unit.get("unit", ""))
        rms_value = imstat(imagename=imagename, region=region)['rms'][0]
        return brightness_to_mjy_per_beam(rms_value, brightness_unit)


def _casa_stat_scalar(statistics, field):
        if field not in statistics:
                raise ValueError(f'CASA imstat result is missing {field!r}')
        values = statistics[field]
        try:
                value = values[0]
        except (IndexError, KeyError, TypeError):
                value = values
        value = float(value)
        if not math.isfinite(value):
                raise ValueError(f'CASA imstat {field!r} is not finite: {value!r}')
        return value


def image_statistics_mjy_per_beam(imagename, region):
        """Return the exact CASA classic statistics used by the final audit."""

        brightness_unit = imhead(imagename=imagename, mode="get", hdkey="bunit")
        if isinstance(brightness_unit, dict):
                brightness_unit = brightness_unit.get(
                    "value", brightness_unit.get("unit", "")
                )
        conversion = brightness_to_mjy_per_beam_factor(brightness_unit)
        statistics = imstat(imagename=imagename, region=region)
        npts = _casa_stat_scalar(statistics, "npts")
        if npts < 0.0 or not npts.is_integer():
                raise ValueError(f'CASA imstat npts is not a nonnegative integer: {npts!r}')
        return {
            "brightness_unit": str(brightness_unit),
            "npts": int(npts),
            "rms_mjy_per_beam": _casa_stat_scalar(statistics, "rms") * conversion,
            "sigma_mjy_per_beam": _casa_stat_scalar(statistics, "sigma") * conversion,
            "mean_mjy_per_beam": _casa_stat_scalar(statistics, "mean") * conversion,
            "scaled_mad_mjy_per_beam": (
                _casa_stat_scalar(statistics, "medabsdevmed")
                * conversion
                * SCALED_MAD_NORMAL_FACTOR
            ),
        }


def _blank_final_noise_audit_record():
        return {field: "" for field in FINAL_NOISE_AUDIT_FIELDS}


def _write_final_noise_audit_records():
        if final_noise_audit_path is None:
                return
        final_noise_audit_path.parent.mkdir(parents=True, exist_ok=True)
        temporary = final_noise_audit_path.with_suffix(
            final_noise_audit_path.suffix + ".tmp"
        )
        with temporary.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=FINAL_NOISE_AUDIT_FIELDS,
                    lineterminator="\n",
                )
                writer.writeheader()
                for product in sorted(final_noise_audit_records):
                        writer.writerow(final_noise_audit_records[product])
        os.replace(temporary, final_noise_audit_path)


def record_final_noise_audit(record):
        """Persist one new complete diagnostic before the hard area gate."""

        if final_noise_audit_path is None:
                return
        if set(record) != set(FINAL_NOISE_AUDIT_FIELDS):
                raise ValueError("Final-noise audit record does not match its exact schema")
        product = str(record["product"]).strip()
        if not product:
                raise ValueError("Final-noise audit product cannot be blank")
        if product in final_noise_audit_records:
                raise ValueError(f'Duplicate final-noise audit product: {product}')
        final_noise_audit_records[product] = dict(record)
        _write_final_noise_audit_records()


def update_final_noise_audit(product, **updates):
        if final_noise_audit_path is None:
                return
        product = portable_release_product(product)
        if product not in final_noise_audit_records:
                raise ValueError(f'No pending final-noise audit record for {product}')
        unknown = set(updates) - set(FINAL_NOISE_AUDIT_FIELDS)
        if unknown:
                raise ValueError(f'Unknown final-noise audit fields: {sorted(unknown)}')
        final_noise_audit_records[product].update(updates)
        _write_final_noise_audit_records()


def record_reviewed_boundary_rejection(product):
        """Finalize one boundary diagnostic only when exact policy admits it."""

        product = portable_release_product(product)
        policy = reviewed_boundary_rejections_by_product.get(product)
        if policy is None:
                raise RuntimeError(
                    f'Unreviewed final-fit boundary rejection: {product}'
                )
        if product in observed_boundary_rejection_products:
                raise RuntimeError(
                    f'Duplicate final-fit boundary rejection: {product}'
                )
        update_final_noise_audit(product, row_disposition="rejected_on_boundary")
        validate_boundary_rejection_record_identity(
                final_noise_audit_records[product], policy
        )
        observed_boundary_rejection_products.add(product)


def _direction_pixel_widths_arcsec(image_header, imagename):
        axis_units = list(image_header.get("axisunits", ()))
        increments = list(image_header.get("incr", ()))
        if len(axis_units) < 2 or len(increments) < 2:
                raise ValueError(f'Image has no two-dimensional direction increment: {imagename}')
        pixel_width = abs(angular_value_to_arcsec(increments[0], axis_units[0]))
        pixel_height = abs(angular_value_to_arcsec(increments[1], axis_units[1]))
        return pixel_width, pixel_height


def _direction_image_shape_pixels(image_header, imagename):
        shape = list(image_header.get("shape", ()))
        if len(shape) < 2:
                raise ValueError(f'Image has no two-dimensional shape: {imagename}')
        width = float(shape[0])
        height = float(shape[1])
        if (
            not math.isfinite(width)
            or not math.isfinite(height)
            or width <= 0.0
            or height <= 0.0
            or not width.is_integer()
            or not height.is_integer()
        ):
                raise ValueError(f'Image has an invalid direction shape: {imagename}')
        return int(width), int(height)


def _restoring_beam_arcsec(image_header, imagename):
        beam = image_header.get("restoringbeam")
        if not isinstance(beam, dict):
                raise ValueError(f'Image has no single restoring beam: {imagename}')
        major = angular_value_to_arcsec(
            beam["major"]["value"], beam["major"]["unit"]
        )
        minor = angular_value_to_arcsec(
            beam["minor"]["value"], beam["minor"]["unit"]
        )
        pa_unit = str(beam["positionangle"]["unit"]).strip().casefold()
        if pa_unit != "deg":
                raise ValueError(
                    f'Restoring-beam position angle is not in degrees for {imagename}'
                )
        return major, minor, float(beam["positionangle"]["value"])


def _native_fitted_direction_degrees(fitted_component, imagename):
        """Return the fitted celestial center in the image's native frame.

        CASA's component result links ``pixelcoords`` and ``shape.direction``
        to the same fitted center.  Persisting both makes the final identity
        audit possible under system Python without reopening FITS data or
        depending on Astropy.  A component/image frame mismatch fails closed.
        """

        try:
                direction = fitted_component["shape"]["direction"]
                component_frame = normalize_celestial_reference_code(
                    direction["refer"]
                )
                ra_deg = angular_value_to_arcsec(
                    direction["m0"]["value"], direction["m0"]["unit"]
                ) / 3600.0
                dec_deg = angular_value_to_arcsec(
                    direction["m1"]["value"], direction["m1"]["unit"]
                ) / 3600.0
        except (KeyError, TypeError) as error:
                raise ValueError(
                    f'Final fitted component has no complete native direction: {imagename}'
                ) from error
        image_frame = image_celestial_frame(imagename)
        if component_frame != image_frame:
                raise ValueError(
                    f'Final fitted direction frame {component_frame} disagrees with '
                    f'image frame {image_frame} for {imagename}'
                )
        ra_deg = normalize_right_ascension_degrees(ra_deg)
        if not 0.0 <= ra_deg < 360.0:
                raise ValueError(
                    f'Final fitted native right ascension normalization failed for '
                    f'{imagename}: {ra_deg!r}'
                )
        if not math.isfinite(dec_deg) or not -90.0 <= dec_deg <= 90.0:
                raise ValueError(
                    f'Final fitted native declination is invalid for {imagename}: '
                    f'{dec_deg!r}'
                )
        return component_frame, ra_deg, dec_deg


def measure_final_publication_noise(
    imagename,
    fitted_component,
    *,
    source,
    telescope,
    frequency_ghz,
    flux_mjy,
):
        """Measure final local statistics and adopt their scaled MAD.

        CASAdocs defines ``centerbox`` dimensions as full side widths but CRTF
        ellipse axes as semi-axes.  The final 3-sigma/FWHM fit square is only
        1.274 BMAJ wide and is therefore swallowed by its own 3-sigma ellipse.
        The noise aperture is instead a documented, source-centered square
        with full width 10 BMAJ, from which the fitted 3-sigma ellipse is
        subtracted.  There is no automatic aperture expansion or fallback.
        """

        header = imhead(imagename=imagename)
        beam_major, beam_minor, beam_pa = _restoring_beam_arcsec(header, imagename)
        pixel_width, pixel_height = _direction_pixel_widths_arcsec(header, imagename)
        image_width, image_height = _direction_image_shape_pixels(header, imagename)

        shape = fitted_component["shape"]
        fit_major = angular_value_to_arcsec(
            shape["majoraxis"]["value"], shape["majoraxis"]["unit"]
        )
        fit_minor = angular_value_to_arcsec(
            shape["minoraxis"]["value"], shape["minoraxis"]["unit"]
        )
        fit_pa_unit = str(shape["positionangle"]["unit"]).strip().casefold()
        if fit_pa_unit != "deg":
                raise ValueError(
                    f'Fitted position angle is not in degrees for {imagename}'
                )
        fit_pa = float(shape["positionangle"]["value"])
        fixed_shape = require_fixed_restoring_beam_shape(
            fit_major,
            fit_minor,
            fit_pa,
            beam_major,
            beam_minor,
            beam_pa,
        )

        center_x = float(fitted_component["pixelcoords"][0])
        center_y = float(fitted_component["pixelcoords"][1])
        if not math.isfinite(center_x) or not math.isfinite(center_y):
                raise ValueError(f'Final fitted center is not finite for {imagename}')
        fit_center_frame, fit_center_ra_deg, fit_center_dec_deg = (
            _native_fitted_direction_degrees(fitted_component, imagename)
        )
        aperture_width = square_bmaj_region_width_arcsec(
            beam_major,
            FINAL_NOISE_APERTURE_BMAJ_MULTIPLIER,
        )
        aperture_coverage = square_aperture_coverage(
            center_x,
            center_y,
            aperture_width,
            pixel_width,
            pixel_height,
            image_width,
            image_height,
        )
        mask_major = fit_major * THREE_SIGMA_FROM_FWHM
        mask_minor = fit_minor * THREE_SIGMA_FROM_FWHM
        aperture_region = (
            f'centerbox[[{center_x}pix,{center_y}pix],'
            f'[{aperture_width}arcsec,{aperture_width}arcsec]]'
        )
        source_mask = (
            f'ellipse[[{center_x}pix,{center_y}pix],'
            f'[{mask_major}arcsec,{mask_minor}arcsec],{fit_pa}deg]'
        )

        aperture_statistics = image_statistics_mjy_per_beam(
            imagename,
            aperture_region,
        )
        noise_region_path = None
        try:
                with tempfile.NamedTemporaryFile(
                    mode="w",
                    encoding="utf-8",
                    suffix=".crtf",
                    prefix="final-publication-noise-",
                    dir=".",
                    delete=False,
                ) as noise_region:
                        noise_region.write(
                            f'#CRTFv0\n{aperture_region}\n-{source_mask}\n'
                        )
                        noise_region_path = noise_region.name
                noise_statistics = image_statistics_mjy_per_beam(
                    imagename,
                    noise_region_path,
                )
        finally:
                if noise_region_path is not None and os.path.exists(noise_region_path):
                        os.remove(noise_region_path)

        pixel_area = pixel_width * pixel_height
        beam_area = gaussian_beam_area_arcsec2(beam_major, beam_minor)
        npts = noise_statistics["npts"]
        beam_areas = independent_beam_areas(
            npts,
            pixel_width,
            pixel_height,
            beam_major,
            beam_minor,
        )
        effective_samples = effective_independent_noise_samples(
            npts,
            pixel_width,
            pixel_height,
            beam_major,
            beam_minor,
        )
        review_flags = final_noise_review_flags(
            noise_statistics["mean_mjy_per_beam"],
            noise_statistics["sigma_mjy_per_beam"],
            noise_statistics["rms_mjy_per_beam"],
            noise_statistics["scaled_mad_mjy_per_beam"],
            aperture_minimum_edge_clearance_pixels=(
                aperture_coverage.minimum_edge_clearance_pixels
            ),
        )
        scaled_mad = noise_statistics["scaled_mad_mjy_per_beam"]
        record = _blank_final_noise_audit_record()
        record.update(
            {
                "schema_version": "1",
                "run_context": summary_prefix,
                "source": source,
                "telescope": telescope,
                "frequency_ghz": f'{float(frequency_ghz):.17g}',
                "product": portable_release_product(imagename),
                "forced_policy": "no",
                "measurement_scope": "publication_final_10bmaj_source_excluded",
                "fit_selection_policy": (
                    "single_initial_fixed_beam_fit_no_widened_retry"
                ),
                "aperture_center_policy": "fitted_component_pixelcoords",
                "aperture_center_x_pixel": repr(center_x),
                "aperture_center_y_pixel": repr(center_y),
                "fit_center_reference_frame": fit_center_frame,
                "fit_center_ra_deg": repr(fit_center_ra_deg),
                "fit_center_dec_deg": repr(fit_center_dec_deg),
                "aperture_width_arcsec": repr(aperture_width),
                "aperture_pre_exclusion_npts": str(
                    aperture_statistics["npts"]
                ),
                # CASA rasterizes the continuous square to selected pixels, so
                # this diagnostic area ratio can be slightly greater than one.
                "aperture_selected_pixel_area_over_requested_area": repr(
                    aperture_statistics["npts"]
                    * pixel_area
                    / (aperture_width * aperture_width)
                ),
                "mask_semantics": "crtf_ellipse_axes_are_semiaxes",
                "mask_major_semiaxis_arcsec": repr(mask_major),
                "mask_minor_semiaxis_arcsec": repr(mask_minor),
                "mask_pa_deg": repr(fit_pa),
                "beam_major_arcsec": repr(beam_major),
                "beam_minor_arcsec": repr(beam_minor),
                "beam_pa_deg": repr(beam_pa),
                "fit_beam_major_delta_arcsec": repr(
                    fixed_shape.major_delta_arcsec
                ),
                "fit_beam_minor_delta_arcsec": repr(
                    fixed_shape.minor_delta_arcsec
                ),
                "fit_beam_pa_delta_deg": repr(fixed_shape.pa_delta_deg),
                "fit_beam_pa_guard_applicable": (
                    "yes" if fixed_shape.pa_guard_applicable else "no"
                ),
                "image_width_pixels": str(image_width),
                "image_height_pixels": str(image_height),
                "aperture_half_width_x_pixels": repr(
                    aperture_coverage.half_width_x_pixels
                ),
                "aperture_half_width_y_pixels": repr(
                    aperture_coverage.half_width_y_pixels
                ),
                "aperture_left_clearance_pixels": repr(
                    aperture_coverage.left_clearance_pixels
                ),
                "aperture_right_clearance_pixels": repr(
                    aperture_coverage.right_clearance_pixels
                ),
                "aperture_bottom_clearance_pixels": repr(
                    aperture_coverage.bottom_clearance_pixels
                ),
                "aperture_top_clearance_pixels": repr(
                    aperture_coverage.top_clearance_pixels
                ),
                "aperture_minimum_edge_clearance_pixels": repr(
                    aperture_coverage.minimum_edge_clearance_pixels
                ),
                "aperture_edge_status": (
                    "truncated_available_valid_intersection"
                    if aperture_coverage.minimum_edge_clearance_pixels < -1e-9
                    else "fully_contained"
                ),
                "npts": str(npts),
                "pixel_width_arcsec": repr(pixel_width),
                "pixel_height_arcsec": repr(pixel_height),
                "pixel_area_arcsec2": repr(pixel_area),
                "valid_unmasked_area_arcsec2": repr(npts * pixel_area),
                "gaussian_beam_area_arcsec2": repr(beam_area),
                "independent_beam_areas": repr(beam_areas),
                "effective_independent_samples": repr(effective_samples),
                "minimum_independent_beams": repr(
                    FINAL_NOISE_MIN_INDEPENDENT_BEAMS
                ),
                "brightness_unit": noise_statistics["brightness_unit"],
                "rms_mjy_per_beam": repr(
                    noise_statistics["rms_mjy_per_beam"]
                ),
                "sigma_mjy_per_beam": repr(
                    noise_statistics["sigma_mjy_per_beam"]
                ),
                "mean_mjy_per_beam": repr(
                    noise_statistics["mean_mjy_per_beam"]
                ),
                "scaled_mad_mjy_per_beam": repr(scaled_mad),
                "adopted_noise_estimator": "scaled_mad",
                "adopted_noise_mjy_per_beam": repr(scaled_mad),
                "mean_over_sigma": repr(
                    abs(noise_statistics["mean_mjy_per_beam"])
                    / noise_statistics["sigma_mjy_per_beam"]
                ),
                "rms_over_scaled_mad": (
                    repr(noise_statistics["rms_mjy_per_beam"] / scaled_mad)
                    if scaled_mad > 0.0
                    else "inf"
                ),
                "review_flags": "+".join(review_flags) or "none",
                "review_disposition": (
                    "diagnostic_flags_recorded_not_classification_overrides"
                    if review_flags
                    else "no_diagnostic_flags"
                ),
                "flux_mjy": repr(float(flux_mjy)),
                "row_disposition": "candidate_pending",
                "guard_status": (
                    "passed"
                    if effective_samples >= FINAL_NOISE_MIN_INDEPENDENT_BEAMS
                    else "failed_insufficient_independent_samples"
                ),
            }
        )
        return scaled_mad, effective_samples, record


def record_forced_sma_noise_audit(
    imagename,
    statistics,
    *,
    source,
    telescope,
    frequency_ghz,
    flux_mjy,
    fit_converged,
    limit_mjy,
    reported_flux_field,
):
        """Record, without changing, the reviewed dirty-SMA one-beam policy."""

        # ``type`` is an inherited module-global plotting selector in this
        # script, so calling it here would shadow Python's builtin ``type``.
        # The upstream diagnostic deliberately normalizes CASA/NumPy boolean
        # scalars to one of the two builtin boolean singletons.
        if fit_converged is not True and fit_converged is not False:
            raise RuntimeError(
                "Forced dirty-SMA fit convergence flag must be boolean"
            )
        if fit_converged and flux_mjy is None:
            raise RuntimeError(
                "Converged forced dirty-SMA fit lacks its diagnostic flux"
            )
        if not fit_converged and flux_mjy is not None:
            raise RuntimeError(
                "Nonconverged forced dirty-SMA fit must not carry a diagnostic flux"
            )
        if fit_converged:
            try:
                fitted_flux_value = float(flux_mjy)
            except (TypeError, ValueError) as error:
                raise RuntimeError(
                    "Converged forced dirty-SMA diagnostic flux is invalid"
                ) from error
            if not math.isfinite(fitted_flux_value):
                raise RuntimeError(
                    "Converged forced dirty-SMA diagnostic flux is non-finite"
                )

        header = imhead(imagename=imagename)
        beam_major, beam_minor, beam_pa = _restoring_beam_arcsec(header, imagename)
        pixel_width, pixel_height = _direction_pixel_widths_arcsec(header, imagename)
        pixel_area = pixel_width * pixel_height
        beam_area = gaussian_beam_area_arcsec2(beam_major, beam_minor)
        npts = statistics["npts"]
        beam_areas = independent_beam_areas(
            npts,
            pixel_width,
            pixel_height,
            beam_major,
            beam_minor,
        )
        effective_samples = effective_independent_noise_samples(
            npts,
            pixel_width,
            pixel_height,
            beam_major,
            beam_minor,
        )
        record = _blank_final_noise_audit_record()
        record.update(
            {
                "schema_version": "1",
                "run_context": summary_prefix,
                "source": source,
                "telescope": telescope,
                "frequency_ghz": f'{float(frequency_ghz):.17g}',
                "product": portable_release_product(imagename),
                "forced_policy": "yes_dirty_sma_nitercal_0",
                "measurement_scope": "reviewed_forced_dirty_sma_one_beam_region",
                "fit_selection_policy": (
                    FORCED_FIT_POLICY_CONVERGED
                    if fit_converged
                    else FORCED_FIT_POLICY_NONCONVERGED
                ),
                "aperture_center_policy": "catalog_or_anchor_coordinate",
                "aperture_width_arcsec": "one_restoring_beam_rotbox",
                "mask_semantics": "none_source_not_excluded_reviewed_forced_policy",
                "beam_major_arcsec": repr(beam_major),
                "beam_minor_arcsec": repr(beam_minor),
                "beam_pa_deg": repr(beam_pa),
                "npts": str(npts),
                "pixel_width_arcsec": repr(pixel_width),
                "pixel_height_arcsec": repr(pixel_height),
                "pixel_area_arcsec2": repr(pixel_area),
                "valid_unmasked_area_arcsec2": repr(npts * pixel_area),
                "gaussian_beam_area_arcsec2": repr(beam_area),
                "independent_beam_areas": repr(beam_areas),
                "effective_independent_samples": repr(effective_samples),
                "minimum_independent_beams": "not_applicable",
                "brightness_unit": statistics["brightness_unit"],
                "rms_mjy_per_beam": repr(statistics["rms_mjy_per_beam"]),
                "sigma_mjy_per_beam": repr(statistics["sigma_mjy_per_beam"]),
                "mean_mjy_per_beam": repr(statistics["mean_mjy_per_beam"]),
                "scaled_mad_mjy_per_beam": repr(
                    statistics["scaled_mad_mjy_per_beam"]
                ),
                "adopted_noise_estimator": "classic_rms_reviewed_forced_policy",
                "adopted_noise_mjy_per_beam": repr(
                    statistics["rms_mjy_per_beam"]
                ),
                "flux_mjy": (
                    repr(fitted_flux_value) if fit_converged else ""
                ),
                "classification": "forced_upper_limit",
                "reported_limit_mjy": repr(float(limit_mjy)),
                "reported_flux_field": reported_flux_field,
                "review_flags": "reviewed_forced_policy",
                "review_disposition": "preapproved_forced_policy",
                "row_disposition": "candidate_pending",
                "guard_status": "not_applicable_reviewed_forced_policy",
            }
        )
        record_final_noise_audit(record)


#print statements and make picturesmakedotcontour
def printrueconv(imfitter,date,make,color,imfitloc,newfilename,printimages,beamsize,freq,oicoords,frequnit,telescopes, tele, file, source, getchoosercoords, fixtochosenfreq, imagefound, prelimchecker, fit_sizenoloc, doublegaussfit,writenoise, *, force_upper_limit=False):
    ordinary_noise_scope = (
        'publication_final_10bmaj_source_excluded'
        if fixtochosenfreq
        else 'intermediate_whole_image'
    )
    if doublegaussfit:
        low=[]
        detectregions=[]
    else:
        low=False
    ordinary_upper_limit_mjy=None
    nelem=imfitter['deconvolved']['nelements']
    counter=-1
    for inelem in range(nelem):
        if doublegaussfit:
            ilow=False
        counter=counter+1
        component=f'component{inelem}'
        fluxun=imfitter['deconvolved'][component]['flux']['unit']
        fluxerr=flux_density_to_mjy(
                imfitter['deconvolved'][component]['flux']['error'][0], fluxun
        )
        fluxval=flux_density_to_mjy(
                imfitter['deconvolved'][component]['flux']['value'][0], fluxun
        )

        #FWHM to sigma from wikipedia
        factorto3sig= 3 / (2 * np.sqrt(2 * np.log(2)))

        fitmajax=imfitter['results'][component]['shape']['majoraxis']['value']
        fitmajax=fitmajax*factorto3sig
        fitmajunit=imfitter['results'][component]['shape']['majoraxis']['unit']
        if fitmajunit=='deg':
            fitmajax=fitmajax*60*60
        elif fitmajunit!='arcsec':
            raise ValueError(f'Unsupported fitted major-axis unit {fitmajunit!r} for {file}')


        fitminjax=imfitter['results'][component]['shape']['minoraxis']['value']
        fitminjax=fitminjax*factorto3sig
        fitminjunit=imfitter['results'][component]['shape']['minoraxis']['unit']
        if fitminjunit=='deg':
            fitminjax=fitminjax*60*60
        elif fitminjunit!='arcsec':
            raise ValueError(f'Unsupported fitted minor-axis unit {fitminjunit!r} for {file}')

        fitpos=imfitter['results'][component]['shape']['positionangle']['value']
        fitposun=imfitter['results'][component]['shape']['positionangle']['unit']
        if fitposun!='deg':
            raise ValueError(f'Unsupported fitted position-angle unit {fitposun!r} for {file}')

        gausscen1=imfitter['results'][component]['pixelcoords'][0]
        gausscen2=imfitter['results'][component]['pixelcoords'][1]

        if not doublegaussfit:
            if force_upper_limit:
                # This reviewed dirty-SMA branch deliberately precedes the
                # ordinary classifier and its final-aperture guard.  Its limit
                # remains three times the historical one-beam-region RMS.
                low=True
                noisesnr=math.nan
            else:
                try:
                    detectregion=f'ellipse[[{gausscen1}pix, {gausscen2}pix], [{fitmajax}arcsec, {fitminjax}arcsec], {fitpos}deg]'
                    if ordinary_noise_scope == 'publication_final_10bmaj_source_excluded':
                        noiseval,effective_samples,audit_record=(
                            measure_final_publication_noise(
                                file,
                                imfitter['results'][component],
                                source=source,
                                telescope=tele,
                                frequency_ghz=freq,
                                flux_mjy=fluxval,
                            )
                        )
                        # Persist the full inventory before enforcing the hard
                        # area gate so a stopped run identifies the exact row.
                        record_final_noise_audit(audit_record)
                        require_minimum_independent_beams(effective_samples)
                    else:
                        # Preserve the reviewed discovery/anchor classifier.
                        # Extended preliminary components can swallow their
                        # fit squares, so these non-publication passes retain
                        # the whole-image-minus-source region exactly.
                        totstats = imstat(imagename=file)
                        blc = totstats['blc']
                        trc = totstats['trc']
                        noise_base_region = (
                            f'box[[{blc[0]}pix, {blc[1]}pix], '
                            f'[{trc[0]}pix, {trc[1]}pix]]'
                        )
                        noiseregfile='noiseregfile.crtf'
                        if os.path.exists(noiseregfile):
                            os.remove(noiseregfile)
                        with open(noiseregfile, 'w') as noise_region_handle:
                            noise_region_handle.write(
                                f'#CRTFv0\n{noise_base_region}\n-{detectregion}\n'
                            )
                        noiseval=image_rms_mjy_per_beam(file, 'noiseregfile.crtf')

                    classification=classify_ordinary_fit(fluxval, noiseval)
                    noisesnr=classification.snr
                    ordinary_upper_limit_mjy=(
                        classification.upper_limit_mjy
                        if ordinary_noise_scope == 'publication_final_10bmaj_source_excluded'
                        else None
                    )
                    low=not classification.is_detection
                    if ordinary_noise_scope == 'publication_final_10bmaj_source_excluded':
                        update_final_noise_audit(
                            file,
                            snr=repr(classification.snr),
                            classification=(
                                "detection"
                                if classification.is_detection
                                else "ordinary_upper_limit"
                            ),
                            ordinary_upper_limit_mjy=repr(
                                classification.upper_limit_mjy
                            ),
                            reported_limit_mjy=(
                                ""
                                if classification.is_detection
                                else repr(classification.upper_limit_mjy)
                            ),
                            reported_flux_field=(
                                format_number(fluxval)
                                if classification.is_detection
                                else f'{format_number(classification.upper_limit_mjy)}*'
                            ),
                        )

                except Exception as e:
                    raise RuntimeError(
                        f'Cannot measure the source-excluded local noise for ordinary fit '
                        f'{file} in scope {ordinary_noise_scope}'
                    ) from e

            fluxval=format_number(fluxval)
            sfluxval=fluxval
            sbeamsize=format_number(beamsize)

            if low!=True:
                if counter==0:
                    if make==True:
                        pixlen=100
                        date=makedotcontour(file,imfitloc,sfluxval,sbeamsize, newfilename, printimages,color,pixlen,oicoords,freq,frequnit,imfitter,telescopes,tele, getchoosercoords, fixtochosenfreq,imagefound, prelimchecker, fit_sizenoloc, doublegaussfit)
                        if date==None:
                            if fixtochosenfreq:
                                update_final_noise_audit(
                                    file, row_disposition="rejected_missing_plot_date"
                                )
                            return(date,date,ordinary_upper_limit_mjy)
                        if date=='ONBOUNDRY':
                            if fixtochosenfreq:
                                record_reviewed_boundary_rejection(file)
                            return(date,date,ordinary_upper_limit_mjy)
                        stotarea=f'{beamsize:.5f}'
                        sfreq=f'{freq:.2f}'
                        ssnr=str(f'{noisesnr:.3f}')
                        summary_product=portable_release_product(file)
                        print(f"{source.ljust(15)}{tele.ljust(15)}{sfreq.ljust(19)}{sfluxval.ljust(16)}{date.ljust(15)}{stotarea.ljust(19)}{ssnr.ljust(15)}{summary_product}")
                        write_fluxerr_row(source, tele, freq, date, fluxerr, summary_product)
                        if fixtochosenfreq:
                            update_final_noise_audit(
                                file,
                                row_disposition=(
                                    "rejected_source17_requires_source18_join"
                                    if source == "NGC4945"
                                    else "generic_table_retained"
                                ),
                            )
                        if make!=True:
                            date='NA'
                        return(low,date,ordinary_upper_limit_mjy)
                    else:
                        input('multcoords function being used but shouldnt be')
            else:
                return(low,date,ordinary_upper_limit_mjy)
        else:
            idetectregion='NA'
            try:
                idetectregion=f'ellipse[[{gausscen1}pix, {gausscen2}pix], [{fitmajax}arcsec, {fitminjax}arcsec], {fitpos}deg]'
            except Exception as e:
                print(f'cannot make noise region for {file}:')
                print(e)
            detectregions.append(idetectregion)





    if doublegaussfit:
        #the first gaussian is fixed to the chosen frequency
        #the second gaussian is free to vary, and centered on the location
        #chosen from the 'noloc' pipeline

        #the second gaussian SHOULD always converge (needs confirmation), and the first one
        #may or may not. if not, it should be considered a nondetection
        #the upperlimit should still be built into the pipeline

        #totdetectregion=f'{detectregions[0]}+{detectregions[1]}'
        #noiseregion=f'{totregion}-{totdetectregion}'

        totstats = imstat(imagename=file)
        blc = totstats['blc']  # bottom-left corner: [x, y]
        trc = totstats['trc']  # top-right corner: [x, y]
        totregion = f'box[[{blc[0]}pix, {blc[1]}pix], [{trc[0]}pix, {trc[1]}pix]]'
        noiseregfile='noiseregfile.crtf'
        if os.path.exists(noiseregfile):
            os.remove(noiseregfile)
        with open(noiseregfile, 'w') as noiseregfile:
            if writenoise==True:
                original_stdoutnoisereg = sys.stdout
                sys.stdout = noiseregfile
            print('#CRTFv0')
            print(totregion)
            print(f'-{detectregions[0]}')
            print(f'-{detectregions[1]}')
            if writenoise==True:
                sys.stdout = original_stdoutnoisereg

        noiseval=image_rms_mjy_per_beam(file, 'noiseregfile.crtf')

        colors=[]
        counter=-1
        sfluxarray=[]
        for inelem in range(nelem):
            counter=counter+1
            component=f'component{inelem}'
            fluxun=imfitter['deconvolved'][component]['flux']['unit']
            fluxerr=flux_density_to_mjy(
                    imfitter['deconvolved'][component]['flux']['error'][0], fluxun
            )
            fluxval=flux_density_to_mjy(
                    imfitter['deconvolved'][component]['flux']['value'][0], fluxun
            )
            noisesnr=background_rms_snr(fluxval, noiseval)
            if noisesnr<3:
                low=True
                print(f'bad double fit for {file}')
                return(low,date,None)
            if make==True:
                imajval = imfitter['results'][component]['shape']['majoraxis']['value']
                iminval = imfitter['results'][component]['shape']['minoraxis']['value']
                imajunit = imfitter['results'][component]['shape']['majoraxis']['unit']
                iminunit = imfitter['results'][component]['shape']['minoraxis']['unit']

                if imajunit!='arcsec':
                    if imajunit=='deg':
                        imajval=imajval*60*60
                        imajunit='arcsec'
                    else:
                        'check units'
                if iminunit!='arcsec':
                    if iminunit=='deg':
                        iminval=iminval*60*60
                        iminunit='arcsec'
                    else:
                        'check units'

                ibeamsize=np.sqrt(imajval*iminval)
                if fluxval<noiseval*3:
                    fluxval=noiseval*3
                    color='r'
                    colors.append('r')
                else:
                    color='b'
                    colors.append('b')
                sfluxval=format_number(fluxval)
                if color=='r':
                    sfluxval=sfluxval+'*'
                sfluxarray.append(sfluxval)
                sbeamsize=format_number(ibeamsize)

                stotarea=f'{ibeamsize:.5f}'

                sfreq=f'{freq:.2f}'

                ssnr=str(f'{noisesnr:.3f}')
                scounter=str(counter)
                summary_product=portable_release_product(file)
                print(f"{source.ljust(15)}{scounter.ljust(5)}{tele.ljust(15)}{sfreq.ljust(19)}{sfluxval.ljust(16)}{date.ljust(15)}{stotarea.ljust(19)}{ssnr.ljust(15)}{summary_product}")
                write_fluxerr_row(source, tele, freq, date, fluxerr, summary_product, component=counter)

        if make==True:
            pixlen=100
            date=makedotcontour(file,imfitloc,sfluxarray,sbeamsize, newfilename, printimages,colors,pixlen,oicoords,freq,frequnit,imfitter,telescopes,tele, getchoosercoords, fixtochosenfreq,imagefound, prelimchecker, fit_sizenoloc, doublegaussfit)
            if date==None:
                return(date,date,None)
            if date=='ONBOUNDRY':
                return(date,date,None)
        if make!=True:
            date='NA'


        return(low,date,None)


# In[ ]:


#function: check case of last letter
def ensure_last_char_uppercase(s):
    # Check if the string is not empty
    if s:
        # Get the last character
        last_char = s[-1]
        # Check if it is an alphabetical letter
        if last_char.isalpha():
            # Ensure it is uppercase
            s = s[:-1] + last_char.upper()
    return s


# In[ ]:


#get FOV and beam
def compute_fov_and_beam(imagename):
      # Get header
      header = imhead(imagename)

      # Pixel dimensions
      nx = header['shape'][0]  # X (RA)
      ny = header['shape'][1]  # Y (Dec)

      # Pixel increment and units
      incr_x = abs(header['incr'][0])
      incr_y = abs(header['incr'][1])
      unit_x = header['axisunits'][0]
      unit_y = header['axisunits'][1]

      # Convert to arcsec
      def to_arcsec(value, unit):
            if unit == 'rad':
                  return value * 206265
            elif unit == 'deg':
                  return value * 3600
            elif unit == 'arcsec':
                  return value
            else:
                  raise ValueError(f"Unrecognized angular unit: '{unit}'")

      dx = to_arcsec(incr_x, unit_x)
      dy = to_arcsec(incr_y, unit_y)

      # Field of view
      fov_x = nx * dx
      fov_y = ny * dy
      fov_diag = (fov_x**2 + fov_y**2)**0.5

      # Beam parameters with unit check
      beam = header['restoringbeam']
      bmaj_val = beam['major']['value']
      bmin_val = beam['minor']['value']
      bmaj_unit = beam['major']['unit']
      bmin_unit = beam['minor']['unit']
      bpa = beam['positionangle']['value']  # degrees

      bmaj = to_arcsec(bmaj_val, bmaj_unit)
      bmin = to_arcsec(bmin_val, bmin_unit)
      bdiag = (bmaj**2 + bmin**2)**0.5  # beam diagonal in arcsec


      scalesize=fov_diag/bdiag

      return scalesize,fov_diag,bdiag



# In[ ]:


#MAIN ENGINE

#now we look for the best image sources. small beam detections
choosercoords_by_source={}
#nofixmultval must be adjust for certain sources
nofixmultval_exc=[['NGC3393', 1,]]
#nofixmultval_exc=[['NGC3393', 1,],['NGC2273',2.5]]
#nofixmultval_exc=[['x',0],['y',0]]

#this fits location to the chosen freq coordinates with a fixed beam size fixed location
#this method sets the fitting region to the size of the beam... fixmultval is the scale of the entire fitting region
#This value goes from FWHM to the 3 sigma of the beam
fixmultval=(3)/(2*np.sqrt(2*np.log(2)))

#too complicated for this data set
#set doublegaussfit to create a 2 componet gaussian fit.  The first component is
#from fixtochosenfreq, which will be fixed and the value of interest. This will appear
#as the first value (the first of 2 gaussians), the second gaussian is allowed to vary
#although the intial guess is that from the fit_sizenoloc fit file
###
#in the case that fixtochosenfreq is a nondetection and fit_sizenoloc is a detection, we will run the double gaussian fit
#there should be no case where fixtochosenfreq is a detection and fit_sizenoloc is a detection
#if both are non detections... we resort to the upper limit in fixtochosenfreq
doublegaussfit=0
###
###
###
###
debugger=0
#debugger=1

yestext=bool(args.enable_notifications)

filechooser=False
#filechooser=True
filechoice='VLA/NGC3079/43.3I0.04_AG0527_1998FEB25_1_346.U26.6S.imfits'

biggerregion='VLA/NGC3079/43.3I0.04_AG0527_1998FEB25_1_346.U26.6S.imfits'

# Reviewed archive exclusions and their immutable FITS hashes are controlled
# solely by this auditable policy; there is no second hard-coded bad-file list.
reviewed_archive_photometry_exclusions = load_reviewed_exclusion_policy(
    Path(__file__).with_name("continuum_product_exclusion_policy.csv")
)
reviewed_boundary_rejection_rows = load_boundary_rejection_policy(
    Path(__file__).with_name("continuum_boundary_rejection_policy.csv"),
    analysis_root=Path("."),
)
reviewed_boundary_rejections_by_product = {
    row["product"]: row for row in reviewed_boundary_rejection_rows
}

#this naming will be from the allsources array format
sourcechooser=False
#sourcechooser=True
#sourcechoice='NGC4258'

printfile=False
#printfile=True

startwith=False
#startwith=True
#sourcestart='NGC4258'
#startval=0

###
###
###
###

#text files are used for creating regions
#fit estimates
writerreg=1
#writerreg=0

#noise region for noise region 
writenoise=1
#writenoise=0

###
###
###
###
#fix the location of the gaussian. 0 for not fix 1 for fix.
###
###
###
###

if prelimchecker+fit_sizenoloc+getchoosercoords+fixtochosenfreq+doublegaussfit!=1:
    input('error in setup')

pipetype='NA'
#freely chosen fit
if prelimchecker==1:
    masterfixer='nofix'
    pipetype='prelim'

#fix the beam size (for all sources)
if fit_sizenoloc==1:
    masterfixer=0
    pipetype='noloc'

#fix the beam size (chosen freq)
if getchoosercoords==1:
    masterfixer=0
    pipetype='getchooser'

#fix the location and beam size
if fixtochosenfreq==1:
    masterfixer=0
    pipetype='fixedtochosen'

#fix the location and beam size
if doublegaussfit==1:
    masterfixer=1
    pipetype='doublegaussfit'
#the first gaussian is fixed to the chosen frequency
#the second gaussian is free to vary, and centered on the location
#chosen from the 'noloc' pipeline

#the second gaussian SHOULD always converge (needs confirmation), and the first one
#may or may not. if not, it should be considered a nondetection
#the upperlimit should still be built into the pipeline

###
###
###
###

from casatools import image
from casatasks import imfit
from casatasks import rmtables

import warnings
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.wcs import FITSFixedWarning
warnings.filterwarnings('ignore', category=FITSFixedWarning)

import os
from contextlib import contextmanager
import warnings
from astropy.wcs.wcs import FITSFixedWarning
from casatasks import casalog

warnings.filterwarnings('ignore', category=FITSFixedWarning)

@contextmanager
def silence_stderr_fd():
        saved_fd = os.dup(2)
        devnull = os.open(os.devnull, os.O_WRONLY)
        try:
                os.dup2(devnull, 2)
                yield
        finally:
                os.dup2(saved_fd, 2)
                os.close(saved_fd)
                os.close(devnull)


if prelimchecker:
    output_plot_filename=f'images/{flagnum}/prelim'
elif fit_sizenoloc:
    output_plot_filename=f'images/{flagnum}/noloc'
elif getchoosercoords:
    output_plot_filename=f'images/{flagnum}/chosenfreq'
elif fixtochosenfreq==1:
    output_plot_filename=f'images/{flagnum}/fixedtochosen'
elif doublegaussfit==1:
    output_plot_filename=f'images/{flagnum}/doublegauss'

if os.path.exists(output_plot_filename):
    shutil.rmtree(output_plot_filename)

txtnamesandcoords='namesandcoords.txt'

ia = image()

source_text='namesandcoords.txt'

#[[sourcename,[[arcsec of tot fit region,ra,dec]]]]
#multgauss=[[[[]]]]
#multgauss=[['NGC3393',[[3.0,'10:48:23.409','-25:09:44.05']]]]
multgauss=[[[[]]]]

#The code will keep only names that come before a dot/decimal: "."
#specify strings to kill
listofkills=['Galaxy']

telescopes=['VLA','ALMA','SMA']
alltelescopes=telescopes.copy()
if telechooser==1:
        telescopes=[telechoice]
indexfiles=[]

if os.path.exists('regrids'):
    shutil.rmtree('regrids')
os.makedirs('regrids')

if os.path.exists('guessregions'):
    shutil.rmtree('guessregions')
os.makedirs('guessregions')

#track file number
index=0
#holds the sources best resolution image
bestcoords=[]
badfiles=[]

os.makedirs('separations', exist_ok=True)
separation_outfile = f'separations/separations_{flagnum}_{pipetype}.txt'
if os.path.exists(separation_outfile):
    os.remove(separation_outfile)



#iterate through the 3 telescopes to get all sources
allsources=[]
oallsources=[]
for  tele in telescopes:
    if tele=='SMA':
        sources, original_sources = get_sma_source_names()
        for original_source, source_name in zip(original_sources, sources):
            oallsources.append(original_source)
            allsources.append(source_name)

    if tele=='VLA':
        sources=os.listdir(tele)
        for i in sources:
            oallsources.append(i)
            if len(i.split('.'))>0:
                i=i.split('.')[0]
            for killstring in listofkills:
                if killstring in i:
                    i=i.replace(killstring,'')
            allsources.append(i)

    if tele=='ALMA':
        sources=os.listdir(f'{tele}')
        for i in sources:
            oallsources.append(i)
            if len(i.split('.'))>0:
                i=i.split('.')[0]
            for killstring in listofkills:
                if killstring in i:
                    i=i.replace(killstring,'')
            allsources.append(i)

allsources=np.unique(allsources)
oallsources=np.unique(oallsources)

if not os.path.exists('fitsumfiles/witherrors'):
    os.makedirs('fitsumfiles/witherrors')
if not os.path.exists('fitsumfiles/fluxerrors'):
    os.makedirs('fitsumfiles/fluxerrors')
summaryfile=f'fitsumfiles/witherrors/{summary_prefix}_fitsummary_{pipetype}_witherrors.txt'
fluxerrfile=f'fitsumfiles/fluxerrors/{summary_prefix}_fluxerr_{pipetype}.txt'
if fixtochosenfreq:
    final_noise_audit_path = Path(
        f'fitsumfiles/noise_audits/'
        f'{summary_prefix}_final_publication_noise_audit.csv'
    )
    final_noise_audit_records.clear()
    if final_noise_audit_path.exists():
        final_noise_audit_path.unlink()

with open(summaryfile, 'w') as sumfile, open(fluxerrfile, 'w') as fluxerr_out:

    fluxerrfile_handle = fluxerr_out
    fluxerr_firstline(doublegaussfit)

    if writesummary==True:
        original_stdoutsum = sys.stdout
        sys.stdout = sumfile

    #print the first line
    firstline()
    newsource=-1
    for source in allsources:
        breaker=0
        if sourcechooser==True:
            if source!=sourcechoice:
                breaker=1

        if startwith==1:
            if source==sourcestart:
                startval=1
            if startval==0:
                breaker=1

        if breaker==1:
            continue
        newsource=newsource+1

        if not os.path.exists(source_text):
            input(f'{source_text} does not exist.  Go to the beginning cells and make it.')

        coordinate_record_count=0
        with open(source_text, 'r') as infile:
            for params in infile:
                params=eval(params)
                for inamecoords in params:
                    oname=inamecoords[0][0]
                    # The alternate-name list may legitimately be empty (for
                    # example, CircinusGalaxy).  Keep the complete list here;
                    # the coordinate lookup below is keyed by ``oname`` and
                    # does not require a first alternate name.
                    nname=inamecoords[1][0]
                    ra=inamecoords[1][1]
                    dec=inamecoords[1][2]
                    if not source in oname:
                        continue
                    icoords=[]
                    icoords.append(ra)
                    icoords.append(dec)
                    oicoords=icoords.copy()
                    coordinate_record_count=coordinate_record_count+1
        if coordinate_record_count!=1:
            raise RuntimeError(
                f'Expected exactly one coordinate record for {source}; '
                f'found {coordinate_record_count}'
            )
        allfile=[source,allfiles(source)]
        selectcoords=[]
        aselectcoords=[]


        #checking to see if bestcoords have already been obtained
        i=0
        for ibest in bestcoords:
            if ibest[0]==source:
                iibest=ibest
                i=i+1
        if i==0:
            namecoordtxt='namesandcoords.txt'
            if not os.path.exists(namecoordtxt):
                input(f'{namecoordtxt} Does Not Exist.  ')
            ibestcoords=[]
            if os.path.exists(namecoordtxt):
                if 0==0:
                    with open(txtnamesandcoords, 'r') as infile:
                        for params in infile:
                            params=eval(params)
                            for inamecoords in params:
                                oname=inamecoords[0]
                                nname=inamecoords[1][0]
                                ra=inamecoords[1][1]
                                dec=inamecoords[1][2]

                                icoord1=ra.split(':')
                                selectcoord1=f'{icoord1[0]}h{icoord1[1]}m{icoord1[2]}s'   
                                icoord2=dec.split(':')
                                selectcoord2=f'{icoord2[0]}d{icoord2[1]}m{icoord2[2]}s'

                                aicoords=au.J2000ToICRS(f'{ra} {dec}',verbose=False).split(', ')
                                aicoord1=aicoords[0].split(':')
                                aselectcoord1=f'{aicoord1[0]}h{aicoord1[1]}m{aicoord1[2]}s'
                                aicoord2=aicoords[1].split(':')
                                aselectcoord2=f'{aicoord2[0]}d{aicoord2[1]}m{aicoord2[2]}s'
                                sendback=[aselectcoord1,aselectcoord2],[selectcoord1,selectcoord2]
                                sendit=0
                                if source in oname:
                                    sendit=1
                                for checksource in oname:
                                    if source in checksource:
                                        sendit=1
                                if sendit==1:
                                    iibest=[]
                                    iibest.append(source)
                                    iibest.append(sendback)
                                    ibestcoords.append(iibest)

            if fixtochosenfreq==0:
                choosefile=0
                choicefile='NA'
                imagefound=0
            if fixtochosenfreq==1:
                choosefile=1
                folder = f'images/chosenfreq'
                imagefound=0
                telefound=0
                for item in os.listdir(folder):
                    if 'pdf' in item:
                        continue
                    if source==item:
                        full_path = os.path.join(folder, item)
                        print(folder)
                        print(item)
                        for item2 in os.listdir(full_path):
                            choicefile=os.path.join(full_path, item2)
                            imagefound=1
                            for findtele in alltelescopes:
                                if findtele.lower() in item2.lower():
                                    chosentele=findtele
                                    telefound=1
            if imagefound==1:
                if telefound!=1:
                    print(telescopes)
                    print(chosentele)
                    input('error in tele')
                # getbestrez returns coherent ICRS and J2000 coordinate pairs.
                #return(sendback)
                newval=[source,getbestrez(allfile,icoords,chosentele,telescopes,choosefile,choicefile)]
                if newval[1]!=None:
                    ibestcoords.append(newval)
                if newval[1]==None:
                    ibestcoords.append(f'{source} has no proper images')
            if imagefound==0:
                ibestcoords.append('no chosen frequency')
        bestcoords.append(ibestcoords)
        #verifying that bestcoords have been obtained
        i=0
        iibest=None
        for ibest in bestcoords:
            if ibest[0][0]==source:
            # if the first entry is the source (the NED query location)
                iibest=ibest
                i=i+1
        if i!=1 or iibest is None:
            raise RuntimeError(
                f'Expected exactly one assembled coordinate record for {source}; '
                f'found {i}'
            )

        #set icoords (the location of the Gaussian Fit) 
        # to the NED query if fixtochosenfreq==0 otherwise fix it to the location of the chosenfreq
        # if fixtochosenfreq==1 but there is no chosenfreq, the location from NED is chosen 
        if fixtochosenfreq==0:
            icoords=iibest[0][1]  
            acoords=icoords[0]
            coords=icoords[1]
        elif fixtochosenfreq==1 and imagefound==1:
            icoords=iibest[1]
            acoords=icoords[1][0]
            coords=icoords[1][1]
        elif fixtochosenfreq==1 and imagefound==0:
            icoords=iibest[0][1]
            acoords=icoords[0]
            coords=icoords[1]


        #fill the files for the proper filename and telescope    
        if newsource!=0:
            if findbest==0:
                #input(f'cannotfind best file for {allsources[newsource-1]}')
                pass
            else:
                #input(f'best file for {allsources[newsource-1]} is {bestfile}')
                pass
        bestfile='NA'
        findbest=0
        for tele in telescopes:
            if telechooser==True:
                if tele!=telechoice:
                    continue

            files=[]
            if yestext==True:
                textme(f'{source}:{tele}','go')

            if tele=='SMA':
                files = get_sma_image_files_for_source(source)

            if tele=='ALMA':
                files = get_archive_image_files('ALMA', source)

            if tele=='VLA':
                files = get_archive_image_files('VLA', source)

            #iterating through all of the files
            for file in files:
                portable_file = portable_release_product(file)
                expected_exclusion_hash = reviewed_archive_photometry_exclusions.get(
                    portable_file
                )
                if expected_exclusion_hash is not None:
                    verify_reviewed_exclusion_hash(file, expected_exclusion_hash)
                    badfiles.append([8, file])
                    continue
                if debugger:
                        print('1')
                #skip bad files flagged from prelimchecker
                breaker=0
                if prelimchecker==1:
                    scalesize=compute_fov_and_beam(file)
                    scalesizeval=scalesize[0]
                    if scalesizeval>=60:
                        badfiles.append(['bigsize',file])
                        plt.close()
                        continue
                else:
                    if fit_sizenoloc==1:
                        badfileversions=['prelim']
                    elif getchoosercoords==1:
                        badfileversions=['prelim','noloc']
                    elif fixtochosenfreq==1:
                        badfileversions=['prelim','noloc']
                    elif doublegaussfit==1:
                        badfileversions=['prelim','noloc']
                    for ibad in badfileversions:
                        badfiletxt=f'bad/masterbadfiles/masterbadfiles_{ibad}.txt'
                        if os.path.exists(badfiletxt):
                            with open(badfiletxt, 'r') as lines:
                                for line in lines:
                                    line=line.replace('\n','')
                                    if line==file:
                                        breaker=1
                if breaker==1:
                    continue

                if printfile==True:
                    print(file)
                breaker=0
                if filechooser==True:
                    if file!=filechoice:
                        breaker=1
                if breaker==1:
                    continue
                if not os.path.exists(file):
                    input(f'{file} doesnt exist')
                #get frequency, date
                #print(imhead(imagename=file)['axisnames'])
                axisnames=imhead(imagename=file)['axisnames']
                i=-1
                for iaxisname in axisnames:
                    i=i+1
                    if iaxisname=='Frequency':
                        findfreq=i
                if i==-1:
                    input('error1')
                freq=imhead(imagename=file)['refval'][findfreq]/10**9
                if getchoosercoords:
                    if not 200<=freq<=400:
                        continue
                if imhead(imagename=file)['axisunits'][findfreq]!='Hz':
                    input('error2')         



                #obtain beam properties for this image
                imheader=imhead(imagename=file)
                celestial_frame = image_celestial_frame(file)
                active_coords = coordinates_for_image(file, acoords, coords)
                if debugger:
                    print('2')
                majunit=imheader['restoringbeam']['major']['unit']
                majval=imheader['restoringbeam']['major']['value']
                minunit=imheader['restoringbeam']['minor']['unit']
                minval=imheader['restoringbeam']['minor']['value']
                rotunit=imheader['restoringbeam']['positionangle']['unit']
                if rotunit!='deg':
                    print(f'Beam Rotation does not use deg for units.  check imhead for: {file}')
                rotval=imheader['restoringbeam']['positionangle']['value']
                if majunit!='arcsec':
                    if majunit=='deg':
                        majval=majval*60*60
                        majunit='arcsec'
                    else:
                        'check units'
                if minunit!='arcsec':
                    if minunit=='deg':
                        minval=minval*60*60
                        minunit='arcsec'
                    else:
                        'check units'
                #imhead appear to give diameter values (FWHM)
                beamsize=np.sqrt(majval*minval)
                if getchoosercoords:
                    if beamsize>=1:
                        continue

                # The position search intentionally uses an isotropic square.
                # Its full side width is derived conservatively from BMAJ; it
                # is not an ellipse and BMIN is not a region-width parameter.
                if fixtochosenfreq==1 and imagefound==1:
                    multipliertot=fixmultval
                elif fixtochosenfreq==1 and imagefound==0:
                    multipliertot=nofixmultval*2
                else:
                    multipliertot=nofixmultval

                minimum_width = 1 if getchoosercoords==1 else None
                maximum_width = (
                    3 if fixtochosenfreq==1 and imagefound==0 else None
                )
                override_width = (
                    3 if fixtochosenfreq==1 and imagefound==0 and file==biggerregion
                    else None
                )
                fit_region_width_arcsec = square_bmaj_region_width_arcsec(
                    majval,
                    multipliertot,
                    minimum_arcsec=minimum_width,
                    maximum_arcsec=maximum_width,
                    override_arcsec=override_width,
                )

                if getchoosercoords==1:
                    #nofixmultval must be adjust for certain sources
                    for inofix in nofixmultval_exc:
                        if inofix[0].lower() in source.lower():
                            fit_region_width_arcsec=square_bmaj_region_width_arcsec(
                                majval, 1, override_arcsec=inofix[1]
                            )


                #if the image is a very high resolution than the coordinates from
                #NED may not be a good guess... in which case we will use the image from the highest obs
                #if getchoosercoords==1:
                    pass
                    #choosefile=1
                    #choicefile=file
                    #xacoords,xcoords,checkthisfile=getbestrez('NA',coords,tele,'NA',choosefile,choicefile)                
                #elif getchoosercoords==1:
                    #choosefile=1
                    #choicefile=file
                    #xacoords,xcoords,checkthisfile=getbestrez('NA',coords,tele,'NA',choosefile,choicefile)                
                #else:
                    #input(f'error for {source} in getchoosercoords')


                #this is the region of the beam                 
                regionsmall=f'rotbox[{coords},[{majval}arcsec,{minval}arcsec],{rotval}{rotunit}]'
                aregionsmall=f'rotbox[{acoords},[{majval}arcsec,{minval}arcsec],{rotval}{rotunit}]'

                namefreq=int(freq)
                frequnit='GHz'
                date=imhead(imagename=file, mode="get", hdkey="date-obs").split('/')

                index=index+1

                #create a new file for different coordinate systems and make a new filename
                newfilename=f'{source}.{str(namefreq)}GHz.{date[0]}.{date[1]}.{date[2]}'
                #mapping from old filename, new filename, and index
                ifile=[index,file,newfilename]
                indexfiles.append(ifile)

                # Route through the image WCS, not the telescope name.  A
                # minority of ALMA FITS products are FK5/J2000.
                if celestial_frame == ICRS:
                    ia.close()
                    ia.open(file)
                    csys = ia.coordsys()
                    apixcoords=csys.topixel(acoords)
                    apixcoord1=[]
                    apixcoord2=[]
                    #input(pixcoords['numeric'][0])
                    apixcoord1.append(apixcoords['numeric'][0])
                    apixcoord2.append(apixcoords['numeric'][1])
                    pixcoord1='NA"'
                    pixcoord2='NA'

                else:
                    ia.close()
                    ia.open(file)
                    csys = ia.coordsys()
                    # J2000/FK5 image, regardless of telescope label.
                    pixcoords=csys.topixel(coords)
                    pixcoord1=[]
                    pixcoord2=[]
                    #input(pixcoords['numeric'][0])
                    pixcoord1.append(pixcoords['numeric'][0])
                    pixcoord2.append(pixcoords['numeric'][1])
                    apixcoord1='NA"'
                    apixcoord2='NA'

                if debugger:
                    print('3')
                filename=file.split('/')
                filename=filename[len(filename)-1]
                force_upper_limit = is_forced_sma_niter0_upper_limit(file, tele, datestr)
                low = force_upper_limit
                ordinary_upper_limit_mjy = None

                source_fit_dir = fit_output_directory(tele, source, file)
                if getchoosercoords:
                    fitfilename=f'images/{flagnum}/chosenfreq/{source}'
                    imfitloc=f'{fitfilename}/{filename}.fit'
                else:
                    fitfilename=source_fit_dir
                    imfitloc=f'{fitfilename}/{filename}_{pipetype}.fit'
                if doublegaussfit:
                    imfitloc_gauss2=f'{fitfilename}/{filename}_noloc.fit'
                    second_gauss_pix=get_multcoords(imfitloc_gauss2, tele, debugger)
                else:
                    imfitloc_gauss2='NA'
                    second_gauss_pix='NA'

                if not os.path.exists(fitfilename):
                    os.makedirs(fitfilename)

                iregion=f'guessregions/{index}region.txt'
                if os.path.exists(iregion):
                    os.remove(iregion)
                Path(iregion).touch()


                #print('\nbeginning imstats\n')
                fixer=masterfixer

                regiontot=f'centerbox[{coords},[{fit_region_width_arcsec}arcsec,{fit_region_width_arcsec}arcsec]]'
                aregiontot=f'centerbox[{acoords},[{fit_region_width_arcsec}arcsec,{fit_region_width_arcsec}arcsec]]'

                worked=write_region(iregion,pixcoord1,pixcoord2,apixcoord1,apixcoord2,tele,file,regiontot,aregiontot,majval,minval,rotval,rotunit,writerreg,fixer,doublegaussfit,second_gauss_pix, debugger, fixtochosenfreq, imagefound, celestial_frame)
                if worked[0]=='bad':
                    badfiles.append([6,file])
                    plt.close()
                    continue
                if debugger:
                    print('4')
                if worked=='skip':
                    continue

                # write_region marks an extended/multiple-component fit with
                # worked[0] == 1; retain its explicit 10-arcsec square.
                if worked[0]==1:
                    fit_region_width_arcsec=square_bmaj_region_width_arcsec(
                        majval, 1, override_arcsec=10
                    )
                    regiontot=f'rotbox[{coords},[{fit_region_width_arcsec}arcsec,{fit_region_width_arcsec}arcsec],{rotval}{rotunit}]'
                    aregiontot=f'rotbox[{acoords},[{fit_region_width_arcsec}arcsec,{fit_region_width_arcsec}arcsec],{rotval}{rotunit}]'

                active_regiontot = aregiontot if celestial_frame == ICRS else regiontot
                active_regionsmall = aregionsmall if celestial_frame == ICRS else regionsmall

                log_file = 'casa_log.txt' 
                casalog.setlogfile(log_file)
                if fit_sizenoloc:
                    try:
                        imstat(file,region=active_regiontot)
                    except Exception as e:
                        badfiles.append([0,file])
                        plt.close()
                        print(f"Error with region {file} ... likely not the location of interest")
                        continue
                breaker=0
                with silence_stderr_fd():
                    try:
                        imfitter=imfit(
                            imagename=file,
                            model=imfitloc,
                            region=active_regiontot,
                            estimates=iregion,
                        )
                    except Exception as e:
                        breaker=1
                if breaker==1:
                    badfiles.append([1,file])
                    plt.close()
                    print(f'imfit does not work for {file} 1')
                if prelimchecker==1:
                    #the following doesn't hold much weight since prelimchecker just throws an error for bad files
                    pos_thresher=0
                    elongthresher=0
                    if breaker!=1:
                        if imfitter['converged'][0]==True:
                            fitaxmaj=imfitter['results']['component0']['shape']['majoraxis']['value']
                            fitaxmajun=imfitter['results']['component0']['shape']['majoraxis']['unit']
                            fitaxmin=imfitter['results']['component0']['shape']['minoraxis']['value']
                            fitaxminun=imfitter['results']['component0']['shape']['minoraxis']['unit']

                            if fitaxmajun=='deg':
                                fitaxmaj=fitaxmaj*60*60
                            if fitaxminun=='deg':
                                fitaxmin=fitaxmin*60*60

                            fitaxang=imfitter['results']['component0']['shape']['positionangle']['value']
                            fitaxangun=imfitter['results']['component0']['shape']['positionangle']['unit']
                            fitaxangerr=imfitter['results']['component0']['shape']['positionangleerror']['value']
                            fitaxangerrun=imfitter['results']['component0']['shape']['positionangleerror']['unit']



                if breaker==1:
                    continue

                conv_test=False
                if imfitter['converged'][0]==True:
                    conv_test=True





                if conv_test==False:
                    if fixtochosenfreq==1 and tele!='SMA':
                        badfiles.append([4,file])
                        breaker=1
                    else:
                        low=True


                if breaker==1:
                    continue

                    '''
                    low=True

                    if fixtochosenfreq==1:
                        multipliertot=1
                    else:
                        multipliertot=10
                    fit_region_width_arcsec=square_bmaj_region_width_arcsec(
                        majval, multipliertot
                    )

                    regiontot=f'rotbox[{coords},[{fit_region_width_arcsec}arcsec,{fit_region_width_arcsec}arcsec],{rotval}{rotunit}]'
                    aregiontot=f'rotbox[{acoords},[{fit_region_width_arcsec}arcsec,{fit_region_width_arcsec}arcsec],{rotval}{rotunit}]'


                    imfitter2='NA'
                    breaker=0
                    if tele=='ALMA':
                        try:
                            imfitter2=imfit(imagename=file,model=imfitloc,region=aregiontot,estimates=iregion)
                        except Exception as e:
                            badfiles.append([2,file])
                            print(f'imfit does not work for {file} 2\n')
                            breaker=1
                    if tele!='ALMA':                  
                        try:
                            imfitter2=imfit(imagename=file,model=imfitloc,region=regiontot,estimates=iregion)
                        except Exception as e:
                            badfiles.append([2,file])
                            print(f'imfit does not work for {file} 2\n')
                            breaker=1

                    if imfitter2['converged'][0]==True: 
                        imfitter=imfitter2
                        conv_test=True

                if breaker==1:
                    continue
                '''

                if conv_test == True:
                    fit_x = imfitter['results']['component0']['pixelcoords'][0]
                    fit_y = imfitter['results']['component0']['pixelcoords'][1]

                    if celestial_frame == ICRS:
                        chosen_x = apixcoord1[0]
                        chosen_y = apixcoord2[0]
                    else:
                        chosen_x = pixcoord1[0]
                        chosen_y = pixcoord2[0]

                    incr = imhead(imagename=file)['incr']
                    dx_deg = (fit_x - chosen_x) * abs(incr[0]) * 180.0 / np.pi
                    dy_deg = (fit_y - chosen_y) * abs(incr[1]) * 180.0 / np.pi
                    sep_deg = np.sqrt(dx_deg**2 + dy_deg**2)

                    if fixtochosenfreq == 1 and imagefound == 1:
                        chosen_image_info = iibest[1][1]
                        if len(chosen_image_info) < 3:
                            raise IndexError(
                                f"Expected chosen image info for {source} to include coordinates and image path, "
                                f"got {chosen_image_info!r}"
                            )
                        chosen_image_path = chosen_image_info[2]
                        chosen_header = imhead(imagename=chosen_image_path)
                        chosen_bmaj_arcsec = chosen_header['restoringbeam']['major']['value']
                        chosen_bmin_arcsec = chosen_header['restoringbeam']['minor']['value']
                        chosen_bmaj_unit = chosen_header['restoringbeam']['major']['unit']
                        chosen_bmin_unit = chosen_header['restoringbeam']['minor']['unit']

                        if chosen_bmaj_unit == 'deg':
                            chosen_bmaj_arcsec = chosen_bmaj_arcsec * 3600.0
                        if chosen_bmin_unit == 'deg':
                            chosen_bmin_arcsec = chosen_bmin_arcsec * 3600.0
                    else:
                        chosen_bmaj_arcsec = majval
                        chosen_bmin_arcsec = minval

                    fit_major_val = imfitter['results']['component0']['shape']['majoraxis']['value']
                    fit_major_unit = imfitter['results']['component0']['shape']['majoraxis']['unit']
                    fit_minor_val = imfitter['results']['component0']['shape']['minoraxis']['value']
                    fit_minor_unit = imfitter['results']['component0']['shape']['minoraxis']['unit']

                    if fit_major_unit == 'deg':
                        fit_major_val = fit_major_val * 3600.0
                    if fit_minor_unit == 'deg':
                        fit_minor_val = fit_minor_val * 3600.0

                    image_path = file

                    append_separation_and_radii(
                            outfile=separation_outfile,
                            source=source,
                            instrument=tele,
                            image_path=image_path,
                            separation_deg=sep_deg,
                            chosen_major_arcsec=chosen_bmaj_arcsec,
                            chosen_minor_arcsec=chosen_bmin_arcsec,
                            interest_major_arcsec=fit_major_val,
                            interest_minor_arcsec=fit_minor_val,
                            chosen_is_limit=False,
                            interest_is_limit=force_upper_limit
                    )

                breaker2=0
                keyword = 'Found no pixels over which to sum for component 0'

                with open(log_file, 'r') as f:
                    log_content = f.read()
                if keyword in log_content:
                    breaker2=1
                os.remove(log_file)
                if breaker2==1:
                    #badfiles.append([8,file])
                    #continue
                    pass
                majorax=f'{majval}{majunit}'
                minorax=f'{minval}{minunit}'
                posnang=f'{rotval}{rotunit}'
                type='orig'
                makedot(active_coords,freq,frequnit,file,imfitter,tele,type)
                makegauss(active_coords,freq,frequnit,file,imfitter,tele,majorax,minorax,posnang)

                #fit a function to the dot to simply get a contour, but color it red
                os.makedirs(source_fit_dir, exist_ok=True)
                imfitloc2=f'{source_fit_dir}/{filename}_g.fit'
                with silence_stderr_fd():
                    imfit(imagename='testg.fits',model=imfitloc2,region=active_regiontot)

                if conv_test==True:
                    writeandstore=1
                    color='b'
                    low,date,ordinary_upper_limit_mjy=printrueconv(imfitter,'NA',writeandstore,color,imfitloc,newfilename,printimages,beamsize,freq,oicoords,frequnit,telescopes, tele, file, source, getchoosercoords, fixtochosenfreq, imagefound, prelimchecker, fit_sizenoloc, doublegaussfit,writenoise, force_upper_limit=force_upper_limit)
                    if not low==True and getchoosercoords==1:
                        anchor_coordinates = getbestrez(
                            allfile,
                            oicoords,
                            tele,
                            telescopes,
                            1,
                            imfitloc,
                        )
                        if anchor_coordinates is None:
                            raise RuntimeError(
                                f'Could not recover the fitted anchor coordinate for {file}'
                            )
                        anchor_frame = image_celestial_frame(imfitloc)
                        icrs_pair = anchor_coordinates[0]
                        j2000_pair = anchor_coordinates[1]
                        native_pair = icrs_pair if anchor_frame == ICRS else j2000_pair
                        portable_input = portable_release_product(file)
                        model_name = f'{tele}.{os.path.basename(file)}.fit'
                        candidate = {
                            'source': source,
                            'telescope': tele,
                            'frequency_ghz': freq,
                            'beam_major_arcsec': majval,
                            'beam_minor_arcsec': minval,
                            'beam_pa_deg': rotval,
                            'beam_geometric_arcsec': beamsize,
                            'anchor_frame': anchor_frame,
                            'anchor_ra_deg': sexagesimal_ra_to_deg(native_pair[0]),
                            'anchor_dec_deg': sexagesimal_dec_to_deg(native_pair[1]),
                            'icrs_ra': icrs_pair[0],
                            'icrs_dec': icrs_pair[1],
                            'j2000_ra': j2000_pair[0],
                            'j2000_dec': j2000_pair[1],
                            'input_product': portable_input,
                            '_input_local_path': file,
                            'model_name': model_name,
                            'model_path': (
                                f'images/{flagnum}/chosenfreq/{source}/{model_name}'
                            ),
                            '_unrenamed_model_path': imfitloc,
                            'run_context': summary_prefix,
                            'selection_rule': CHOSEN_FREQUENCY_SELECTION_RULE,
                            'scientific_disposition': (
                                'rejected_component_seed_source18_fixed_nucleus_rerun_required'
                                if source == 'NGC4945'
                                else 'generic_anchor_seed'
                            ),
                            '_plot_stem': newfilename,
                        }
                        candidate_key = (
                            float(beamsize),
                            portable_input,
                            tele,
                            model_name,
                        )
                        current = choosercoords_by_source.get(source)
                        if current is None or candidate_key < current['_selection_key']:
                            candidate['_selection_key'] = candidate_key
                            choosercoords_by_source[source] = candidate



                    if date=='ONBOUNDRY':
                        if fixtochosenfreq==1:
                            badfiles.append([5,file])
                            plt.close()
                            print(f'ON BOUNDRY: {file}')
                            continue
                if low==True:
                    if fixtochosenfreq and conv_test and not doublegaussfit and not force_upper_limit and ordinary_upper_limit_mjy is None:
                        raise RuntimeError(
                            f'Ordinary fit has no reusable local source-excluded upper limit: {file}'
                        )
                    if fixtochosenfreq==True and tele!='SMA':
                        #scalesize,fov_diag,bdiag
                        scalesize=compute_fov_and_beam(file)
                        if scalesize[1]<10:
                            update_final_noise_audit(
                                file, row_disposition="rejected_field_of_view"
                            )
                            badfiles.append([9,file])
                            plt.close()
                            continue
                        else:
                            upbound=ordinary_upper_limit_mjy
                            fluxval=upbound
                            fluxval=format_number(fluxval)
                            sfluxval=f'{fluxval}*'
                            sbeamsize=format_number(beamsize)
                            color='red'
                            pixlen=100
                            date=makeonlydotcontour(file,imfitloc,sfluxval,sbeamsize, newfilename, printimages,color,pixlen,oicoords,freq,frequnit,imfitter,telescopes,tele, getchoosercoords, fixtochosenfreq,imagefound, prelimchecker, fit_sizenoloc, doublegaussfit)

                    else:                            
                        if force_upper_limit:
                            # Preserve the reviewed dirty-SMA policy: its
                            # limit is three times the one-beam-region RMS,
                            # separate from the ordinary-fit classifier.
                            forced_statistics=image_statistics_mjy_per_beam(
                                file, active_regionsmall
                            )
                            noise=forced_statistics["rms_mjy_per_beam"]
                            upbound=noise*3
                        elif ordinary_upper_limit_mjy is not None:
                            upbound=ordinary_upper_limit_mjy
                        else:
                            noise=image_rms_mjy_per_beam(file, active_regionsmall)
                            upbound=noise*3
                        fluxval=upbound
                        fluxval=format_number(fluxval)
                        sfluxval=f'{fluxval}*'
                        if fixtochosenfreq and force_upper_limit:
                            forced_fit=forced_fit_diagnostic(imfitter)
                            fitted_flux_mjy=None
                            if forced_fit.converged:
                                fitted_flux_mjy=flux_density_to_mjy(
                                    forced_fit.flux_value,
                                    forced_fit.flux_unit,
                                )
                            record_forced_sma_noise_audit(
                                file,
                                forced_statistics,
                                source=source,
                                telescope=tele,
                                frequency_ghz=freq,
                                flux_mjy=fitted_flux_mjy,
                                fit_converged=forced_fit.converged,
                                limit_mjy=upbound,
                                reported_flux_field=sfluxval,
                            )
                        sbeamsize=format_number(beamsize)
                        color='red'
                        pixlen=100
                        date=makedotcontour(file,imfitloc,sfluxval,sbeamsize, newfilename, printimages,color,pixlen,oicoords,freq,frequnit,imfitter,telescopes,tele, getchoosercoords, fixtochosenfreq,imagefound, prelimchecker, fit_sizenoloc, doublegaussfit)

                    if date=='ONBOUNDRY':
                        if fixtochosenfreq:
                            record_reviewed_boundary_rejection(file)
                        print(f'{file} on boundry')
                        continue 
                    #date=makeonlydotcontour(file, sbeamsize, sfluxval, newfilename, printimages)
                    blank='NA'
                    stotarea=f'{beamsize:.5f}'
                    sfreq=f'{freq:.2f}'
                    summary_product=portable_release_product(file)
                    print(f"{source.ljust(15)}{tele.ljust(15)}{sfreq.ljust(19)}{sfluxval.ljust(16)}{date.ljust(15)}{stotarea.ljust(19)}{blank.ljust(15)}{summary_product}")
                    write_fluxerr_row(source, tele, freq, date, None, summary_product)
                    if fixtochosenfreq:
                        update_final_noise_audit(
                            file,
                            row_disposition=(
                                "rejected_source17_requires_source18_join"
                                if source == "NGC4945"
                                else "generic_table_retained"
                            ),
                        )

    if writesummary==True:
        sys.stdout = original_stdoutsum
    fluxerrfile_handle = None


# In[ ]:


if fixtochosenfreq:
    expected_boundary_rejections = set(
        reviewed_boundary_rejections_by_product
    )
    if observed_boundary_rejection_products != expected_boundary_rejections:
        raise RuntimeError(
            "Reviewed boundary-rejection outcome inventory differs from policy: "
            f"missing={sorted(expected_boundary_rejections - observed_boundary_rejection_products)}, "
            f"extra={sorted(observed_boundary_rejection_products - expected_boundary_rejections)}"
        )


#check badfiles
for ibad in badfiles:
    if not 'bigsize' == ibad[0]:
        print(ibad)

# Fixed-to-chosen bad-file diagnostics include several historical reason
# classes. The exact converged half-maximum-contour boundary rejections are
# separately fail-closed by continuum_boundary_rejection_policy.csv.


# In[ ]:


# make badfiles txts and images
from casatasks import exportfits
import matplotlib

write=True


# make badfiles txts and images
if len(badfiles)!=0:
    if printimages==True:
        import matplotlib.pyplot as plt
    else:
        import matplotlib
        matplotlib.use('Agg')      # must come before pyplot is imported
        import matplotlib.pyplot as plt
    #check bad files      


    if prelimchecker==1:
        badfileversion='prelim'

    if fit_sizenoloc==1:
        badfileversion='noloc'

    if getchoosercoords==1:
        badfileversion='getchooser'
        if len(badfiles)!=0:
            #input('bad files for getchoosercoords')
            pass

    if fixtochosenfreq==1:
        badfileversion='fixedtochosen'
        if len(badfiles)!=0:
            #input('bad files for getchoosercoords')
            pass

    if os.path.exists(f'bad/badfiles/{badfileversion}'):
        shutil.rmtree(f'bad/badfiles/{badfileversion}')
    os.makedirs(f'bad/badfiles/{badfileversion}')
    if os.path.exists((f'bad/badimages/bad_location/{badfileversion}')):
        shutil.rmtree((f'bad/badimages/bad_location/{badfileversion}'))
    if os.path.exists((f'bad/badimages/bad_obs/{badfileversion}')):
        shutil.rmtree((f'bad/badimages/bad_obs/{badfileversion}'))
    if os.path.exists((f'bad_fixed/{badfileversion}')):
        shutil.rmtree((f'bad_fixed/{badfileversion}'))


    badfiletxt=f'bad/badfiles/{badfileversion}/badfiles_{badfileversion}.txt'


    with open(badfiletxt, 'w') as sumfile:
        if write==True:
            original_stdout = sys.stdout
            sys.stdout = sumfile
        for ibad in badfiles:
            output_badfilename='NA'
            ibad1val=ibad[1]
            onlyname=str(ibad1val).split('/')[len(ibad1val.split('/'))-1]
            printer=False



            if ibad[0]==0:
                getimview(ibad1val)
                print(ibad1val)
                if not os.path.exists(f'bad/badimages/bad_location/{badfileversion}'):
                    os.makedirs(f'bad/badimages/bad_location/{badfileversion}')  
                output_badfilename=f'bad/badimages/bad_location/{badfileversion}/badimages_location_{onlyname}.pdf'
                printer=True

            if ibad[0]!=0:
                if prelimchecker==1:
                    scalesize=compute_fov_and_beam(ibad1val)
                    scalesizeval=scalesize[0]
                    if scalesizeval<60:
                        #many files seem just fine... just print them off and record which ones are bad
                        getimview(ibad1val)
                        print(ibad1val)
                        if not os.path.exists(f'bad/badimages/bad_obs/{badfileversion}'):
                            os.makedirs(f'bad/badimages/bad_obs/{badfileversion}')  
                        output_badfilename=f'bad/badimages/bad_obs/{badfileversion}/badimages_beam_FOV_{onlyname}.pdf'
                        printer=True

                    if scalesizeval>=60:
                        #prelim checker has no prior so any large image will be flagged... so 
                        #these files are likely just fine
                        if printbadfix==1:
                            getimview(ibad1val)
                            if not os.path.exists(f'bad_fixed/{badfileversion}'):
                                os.makedirs(f'bad_fixed/{badfileversion}')  
                            output_badfilename=f'bad_fixed/{badfileversion}/large_beam_FOV_{onlyname}.pdf'
                        else:
                            printer=False
                else:
                    getimview(ibad1val)
                    print(ibad1val)
                    if not os.path.exists(f'bad/badimages/bad_obs/{badfileversion}'):
                        os.makedirs(f'bad/badimages/bad_obs/{badfileversion}')  
                    output_badfilename=f'bad/badimages/bad_obs/{badfileversion}/badimages_beam_FOV_{onlyname}.pdf'
                    printer=True

            if printer==True:
                # Open the FITS file (the input image)
                if '/Gaussian_Fit_Pipeline/230316_128/images/flag2_cal1/fields/science/' in ibad1val:
                    sourcename=onlyname.split('_test.image')[0]
                    badsmadir=f'analysis/bad_sma_fits/{sourcename}'
                    if not os.path.exists(badsmadir):
                        os.makedirs(badsmadir)
                    fitsname=f'{sourcename}.fits'
                    newsmabadfit_loc=os.path.join(badsmadir,fitsname)
                    exportfits(imagename=ibad1val, fitsimage=fitsname, overwrite=True)
                    ibad1val=fitsname


                hdul_main = fits.open(ibad1val)
                wcs_main = WCS(hdul_main[0].header, naxis=2)

                # Set up WCS-based plot
                fig = plt.figure(figsize=(8, 6))
                ax = plt.subplot(111, projection=wcs_main)

                # Plot the image data
                ax.imshow(hdul_main[0].data.squeeze(), origin='lower', cmap='viridis')

                # Save and show
                plt.savefig(output_badfilename, dpi=100, bbox_inches='tight')

            plt.close()
        if write==True:
            sys.stdout = original_stdout


# In[ ]:


#set bad files to master bad list
write=True

if not os.path.exists(f'bad/masterbadfiles'):
    os.makedirs(f'bad/masterbadfiles')  

#for fit_sizenoloc
#all files from badfiletxt
#these runtime fit failures will be added to the masterbadfiles folder;
#reviewed exclusions are handled separately by the CSV policy above
if fit_sizenoloc==1:
    badfileversion='noloc'
    badfiletxt=f'bad/badfiles/{badfileversion}/badfiles_{badfileversion}.txt'

    #check band files
    masterbadfiletxt=f'bad/masterbadfiles/masterbadfiles_{badfileversion}.txt'

    shutil.copy(badfiletxt, masterbadfiletxt)




# In[ ]:


#adjust fitsummaries to skip error files
import os
import sys

textin=f'fitsumfiles/witherrors/{flagnum}_fitsummary_{pipetype}_witherrors.txt'
if not os.path.exists('fitsumfiles'):
    os.makedirs('fitsumfiles/noerrors')
textin=f'fitsumfiles/witherrors/{summary_prefix}_fitsummary_{pipetype}_witherrors.txt'
textout=f'fitsumfiles/noerrors/{summary_prefix}_fitsummary_{pipetype}_noerrors.txt'

if os.path.exists(textout):
    os.remove(textout)

os.makedirs('fitsumfiles/noerrors', exist_ok=True)

write=True
#write=False
with open(textin, 'r') as infile:
    with open(textout, 'w') as outfile:
        for line in infile:
            if 'Error' in line:
                continue
            if 'pixels' in line:
                continue
            if 'did not converge' in line:
                continue
            if 'on boundry' in line.lower():
                continue
            if 'imfit does not work' in line:
                continue
            if len(line)==1:
                continue
            original_stdout = sys.stdout
            sys.stdout = outfile
            print(line.strip())
            sys.stdout = original_stdout


# In[ ]:

import numpy as np
import sys
import os

all_lines=[]
textin=f'fitsumfiles/noerrors/{summary_prefix}_fitsummary_{pipetype}_noerrors.txt'

if not os.path.exists('fitsumfiles'):
    os.makedirs('fitsumfiles/sorted')
textout=f'fitsumfiles/sorted/{summary_prefix}_fitsummary_{pipetype}_sorted.txt'

os.makedirs('fitsumfiles/sorted', exist_ok=True)

write=True
#write=False

totex=0
with open(textout, 'w') as file:
    if write==True:
        original_stdout = sys.stdout
        sys.stdout = file
    all_lines=[]
    with open(textin, 'r') as infile:
        tline=[]


        for line in infile:
            line = line.split(' ')
            line2 = []
            for i in line:
                if i != '':
                    line2.append(i)

            if len(line2) < 6:
                continue

            all_lines.append(line2)



    source_order=[]
    doneindex=[]
    n = len(all_lines)
    for i in range(n):
        if i==0:
            continue
        breaker=False
        for idoneindex in doneindex:
            if i==idoneindex:
                breaker=True
        if breaker==True:
            continue


        isource_order=[]
        isource_order.append(all_lines[i])
        for j in range(n):
            if j==1:
                continue
            if j!=i:
                if all_lines[i][0]==all_lines[j][0]:
                    #so that while iterating through i the source is not repeated
                    isource_order.append(all_lines[j])
                    doneindex.append(j)
            n2 = len(isource_order)


        isource_order = sorted(isource_order, key=lambda x: (
            float(x[2]),
            calendar_date_sort_key(x[4]),
            (-1 if x[5] == 'NA' else float(x[5])),
            portable_release_product(x[7]),
        ))


        for i2source_order in isource_order:
            print(i2source_order)
        print('\n')
    if write==True:
        sys.stdout = original_stdout


# In[ ]:


# Keep one deterministic minimum-beam candidate per source, then materialize
# the exact CASA model trees needed by the final pass.
if getchoosercoords==1:
    for source in sorted(choosercoords_by_source):
        ichoose = choosercoords_by_source[source]
        base = os.path.basename(ichoose['_input_local_path'])
        tele=ichoose['telescope']
        
        fit_filename = f'images/{flagnum}/chosenfreq/{source}/{base}.fit'
        files = os.listdir(f'images/{flagnum}/chosenfreq/{source}')

        plot_filename_out = f"images/{flagnum}/chosenfreq/{source}/{tele}.{ichoose['_plot_stem']}.pdf"
        fit_filename_out = f'images/{flagnum}/chosenfreq/{source}/{tele}.{base}.fit'

        counter=0
        for file in files:
          if ichoose['_plot_stem'] in file:
              counter=counter+1
              plot_filename = f'images/{flagnum}/chosenfreq/{source}/{file}'
          if file in fit_filename:
              counter=counter+1
        if counter!=2:
           raise RuntimeError(
               f'Expected one plot and one CASA model for {source}, found {counter}'
           )
        else:

            shutil.move(plot_filename,plot_filename_out)
            shutil.move(fit_filename,fit_filename_out)
        files = os.listdir(f'images/{flagnum}/chosenfreq/{source}')
        for file in files:
          if 'pdf' in file:
             pdfin = f'images/{flagnum}/chosenfreq/{source}/{file}'
             pdfout = f'images/{flagnum}/chosenfreq/{file}'
             shutil.move(pdfin,pdfout)
          elif not file in fit_filename_out:
             shutil.rmtree(f'images/{flagnum}/chosenfreq/{source}/{file}')
             pass
          else:
              pass

chosenfreq_dir = f'images/{flagnum}/chosenfreq'
remove_unselected_chosen_frequency_source_folders(chosenfreq_dir)


# In[ ]:


if getchoosercoords==1:
    candidate_rows = []
    for source in sorted(choosercoords_by_source):
        candidate = dict(choosercoords_by_source[source])
        input_path = Path(candidate.pop('_input_local_path'))
        candidate.pop('_unrenamed_model_path')
        candidate.pop('_selection_key')
        candidate.pop('_plot_stem')
        model_path = Path(candidate['model_path'])
        if not model_path.is_dir():
            raise FileNotFoundError(
                f'Chosen-frequency CASA model was not retained for {source}: {model_path}'
            )
        candidate['input_sha256'] = path_sha256(input_path)
        candidate['model_tree_sha256'] = tree_sha256(model_path)
        candidate_rows.append(candidate)

    candidate_manifest = Path(
        f'chosen_frequency/candidates/{summary_prefix}_getchooser.csv'
    )
    write_chosen_frequency_manifest(candidate_manifest, candidate_rows)
    consolidate_chosen_frequency_anchors(
        [candidate_manifest],
        analysis_root=Path('.'),
        output_manifest=Path('chosen_frequency/chosen_frequency_anchors.csv'),
        models_output=Path('images/chosenfreq'),
        component_policy=Path('chosen_frequency_component_policy.csv'),
        replace_models=True,
    )


# In[ ]:
if moveniter:
        import os
        import shutil

        def merge_tree(src, dst):
                    os.makedirs(dst, exist_ok=True)
                    for item in os.listdir(src):
                                src_item = os.path.join(src, item)
                                dst_item = os.path.join(dst, item)
                                if os.path.isdir(src_item):
                                            merge_tree(src_item, dst_item)
                                else:
                                            if os.path.exists(dst_item):
                                                        os.remove(dst_item)
                                            shutil.move(src_item, dst_item)
                    if os.path.isdir(src) and not os.listdir(src):
                                os.rmdir(src)

        niter_folder = str(Path(__file__).resolve().parents[1] / "sma_calibration" / f"nitercal_{niternum}")
        os.makedirs(niter_folder, exist_ok=True)

        for name in ["fitsumfiles", "bad", "separations"]:
                    if os.path.exists(name):
                                dst = os.path.join(niter_folder, name)
                                if os.path.isdir(name):
                                            merge_tree(name, dst)
                                else:
                                            if os.path.exists(dst):
                                                        os.remove(dst)
                                            shutil.move(name, dst)

        if os.path.exists("images"):
                    dst_images = os.path.join(niter_folder, "images")
                    os.makedirs(dst_images, exist_ok=True)

                    for item in os.listdir("images"):
                                if item == "chosenfreq":
                                            continue

                                src = os.path.join("images", item)
                                dst = os.path.join(dst_images, item)

                                if os.path.exists(dst):
                                            if os.path.isdir(dst):
                                                        shutil.rmtree(dst)
                                            else:
                                                        os.remove(dst)

                                shutil.move(src, dst)
