#!/usr/bin/env python
# coding: utf-8

import csv
import os
from pathlib import Path

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
parser.add_argument("--moveniter", type=int, default=0)
parser.add_argument("--niternum", type=int, default=0)
parser.add_argument("--datestr", type=str, default="")

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

def normalize_source_name(source_name, killstrings=None):
    if source_name is None:
        return ""

    cleaned_source = str(source_name).strip()
    if len(cleaned_source.split('.')) > 0:
        cleaned_source = cleaned_source.split('.')[0]

    if killstrings is None:
        killstrings = globals().get("listofkills", ["Galaxy"])

    for killstring in killstrings:
        if killstring in cleaned_source:
            cleaned_source = cleaned_source.replace(killstring, '')

    if cleaned_source:
        last_char = cleaned_source[-1]
        if last_char.isalpha():
            cleaned_source = cleaned_source[:-1] + last_char.upper()

    return cleaned_source

def load_selected_clean_manifest(manifest_path, date_code="", flag_number=None):
    selected_lookup = {}
    selected_original_sources = []

    if not manifest_path or not os.path.exists(manifest_path):
        return selected_lookup, selected_original_sources

    with open(manifest_path, newline="") as manifest_file:
        reader = csv.DictReader(manifest_file)
        for row in reader:
            image_path = (row.get("image_path") or "").strip()
            if not image_path:
                continue

            row_date = (row.get("date_code") or "").strip()
            if date_code and row_date and row_date != date_code:
                continue

            row_flagnum = (row.get("flagnum") or "").strip()
            if flag_number is not None and row_flagnum:
                try:
                    if int(float(row_flagnum)) != int(flag_number):
                        continue
                except ValueError:
                    continue

            source_name = (row.get("source") or "").strip()
            if not source_name:
                source_name = os.path.basename(image_path).replace("_test.image", "")

            normalized_source = normalize_source_name(source_name)
            if not normalized_source:
                continue

            selected_lookup[normalized_source] = image_path
            selected_original_sources.append(
                os.path.basename(image_path).replace("_test.image", "")
            )

    return selected_lookup, list(dict.fromkeys(selected_original_sources))

selected_clean_image_lookup, selected_clean_original_sources = load_selected_clean_manifest(
    selected_clean_manifest,
    datestr,
    flagnum,
)

def get_selected_sma_image_path(source):
    normalized_source = normalize_source_name(source)
    image_path = selected_clean_image_lookup.get(normalized_source)
    if image_path and os.path.exists(image_path):
        return image_path
    return None

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


def getbestrez(getfiles,icoords,tele,alltele,choosefile,choicefile):
    while 0==0:
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
            if besttele=='ALMA':
                if tele=='NA':
                    aicoords=[f"{imstat(bestrezfile,region=aregion)['maxposf'].split(',')[0]}", f"{imstat(bestrezfile,region=aregion)['maxposf'].split(',')[1].replace('.',':',2)}"]
                if tele!='NA':
                    aicoords=[f"{imstat(bestrezfile)['maxposf'].split(',')[0]}", f"{imstat(bestrezfile)['maxposf'].split(',')[1].replace('.',':',2)}"]
                aicoord1=aicoords[0].split(':')
                aselectcoord1=f'{aicoord1[0]}h{aicoord1[1]}m{aicoord1[2]}s'
                aicoord2=aicoords[1].split(':')
                aselectcoord2=f'{aicoord2[0]}d{aicoord2[1]}m{aicoord2[2]}s'

                icoords=au.ICRSToJ2000(f'{aicoords[0]} {aicoords[1]}',verbose=False).split(', ')
                icoord1=icoords[0].split(':')
                selectcoord1=f'{icoord1[0]}h{icoord1[1]}m{icoord1[2]}s'
                icoord2=icoords[1].split(':')
                selectcoord2=f'{icoord2[0]}d{icoord2[1]}m{icoord2[2]}s'
                break
            if besttele!='ALMA':
                if tele=='NA':
                    icoords=[f"{imstat(bestrezfile,region=region)['maxposf'].split(',')[0]}", f"{imstat(bestrezfile,region=region)['maxposf'].split(',')[1].replace('.',':',2)}"]
                if tele!='NA':
                    icoords=[f"{imstat(bestrezfile)['maxposf'].split(',')[0]}", f"{imstat(bestrezfile)['maxposf'].split(',')[1].replace('.',':',2)}"]
                icoord1=icoords[0].split(':')
                selectcoord1=f'{icoord1[0]}h{icoord1[1]}m{icoord1[2]}s'   
                icoord2=icoords[1].split(':')
                selectcoord2=f'{icoord2[0]}d{icoord2[1]}m{icoord2[2]}s'

                aicoords=au.J2000ToICRS(f'{icoords[0]} {icoords[1]}',verbose=False).split(', ')
                aicoord1=aicoords[0].split(':')
                aselectcoord1=f'{aicoord1[0]}h{aicoord1[1]}m{aicoord1[2]}s'
                aicoord2=aicoords[1].split(':')
                aselectcoord2=f'{aicoord2[0]}d{aicoord2[1]}m{aicoord2[2]}s'
                break
        except Exception as e:
            i=-1
            braker=-1
            if choosefile==0:
                for telescopes in getfiles:
                    i=i+1
                    #iterate past the sourcename
                    if i==0:
                        continue
                    k=-1
                    for telescope in telescopes:
                        k=k+1
                        j=-1
                        for file in telescope:
                            j=j+1
                            if file==bestrezfile:
                                getfiles[i][k][j]=None
                                braker=1
                                break
                        if braker==1:
                            break
                    if braker==1:
                        break


    sendback=[aselectcoord1,aselectcoord2],[selectcoord1,aselectcoord2],bestrezfile
    return(sendback)




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
def write_region(iregion,pixcoord1,pixcoord2,apixcoord1,apixcoord2,tele,file,region,aregion,majval,minval,rotval,rotunit,write,ofixer,doublegaussfit,second_gauss_pix, debugger, fixtochosenfreq, imagefound):
    worked1=0
    with open(iregion, 'w') as file2:
        if write==1:
            original_stdout2 = sys.stdout
            sys.stdout = file2

        if tele!='ALMA':
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

        if tele=='ALMA':
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
    #There doesn't appear to be a difference if the direction is J2000 or ICRS
    if tele=='ALMA':
        direction = f"ICRS {ra} {dec}"
    if tele!='ALMA':
        direction = f"J2000 {ra} {dec}"
    cl.done()

    cl.addcomponent(dir=direction, flux=1.0, fluxunit='Jy', freq=f'{freq}{frequn}', shape="Gaussian", majoraxis=majorax, minoraxis=minorax, positionangle=posang)

    aspp=f'{str(aspp)}arcsec'
    if os.path.exists("testg.im"):
        shutil.rmtree("testg.im")
    ia.fromshape("testg.im",[fulpix1,fulpix2],overwrite=True)
    cs=ia.coordsys()
    cs.setunits(['rad','rad'])
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
    #There doesn't appear to be a difference if the direction is J2000 or ICRS
    if tele=='ALMA':
        direction = f"ICRS {ra} {dec}"
    if tele!='ALMA':
        direction = f"J2000 {ra} {dec}"
    cl.done()

    cl.addcomponent(dir=direction, flux=1.0, fluxunit='Jy', freq=f'{freq}{frequn}', shape="Point")
    aspp=f'{str(aspp)}arcsec'
    if os.path.exists(f"{pointname}.im"):
        shutil.rmtree(f"{pointname}.im")
        os.mkdir(f"{pointname}.im")
    ia.fromshape(f"{pointname}.im",[fulpix1,fulpix2],overwrite=True)
    cs=ia.coordsys()
    cs.setunits(['rad','rad'])
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
        data_main *= 1000  # Multiply intensity by 1000 to get into mJy
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
            if tele=='ALMA':
                makedot(acoords,freq2,frequnit,file,imfitter,tele,type)
            if tele!='ALMA':
                makedot(coords,freq2,frequnit,file,imfitter,tele,type)

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
            if tele=='ALMA':
                makedot(acoords,freq2,frequnit,file,imfitter,tele,type)
            if tele!='ALMA':
                makedot(coords,freq2,frequnit,file,imfitter,tele,type)

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
        data_main *= 1000  # Multiply intensity by 1000 to get into mJy
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
        # Case 1: If the number is less than 1 (e.g., 0.004567)
        if number < 0:
                # Format the positive part of the number and prepend the negative sign
                return "-" + format_number(abs(number))
        if number < 1:
                # Convert the number to a string with high precision
                num_str = f"{number:.16g}"

                # Identify the leading zeros and decimal point
                leading_part = []
                for char in num_str:
                        if char == '0' or char == '.':
                                leading_part.append(char)
                        else:
                                break

                # Remove leading zeros and the decimal point for significant digits
                significant_digits = "".join(char for char in num_str if char.isdigit() and char != "0")

                # Extract the first two non-zero digits
                if len(significant_digits) >= 2:
                        first_two = significant_digits[:2]
                elif len(significant_digits) == 1:
                        first_two = significant_digits + "0"  # Pad with zero if only one significant digit exists
                else:
                        first_two = "00"  # Handle edge case like 0

                # Combine the leading zeros and the first two non-zero digits
                return "".join(leading_part) + first_two

        # Case 2: If the number is greater than or equal to 1
        else:
                # Format to two decimal places
                return f"{number:.2f}"


# In[ ]:


#print statements and make picturesmakedotcontour
def printrueconv(imfitter,date,make,color,imfitloc,newfilename,printimages,beamsize,freq,oicoords,frequnit,telescopes, tele, file, source, getchoosercoords, fixtochosenfreq, imagefound, prelimchecker, fit_sizenoloc, doublegaussfit,writenoise, force_upper_limit=False):
    if doublegaussfit:
        low=[]
        detectregions=[]
    else:
        low=False
    nelem=imfitter['deconvolved']['nelements']
    counter=-1
    for inelem in range(nelem):
        if doublegaussfit:
            ilow=False
        counter=counter+1
        component=f'component{inelem}'
        fluxerr=imfitter['deconvolved'][component]['flux']['error'][0]
        fluxval=imfitter['deconvolved'][component]['flux']['value'][0]
        fluxval=fluxval*1000
        fluxun=imfitter['deconvolved'][component]['flux']['unit']

        #FWHM to sigma from wikipedia
        factorto3sig= 3 / (2 * np.sqrt(2 * np.log(2)))

        fitmajax=imfitter['results'][component]['shape']['majoraxis']['value']
        fitmajax=fitmajax*factorto3sig
        fitmajunit=imfitter['results'][component]['shape']['majoraxis']['unit']
        if fitmajunit=='deg':
            fitmajax=fitmajax*60*60


        fitminjax=imfitter['results'][component]['shape']['minoraxis']['value']
        fitminjax=fitminjax*factorto3sig
        fitminjunit=imfitter['results'][component]['shape']['minoraxis']['unit']
        if fitminjunit=='deg':
            noisefitminjax=fitminjax*60*60

        fitpos=imfitter['results'][component]['shape']['positionangle']['value']
        fitposun=imfitter['results'][component]['shape']['positionangle']['unit']
        if fitposun!='deg':
            input('pos not in deg')

        totstats = imstat(imagename=file)
        blc = totstats['blc']  # bottom-left corner: [x, y]
        trc = totstats['trc']  # top-right corner: [x, y]
        totregion = f'box[[{blc[0]}pix, {blc[1]}pix], [{trc[0]}pix, {trc[1]}pix]]'

        gausscen1=imfitter['results'][component]['pixelcoords'][0]
        gausscen2=imfitter['results'][component]['pixelcoords'][1]

        if not doublegaussfit:
            try:
                detectregion=f'ellipse[[{gausscen1}pix, {gausscen2}pix], [{fitmajax}arcsec, {fitminjax}arcsec], {fitpos}deg]'
                #noiseregion=f'{totregion}-{detectregion}'
                noiseregfile='noiseregfile.crtf'
                if os.path.exists(noiseregfile):
                    os.remove(noiseregfile)
                with open(noiseregfile, 'w') as noiseregfile:
                    if writenoise==True:
                        original_stdoutnoisereg = sys.stdout
                        sys.stdout = noiseregfile
                    print('#CRTFv0')
                    print(totregion)
                    print(f'-{detectregion}')
                    if writenoise==True:
                        sys.stdout = original_stdoutnoisereg
                #noisesnr=fluxval/imstat(imagename=file,region='noiseregfile.crtf')['rms'][0]
                noiseval=imstat(imagename=file,region='noiseregfile.crtf')['rms'][0]*1000
                if math.isnan(fluxerr):
                    tot_err=noiseval
                else:
                    fluxerr=fluxerr*1000
                    tot_err=np.sqrt(fluxerr**2 + noiseval**2)
                noisesnr=fluxval/tot_err

            except Exception as e:
                print(f'cannot make noise region for {file}:')
                print(e)

            if noisesnr<3:
                low=True

            # Treat fallback SMA nitercal_0 images as upper limits even when the
            # fitted SNR is formally above threshold.
            if force_upper_limit:
                low=True

            fluxval=format_number(fluxval)
            sfluxval=fluxval
            sbeamsize=format_number(beamsize)

            if low!=True:
                if counter==0:
                    if make==True:
                        pixlen=100
                        date=makedotcontour(file,imfitloc,sfluxval,sbeamsize, newfilename, printimages,color,pixlen,oicoords,freq,frequnit,imfitter,telescopes,tele, getchoosercoords, fixtochosenfreq,imagefound, prelimchecker, fit_sizenoloc, doublegaussfit)
                        if date==None:
                            return(date,date)
                        if date=='ONBOUNDRY':
                            return(date,date)
                        stotarea=f'{beamsize:.5f}'
                        sfreq=f'{freq:.2f}'
                        ssnr=str(f'{noisesnr:.3f}')
                        if color=='r':
                            sfluxval=fluxval+'*'
                        print(f"{source.ljust(15)}{tele.ljust(15)}{sfreq.ljust(19)}{sfluxval.ljust(16)}{date.ljust(15)}{stotarea.ljust(19)}{ssnr.ljust(15)}{file}")
                        write_fluxerr_row(source, tele, freq, date, fluxerr, file)
                        if make!=True:
                            date='NA'
                        return(low,date)
                    else:
                        input('multcoords function being used but shouldnt be')
            else:
                return(low,date)
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

        noiseval=imstat(imagename=file,region='noiseregfile.crtf')['rms'][0]*1000

        colors=[]
        counter=-1
        sfluxarray=[]
        for inelem in range(nelem):
            counter=counter+1
            component=f'component{inelem}'
            fluxerr=imfitter['deconvolved'][component]['flux']['error'][0]
            fluxval=imfitter['deconvolved'][component]['flux']['value'][0]
            fluxun=imfitter['deconvolved'][component]['flux']['unit']
            tot_err=np.sqrt(fluxerr**2 + noiseval**2)
            noisesnr=fluxval/tot_err
            tot_err=np.sqrt(fluxerr**2 + noiseval**2)
            noisesnr=fluxval/tot_err
            if noisesnr<3:
                low=True
                print(f'bad double fit for {file}')
                return(low,date)
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
                fluxval=fluxval*1000
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
                print(f"{source.ljust(15)}{scounter.ljust(5)}{tele.ljust(15)}{sfreq.ljust(19)}{sfluxval.ljust(16)}{date.ljust(15)}{stotarea.ljust(19)}{ssnr.ljust(15)}{file}")
                write_fluxerr_row(source, tele, freq, date, fluxerr * 1000, file, component=counter)

        if make==True:
            pixlen=100
            date=makedotcontour(file,imfitloc,sfluxarray,sbeamsize, newfilename, printimages,colors,pixlen,oicoords,freq,frequnit,imfitter,telescopes,tele, getchoosercoords, fixtochosenfreq,imagefound, prelimchecker, fit_sizenoloc, doublegaussfit)
            if date==None:
                return(date,date)
            if date=='ONBOUNDRY':
                return(date,date)
        if make!=True:
            date='NA'


        return(low,date)


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
choosercoordsarray=[]
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
#this sets the region total for the double guassian fitted image
dubmultval=3

###
###
###
###
debugger=0
#debugger=1

yestext=True
#yestext=False

filechooser=False
#filechooser=True
filechoice='VLA/NGC3079/43.3I0.04_AG0527_1998FEB25_1_346.U26.6S.imfits'

biggerregion='VLA/NGC3079/43.3I0.04_AG0527_1998FEB25_1_346.U26.6S.imfits'

lowfiles=['VLA/NGC4258/15.3I0.14_BR0043_1997JAN04_1_351.U42.7S.imfits', 'VLA/NGC4258/22.4I0.07_AH0594_1996DEC31_1_294.U42.8S.imfits']

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

#noise region for double gaussian fits
writerreglowdub=1
#writerreglowdub=0

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
    masterfixerlow='nofix'
    pipetype='prelim'

#fix the beam size (for all sources)
if fit_sizenoloc==1:
    masterfixer=0
    masterfixerlow=0
    pipetype='noloc'

#fix the beam size (chosen freq)
if getchoosercoords==1:
    masterfixer=0
    masterfixerlow=0
    pipetype='getchooser'

#fix the location and beam size
if fixtochosenfreq==1:
    masterfixer=0
    masterfixerlow=0
    pipetype='fixedtochosen'

#fix the location and beam size
if doublegaussfit==1:
    masterfixer=1
    masterfixerlow=1
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

        nonames=1
        with open(source_text, 'r') as infile:
            for params in infile:
                params=eval(params)
                for inamecoords in params:
                    oname=inamecoords[0][0]
                    nname=inamecoords[1][0][0]
                    ra=inamecoords[1][1]
                    dec=inamecoords[1][2]
                    if not source in oname:
                        nonames=0
                        continue
                    icoords=[]
                    icoords.append(ra)
                    icoords.append(dec)
                    oicoords=icoords.copy()
        if nonames==1:
            input(f'no names for {source}')
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
                #sendback=[aselectcoord1,aselectcoord2],[selectcoord1,aselectcoord2],bestrezfile
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
        iibest=='NA'
        for ibest in bestcoords:
            if ibest[0][0]==source:
            # if the first entry is the source (the NED query location)
                iibest=ibest
                i=i+1

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

                #for a single gaussian fit
                #consider a region "multiplier" times the size of the beam for a compact source
                #minval and majval which are FWHM (diamaters)
                if fixtochosenfreq==1 and imagefound==1:
                    multipliertot=fixmultval
                elif fixtochosenfreq==1 and imagefound==0:
                    multipliertot=nofixmultval*2
                else:
                    multipliertot=nofixmultval

                offset1tot=majval*multipliertot
                if fixtochosenfreq==1 and imagefound==0:
                    offset2tot=offset1tot
                else:
                    offset2tot=minval*multipliertot

                if fixtochosenfreq==1 and imagefound==0:
                    if file==biggerregion:
                        offset1tot=3
                        offset2tot=3
                    if offset2tot>3:
                        offset1tot=3
                        offset2tot=3 
                    if offset1tot>3:
                        offset1tot=3
                        offset2tot=3 

                if getchoosercoords==1 and offset1tot<1:
                    offset1tot=1
                if getchoosercoords==1 and offset2tot<1:
                    offset2tot=1

                if getchoosercoords==1:
                    #nofixmultval must be adjust for certain sources
                    for inofix in nofixmultval_exc:
                        if inofix[0].lower() in source.lower():
                            offset1tot=inofix[1]


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

                #for ALMA
                if tele=='ALMA':
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

                #for others
                else:
                    ia.close()
                    ia.open(file)
                    csys = ia.coordsys()
                    #for VLA and SMA using original imagefile
                    ia.close()
                    ia.open(file)
                    csys = ia.coordsys()
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

                for oi in oallsources:
                    if source in oi:
                        if getchoosercoords:
                            fitfilename=f'images/{flagnum}/chosenfreq/{source}'
                            imfitloc=f'images/{flagnum}/chosenfreq/{source}/{filename}.fit' 
                        else:
                            fitfilename=f'{tele}/{oi}/fitfiles'
                            imfitloc=f'{fitfilename}/{filename}_{pipetype}.fit'
                        if doublegaussfit:
                            imfitloc_gauss2=f'{fitfilename}/{filename}_noloc.fit'
                            #output will be in same coordinate system as the current file... since it is a fit from the same telescope
                            if debugger:
                                print('3.0.0')
                            second_gauss_pix=get_multcoords(imfitloc_gauss2, tele, debugger)
                            if debugger:
                                print('3.0.1')
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

                regiontot=f'centerbox[{coords},[{offset1tot}arcsec,{offset1tot}arcsec]]'
                aregiontot=f'centerbox[{acoords},[{offset1tot}arcsec,{offset1tot}arcsec]]'

                worked=write_region(iregion,pixcoord1,pixcoord2,apixcoord1,apixcoord2,tele,file,regiontot,aregiontot,majval,minval,rotval,rotunit,writerreg,fixer,doublegaussfit,second_gauss_pix, debugger, fixtochosenfreq, imagefound)
                if worked[0]=='bad':
                    badfiles.append([6,file])
                    plt.close()
                    continue
                if debugger:
                    print('4')
                if worked=='skip':
                    continue

                #for extended sources 10" should be used
                #otherwise a region twice the beam size should suffice (compact)
                #for a multiple gaussian fit
                if worked[0]==1:
                    offset1tot=10
                    offset2tot=10
                    regiontot=f'rotbox[{coords},[{offset1tot}arcsec,{offset1tot}arcsec],{rotval}{rotunit}]'
                    aregiontot=f'rotbox[{acoords},[{offset1tot}arcsec,{offset1tot}arcsec],{rotval}{rotunit}]'

                log_file = 'casa_log.txt' 
                casalog.setlogfile(log_file)
                if fit_sizenoloc:
                    try:
                        imstat(file,region=regiontot)
                    except Exception as e:
                        badfiles.append([0,file])
                        plt.close()
                        print(f"Error with region {file} ... likely not the location of interest")
                        continue
                breaker=0
                if tele=='ALMA':
                    with silence_stderr_fd():
                            try:
                                imfitter=imfit(imagename=file,model=imfitloc,region=aregiontot,estimates=iregion)
                            except Exception as e:
                                breaker=1
                    if breaker==1:
                        badfiles.append([1,file])
                        plt.close()
                        print(f'imfit does not work for {file} 1')
                if tele!='ALMA':       
                    with silence_stderr_fd():        
                          try:
                              imfitter=imfit(imagename=file,model=imfitloc,region=regiontot,estimates=iregion)
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


                if fixtochosenfreq==1:
                    for ilow in lowfiles:
                        if ilow==file:  
                          low=True
                          conv_test=False

                if breaker==1:
                    continue

                    '''
                    low=True

                    if fixtochosenfreq==1:
                        multipliertot=1
                    else:
                        multipliertot=10
                    offset1tot=majval*multipliertot
                    offset2tot=minval*multipliertot

                    regiontot=f'rotbox[{coords},[{offset1tot}arcsec,{offset1tot}arcsec],{rotval}{rotunit}]'
                    aregiontot=f'rotbox[{acoords},[{offset1tot}arcsec,{offset1tot}arcsec],{rotval}{rotunit}]'


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

                    if tele == 'ALMA':
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
                if tele=='ALMA':
                    #Here you can see the difference if you change coords to acoords
                    type='orig'
                    makedot(acoords,freq,frequnit,file,imfitter,tele,type)
                    makegauss(acoords,freq,frequnit,file,imfitter,tele,majorax,minorax,posnang)
                else:
                    type='orig'
                    makedot(coords,freq,frequnit,file,imfitter,tele,type)
                    makegauss(coords,freq,frequnit,file,imfitter,tele,majorax,minorax,posnang)

                #fit a function to the dot to simply get a contour, but color it red
                for oi in oallsources:
                    if source in oi:
                        imfitloc2=f'{tele}/{oi}/fitfiles/{filename}_g.fit'
                if tele=='ALMA':
                    with silence_stderr_fd():        
                        imfit(imagename='testg.fits',model=imfitloc2,region=aregiontot)
                if tele!='ALMA':
                    with silence_stderr_fd():  
                        imfit(imagename='testg.fits',model=imfitloc2,region=regiontot)

                if conv_test==True:
                    writeandstore=1
                    color='b'
                    low,date=printrueconv(imfitter,'NA',writeandstore,color,imfitloc,newfilename,printimages,beamsize,freq,oicoords,frequnit,telescopes, tele, file, source, getchoosercoords, fixtochosenfreq, imagefound, prelimchecker, fit_sizenoloc, doublegaussfit,writenoise, force_upper_limit=force_upper_limit)
                    if not low==True and getchoosercoords==1:
                        foundsource=0
                        if len(choosercoordsarray) != 0:
                          i=-1
                          for ichoose in choosercoordsarray:
                              i=i+1
                              if source==ichoose[0]:
                                  foundsource=1

                                  if beamsize<ichoose[1]:
                                      choosercoordsarray[i]=[source,beamsize,file,newfilename, tele]
                        if foundsource==0:
                            choosercoordsarray.append([source,beamsize,file,newfilename, tele])



                    if date=='ONBOUNDRY':
                        if fixtochosenfreq==1:
                            badfiles.append([5,file])
                            plt.close()
                            print(f'ON BOUNDRY: {file}')
                            continue
                    if worked[0]==1:
                        if low==True:
                            #The double convolution may not be needed and should not be used.  Try a single imfitter image
                            #consider a region "multiplier" times the size of the beam for a compact source
                            multipliertot=dubmultval
                            offset1tot=majval*multipliertot
                            offset2tot=minval*multipliertot
                            regiontot=f'rotbox[{coords},[{offset1tot}arcsec,{offset1tot}arcsec],{rotval}{rotunit}]'
                            aregiontot=f'rotbox[{acoords},[{offset1tot}arcsec,{offset1tot}arcsec],{rotval}{rotunit}]'

                            iregion=f'guessregions/{index}region.txt'
                            if os.path.exists(iregion):
                                os.remove(iregion)
                            Path(iregion).touch()

                            fixer=masterfixerlow
                            worked=write_region(iregion,[pixcoord1[0]],[pixcoord2[0]],[apixcoord1[0]],[apixcoord2[0]],tele,file,regiontot,aregiontot,majval,minval,rotval,rotunit,writerreglowdub,fixer,doublegaussfit,second_gauss_pix, debugger, fixtochosenfreq, imagefound)
                            if worked[0]=='bad':
                                badfiles.append([7,file])
                                plt.close()
                                continue
                            with silence_stderr_fd():  
                                try:
                                    imfitter=imfit(imagename=file,model=imfitloc,region=regiontot,estimates=iregion)
                                except Exception as e:
                                    breaker=1
                            if breaker==1:
                                badfiles.append([3,file])
                                plt.close()
                                print(f'imfit does not work for {file} 3')
                            if imfitter['converged'][0]==True:
                                conv_test=True
                            if conv_test==True:
                                writeandstore=1
                                color='r'
                                low2,date2=printrueconv(imfitter,'NA',writeandstore,color,imfitloc,newfilename,printimages,beamsize,freq,oicoords,frequnit,telescopes, tele, file, source, getchoosercoords, fixtochosenfreq,imagefound, prelimchecker, fit_sizenoloc, doublegaussfit,writenoise, force_upper_limit=force_upper_limit)
                                if date2=='ONBOUNDRY':
                                    print(f'{file} on boundry')
                                    continue
                                #low=True
                if low==True:
                    if fixtochosenfreq==True and tele!='SMA':
                        #scalesize,fov_diag,bdiag
                        scalesize=compute_fov_and_beam(file)
                        if scalesize[1]<10:
                            badfiles.append([9,])
                            plt.close()
                            continue
                        else:
                            noise=imstat(imagename=file,region=regiontot)['rms'][0]
                            upbound=noise*3
                            fluxval=upbound*1000
                            fluxval=format_number(fluxval)
                            sfluxval=f'{fluxval}*'
                            sbeamsize=format_number(beamsize)
                            color='red'
                            pixlen=100
                            date=makeonlydotcontour(file,imfitloc,sfluxval,sbeamsize, newfilename, printimages,color,pixlen,oicoords,freq,frequnit,imfitter,telescopes,tele, getchoosercoords, fixtochosenfreq,imagefound, prelimchecker, fit_sizenoloc, doublegaussfit)

                    else:                            
                        noise=imstat(imagename=file,region=regionsmall)['rms'][0]
                        upbound=noise*3
                        fluxval=upbound*1000
                        fluxval=format_number(fluxval)
                        sfluxval=f'{fluxval}*'
                        sbeamsize=format_number(beamsize)
                        color='red'
                        pixlen=100
                        date=makedotcontour(file,imfitloc,sfluxval,sbeamsize, newfilename, printimages,color,pixlen,oicoords,freq,frequnit,imfitter,telescopes,tele, getchoosercoords, fixtochosenfreq,imagefound, prelimchecker, fit_sizenoloc, doublegaussfit)

                    if date=='ONBOUNDRY':
                        print(f'{file} on boundry')
                        continue 
                    #date=makeonlydotcontour(file, sbeamsize, sfluxval, newfilename, printimages)
                    blank='NA'
                    stotarea=f'{beamsize:.5f}'
                    sfreq=f'{freq:.2f}'
                    print(f"{source.ljust(15)}{tele.ljust(15)}{sfreq.ljust(19)}{sfluxval.ljust(16)}{date.ljust(15)}{stotarea.ljust(19)}{blank.ljust(15)}{file}")
                    write_fluxerr_row(source, tele, freq, date, None, file)

    if writesummary==True:
        sys.stdout = original_stdoutsum
    fluxerrfile_handle = None


# In[ ]:


#check badfiles
for ibad in badfiles:
    if not 'bigsize' == ibad[0]:
        print(ibad)

#for fixtochenfreq, the badfiles are all error 4 which means they didnt converge
#none of them are very important files either. We can investigate them more,
#potentially even just set their values as lower limits
#for now we exclude them from the analysis.


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

#for prelimchecker
if prelimchecker==1:
    bad1='VLA/NGC4258/22.4I1.92_AH0847_2004MAY22_1_2.38M55.7S.imfits'

    totmasterbad=[]
    totmasterbad.append(bad1)

    #check band files
    badfileversion='prelim'
    masterbadfiletxt=f'bad/masterbadfiles/masterbadfiles_{badfileversion}.txt'

    with open(masterbadfiletxt, 'w') as sumfile:
        if write==True:
            original_stdout = sys.stdout
            sys.stdout = sumfile

        for ibad in totmasterbad:
            print(ibad)

        if write==True:
            sys.stdout = original_stdout

#for fit_sizenoloc
#all files from badfiletxt
#all files are from bad location... which consist of 2 files
#these files will be added to the masterbadfilesfolder
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
            float(x[2]),  # Primary sorting key remains unchanged.
            # Handle 'NA' in x[3] by using a tuple where 'NA' sorts as a small number (e.g., -1).
            (-1 if x[4] == 'NA' else float(x[4].split('/')[0]) +
            float(x[4].split('/')[1]) / 12 +
            float(x[4].split('/')[2]) / 365),
            # Handle 'NA' in x[4] similarly, making it sort first if present.
            (-1 if x[5] == 'NA' else float(x[5]))
        ))


        for i2source_order in isource_order:
            print(i2source_order)
        print('\n')
    if write==True:
        sys.stdout = original_stdout


# In[ ]:


#setup chosenfreq image
#choosercoordsarray.append([source,beamsize,file,newfilename, tele])
if getchoosercoords==1:
    for ichoose in choosercoordsarray:
        #plot_filename = f'images/chosenfreq/{ichoose[0]}/{ichoose[3]}.pdf'

        base = os.path.basename(ichoose[2])
        tele=ichoose[4]
        
        fit_filename = f'images/{flagnum}/chosenfreq/{ichoose[0]}/{base}.fit'
        files = os.listdir(f'images/{flagnum}/chosenfreq/{ichoose[0]}')

        plot_filename_out = f'images/{flagnum}/chosenfreq/{ichoose[0]}/{tele}.{ichoose[3]}.pdf'
        fit_filename_out = f'images/{flagnum}/chosenfreq/{ichoose[0]}/{tele}.{base}.fit'

        counter=0
        for file in files:
          if ichoose[3] in file:
              counter=counter+1
              plot_filename = f'images/{flagnum}/chosenfreq/{ichoose[0]}/{file}'
          if file in fit_filename:
              counter=counter+1
        if counter!=2:
           pass
        else:

            shutil.move(plot_filename,plot_filename_out)
            shutil.move(fit_filename,fit_filename_out)
        files = os.listdir(f'images/{flagnum}/chosenfreq/{ichoose[0]}')
        for file in files:
          if 'pdf' in file:
             pdfin = f'images/{flagnum}/chosenfreq/{ichoose[0]}/{file}'
             pdfout = f'images/{flagnum}/chosenfreq/{file}'
             shutil.move(pdfin,pdfout)
          elif not file in fit_filename_out:
             shutil.rmtree(f'images/{flagnum}/chosenfreq/{ichoose[0]}/{file}')
             pass
          else:
              #print(fit_filename_out)
              pass

chosenfreq_dir = f'images/{flagnum}/chosenfreq'
if os.path.isdir(chosenfreq_dir):
    folders = os.listdir(chosenfreq_dir)
else:
    folders = []

for ifolder in folders:
        if 'pdf' in ifolder:
               continue
        files = os.listdir(f'images/{flagnum}/chosenfreq/{ifolder}')
        for ifile in files:
             if 'pdf' in ifile:
                    shutil.rmtree(f'images/{flagnum}/chosenfreq/{ifolder}')


# In[ ]:


#run get_chosenfreq.ipynb
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

if getchoosercoords==1:
    with open("get_chosenfreq.ipynb") as f:
        nb = nbformat.read(f, as_version=4)

    ep = ExecutePreprocessor(timeout=600)
    ep.preprocess(nb, {"metadata": {"path": "./"}})


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
