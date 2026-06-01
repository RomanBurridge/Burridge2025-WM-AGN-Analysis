#!/usr/bin/env python
# coding: utf-8

# In[1]:

import argparse
import json
import sys
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--flagnum", type=int, required=True)
parser.add_argument("--sciniter", type=int, required=True)
parser.add_argument("--date", type=int, required=True)
args = parser.parse_args()

flagnum = args.flagnum
sciniter = args.sciniter
date = args.date


# sciniter=0

#This flagnum must match the flagnum created using .copy_folder.py
# flag1 refers not flagging anything
# flag2 apply line flags with both sidebands
# flag3: USB image [6-8,9-11]
# flag4: LSB image [0-2,3-5]
# flagnum=2


# In[2]:


fundir2 = str(Path(__file__).resolve().parent)
binsize=128


# In[3]:


#excecutions
namer=f'{date}_{binsize}'

imsizer=binsize

#set the edges to be flagged
flagger=True 
#stop at each calibration table to check the logger.  Ensure solutions are found
inputs=False
#The calibration method can be changed in .calibrate.py
calnum=1 

#This determines the entire region by which the imfitter considers while conducting the fit.  It multiplies the size 
#of the major axis beam and is located in .obs_script_to_fields.py
multiplier=1.5

#This allows the images to be viewed in midst of the pipeline and is located in .fits.py
imviewer=False

#This will allow one science target to be chosen
chooser=True
chooser=False
choose_science='NGC4258'

if date==220131:
  bpchooser=True
  choose_bp='3c279'
  fluxchooser=True
  fluxchoice='mwc349a'

elif date==230316:
  bpchooser=True
  choose_bp='3c279'
  fluxchooser=True
  fluxchoice='Uranus'


# In[4]:


#make functions to clear tables and whatnot
from casatasks import tclean, rmtables
import os, glob, shutil

if flagnum == 1:
    print('Skipping clean run for flag1 because clean selections only use detected non-raw flags.')
    sys.exit(0)

def _is_selected_clean_detection(science_field: str) -> bool:
    # Clean every science target for the requested flag/date unless a manual
    # chooser is restricting the run to a single target.
    return (not chooser) or science_field == choose_science

def _rm_tree(path: str):
    try:
        shutil.rmtree(path, ignore_errors=True)
    except Exception:
        pass

def _remove_mask_tables(imagename: str):
    """Remove any existing CASA mask tables for this image."""
    for p in list(glob.glob(f"{imagename}.mask")) + list(glob.glob(f"{imagename}.mask*")):
        try:
            rmtables(p)
        except Exception:
            pass
        if os.path.exists(p):
            _rm_tree(p)

def _nuke_image_family(base: str):
    """Delete all CASA sidecar tables for an image family."""
    try:
        rmtables(base)
        rmtables(base + '.*')
    except Exception:
        pass
    _rm_tree(f'{base}.image')
    _rm_tree(f'{base}.model')
    _rm_tree(f'{base}.psf')
    _rm_tree(f'{base}.pb')
    _rm_tree(f'{base}.residual')
    _rm_tree(f'{base}.sumwt')
    _rm_tree(f'{base}.weight')
    _rm_tree(f'{base}.mask')
    _rm_tree(f'{base}.workdirectory')
    _rm_tree(f'{base}.flux')
    _rm_tree(f'{base}.image.pbcor')

def tclean_safe(*, imagename: str, **kwargs):
    """
    Safe tclean wrapper:
      • If mask=... is provided, remove any old {imagename}.mask first.
      • If no mask is provided but {imagename}.mask exists, reuse it.
      • On common RuntimeErrors (mask conflict, restfreq mismatch, etc.)
        delete stale products and retry once.
    """
    user_mask = kwargs.get('mask', None)
    user_usemask = kwargs.get('usemask', None)

    # Remove stale mask if user passes a new region
    if user_mask not in (None, ''):
        _remove_mask_tables(imagename)
        if user_usemask is None:
            kwargs['usemask'] = 'user'

    # Reuse existing mask if no new one supplied
    if (user_mask in (None, '')) and os.path.exists(f'{imagename}.mask'):
        kwargs['usemask'] = kwargs.get('usemask', 'user')
        kwargs['mask'] = ''

    try:
        return tclean(imagename=imagename, **kwargs)

    except RuntimeError as e:
        msg = str(e).lower()

        # Handle explicit mask conflict
        if 'mask image' in msg and 'also been supplied' in msg:
            _remove_mask_tables(imagename)
            if kwargs.get('mask') in (None, ''):
                kwargs['usemask'] = 'user'
                kwargs['mask'] = ''
            return tclean(imagename=imagename, **kwargs)

        # Handle coordinate/restfreq/stale image conflicts
        triggers = (
            'coordinate system mismatch',
            'active rest frequencies',
            'cannot open existing image',
            'image exists and is not compatible',
        )
        if any(t in msg for t in triggers):
            _nuke_image_family(imagename)
            if kwargs.get('mask', None) not in (None, ''):
                _remove_mask_tables(imagename)
                kwargs['usemask'] = 'user'
            return tclean(imagename=imagename, **kwargs)

        raise


# In[5]:


#obs_script_to_fields
import os
import sys
from collections import defaultdict
from casatasks import clearcal
from casatasks import listobs
from casatasks import tclean
from casatasks import imhead
from casatasks import imstat
from casatasks import rmtables

os.chdir(fundir2)

clear=True
write=False
printer=True
tcleaner=True

input_folder=f'{fundir2}/sma_data/{namer}_master.ms/'
new_folder_name = input_folder.replace("master", "")
new_folder_name = new_folder_name.replace("sma_data/", "")
parent_folder_name = new_folder_name.replace("_.ms", "")


prefix='uncal_data'
if not os.path.exists(parent_folder_name):
    os.mkdir(str(parent_folder_name))
if not os.path.exists(f'{parent_folder_name}/uncal_data'):
    os.mkdir(f'{parent_folder_name}/{prefix}')

outfile=f'{parent_folder_name}{date}_summary.txt'
source_text =f'{fundir2}/txt_{date}.txt'

with open(outfile, 'w') as file:

    if write==True:
        original_stdout = sys.stdout
        sys.stdout = file
    myvis=input_folder
    if not os.path.exists(myvis):
        input(f'{myvis} does not exist.  please make it')
    if clear==True:
        clearcal(myvis)
    listobs_output = listobs(myvis)

    fields=[]
    for i in range(listobs_output['nfields']):
        fields.append(listobs_output[f'field_{i}']['name'])

    #Finding largest and smallest scanum
    # Initialize variables to store the highest and lowest scanum
    highest_scanum = float('-inf')  # Initialize to negative infinity, suitable for finding maximum
    lowest_scanum = float('inf')    # Initialize to positive infinity, suitable for finding minimum

    # Iterate over scans in listobs_output
    for scan in listobs_output:
        # Split the current scan using 'scan_' as the separator
        split_item = scan.split('scan_')

        # Check if the split resulted in two parts
        if len(split_item) == 2:
            # Extract the scanum from the second part after 'scan_'
            scanum = int(split_item[1])

            # Update the highest and lowest scanum accordingly
            highest_scanum = int(max(highest_scanum, scanum))
            lowest_scanum = int(min(lowest_scanum, scanum))

    # Print the highest and lowest scanum
    #print("Highest scanum:", highest_scanum)
    #print("Lowest scanum:", lowest_scanum)

    scans=[]
    scanames=[]
    for i in range(lowest_scanum,highest_scanum):
        #print(listobs_output[f"scan_{i}"].keys())
        if set(listobs_output[f"scan_{i}"].keys()) != {'0'}:
            if printer==True:   print(listobs_output[f"scan_{i}"].keys())
            input(f'scan number {i} has a dictionary key other than 0...')

        # Get the FieldName for the current scan
        iscan=listobs_output[f"scan_{i}"]['0']
        i2scan=[]
        field_name = iscan['FieldName']
        scanames.append(field_name)
        i2scan.append(field_name)
        i2scan.append(f"{iscan['IntegrationTime']:.1f}")

        i2scan.append(iscan['scanId'])

        # Check if the FieldName exists in the 'fields' list
        if field_name in fields:
            scans.append(i2scan)
        else:
            if printer==True:   print(f"Field name '{field_name}' not found in the 'fields' list.")

    #Now we will extract the science fields from the observing script

    tfields=[]
    cfields=[]
    bpfields=[]
    fluxfields=[]
    tcoords=[]
    pcoords=[]
    try:
        with open(source_text, 'r') as infile:
            for line in infile:
                line = line.strip()
                if '=' not in line:
                    continue  # Skip lines without an equal sign

                key, value = line.split('=', 1)
                value = value.split('"')[1]

                if 'targ' in key:
                    science=[key.split('targ')[1].split('="')[0], value.split(' -r')[0]]
                    breaker=0
                    if chooser==True and science[1]!=choose_science:
                        breaker=1
                    if breaker==1:
                        continue
                    icoords=[]
                    icoords.append(((value.split('-r'))[1].split('-d')[0]).replace(" ",""))
                    icoords.append((((value.split('-r'))[1].split('-d')[1]).split(" -e")[0]).replace(" ",""))
                    tcoords.append(icoords)
                    tfields.append(science)
                elif 'cal' in key:
                    cfields.append([key.replace('$cal',''), value])
                elif 'bpass' in key and bpchooser==False:
                    bp=[key.replace('$bpass',''), value]
                    breaker=0
                    bpfields.append(bp)
                elif 'flux' in key and fluxchooser==False:
                    flux=[key.replace('$flux',''), value]
                    breaker=0
                    fluxfields.append(flux)

    except FileNotFoundError:
        if printer==True:   print("File not found or path is incorrect.")
    except IOError:
        if printer==True:   print("An error occurred while reading the file.")
    if printer==True:   
        if len(tcoords)==len(tfields):  print('target coordinates same length as coordinates')
        else:   print('error')


    if bpchooser==True:
        bpfields.append(['0',choose_bp])
    if fluxchooser==True:
        fluxfields.append(['0',fluxchoice])

    # Assuming tfields and cfields are defined
    # Create dictionaries to store the data in a more organized manner
    tfields_dict = {field[0]: field[1] for field in tfields}
    cfields_dict = {}
    for icfield in cfields:
        cdub = icfield[0].split("_")
        for icdub in cdub:
            cfields_dict[icdub] = icfield[1]

    # Matching up Fields and Targets
    tc_match = []
    if printer==True:   print("Unfiltered Fields\n  - target         calibrator")
    #tfield must be nested outside of cfield for this search so that coordinate order is preserved
    for tfield_key, tfield_value in tfields_dict.items():
        if tfield_key in cfields_dict:
            itc_match = [tfield_value, cfields_dict[tfield_key]]
            tc_match.append(itc_match)
    for itc_match in tc_match: 
        if printer==True:   print(f"  - {itc_match[0]}        {itc_match[1]}")

    # Now to look through scans
    missing_in_scans = [field for field in fields if field not in scanames]

    if missing_in_scans:
        print("\nFields found in 'fields' array but not in 'scans' array:")
        for missing_field in missing_in_scans:
            print(f"  - {missing_field}")
    else:
        if printer==True:   print("\nAll fields in the 'fields' array are found in the 'scans' array.")

    scanfields = []
    scanids_by_field = defaultdict(list)

    for iscan in scans:
        ifield = iscan[0]
        scanids_by_field[ifield].append(iscan[2])

    for ifield in fields:
        num_scans = 0
        integration_value = 'random number'
        scanids = scanids_by_field[ifield]

        for iscan in scans:
            if iscan[0] == ifield:
                num_scans += 1
                if integration_value != iscan[1] and integration_value != 'random number':
                    input(f"{ifield} has multiple integration values")
                integration_value = iscan[1]

        iscanfield = [ifield, num_scans, integration_value]
        if scanids:
            iscanfield.append(scanids)
        scanfields.append(iscanfield)

    # Convert the observed field lists to sets for faster lookups
    observed_fields = set()
    for observed_list in (tfields, cfields, bpfields, fluxfields):
        observed_fields.update(field[1] for field in observed_list)

    # Find fields that were observed but not in the observing script
    not_observed_fields = set(fields) - observed_fields

    if not_observed_fields:
        if printer==True:   print(f"\nFields that were observed but not in the observing script:\n  Field     TotNum    Length      Num")
        for field in not_observed_fields:
            for iscanfield in scanfields:
                if iscanfield[0]==field:
                    count=iscanfield[1]
                    scans_time=iscanfield[2]
                    if scans_time!='random number':
                        scans=iscanfield[3]
                    else:
                        scans='NA'
            if printer==True:   print(f"  - {field}     {count}     {scans_time}     {scans}")
    else:
        if printer==True:   print("\nAll fields in the observing script were found in observed data.")

    # Extract the second element (field name) from each tuple in the lists
    combined_fields = [field[1] for field in tfields + cfields + bpfields + fluxfields]

    # Check if there are any fields in the combined list that are not in the 'fields' array
    missing_in_fields = [field for field in combined_fields if field not in fields]

    if missing_in_fields:
        if printer==True:   print("\nFields found in observing script but not in observed data:")
        for missing_field in missing_in_fields:
            if printer==True:   print(f"  - {missing_field}")
    else:
        if printer==True:   print("\nAll fields found in observing script were found in observed data")

    targets_not_observed = [field for field in combined_fields if field not in fields and field in [match[0] for match in tc_match]]
    calibrators_not_observed = [field for field in combined_fields if field not in fields and field not in targets_not_observed]

    if printer==True:   print("\nTargets found in observing script but not observed:")
    for missing_target in targets_not_observed:
        if printer==True:   print(f"  - {missing_target}")

    if printer==True:   print("\nCalibrators found in observing script but not observed:")
    for missing_calibrator in calibrators_not_observed:
        if printer==True:   print(f"  - {missing_calibrator}")

    if printer==True:   print("\nCalibrators found in observing script but not observed, with corresponding targets:")
    for missing_calibrator in calibrators_not_observed:
        corresponding_target = next((match[0] for match in tc_match if match[1] == missing_calibrator), "None")
        count=0
        for iscanfield in scanfields:
            if iscanfield[0]==corresponding_target:
                count=iscanfield[1]
                scans_time=iscanfield[2]
                scans=iscanfield[3]
        if count==0 and printer==True:   print(f"  - Calibrator: {missing_calibrator} | Target: {corresponding_target} was also unobserved")
        elif printer==True:print(f"  - Calibrator: {missing_calibrator} | Target: {corresponding_target}\n       {corresponding_target}: scanum:{count}, scanlength:{scans_time}, scan:{scans}")
    # Filter tc_match to include only entries where both target and calibrator were observed
    f_tc_match = [match for match in tc_match if match[0] in fields and match[1] in fields]
    f_bp_match = [match for match in bpfields if match[1] in fields]
    f_flux_match = [match for match in fluxfields if match[1] in fields]

    unobserved_targets = [pair[0] for pair in tc_match if pair not in f_tc_match]
    unobserved_targets = [pair for pair in tc_match if pair not in f_tc_match]
    unobserved_bp = [singf for singf in bpfields if singf not in f_bp_match]
    unobserved_flux = [singf for singf in fluxfields if singf not in f_flux_match]

    # Printing the unobserved targets
    if printer==True:   print("\nUnobserved targets:")
    for target in unobserved_targets:
        if printer==True:   print(f"  - {target[0]}")

    if len(unobserved_bp)!=0:
        for iscanfield in scanfields:
            if printer==True:   print("\nUnobserved bandpass:")
            for bp in unobserved_bp:
                if printer==True:   print(f"  - {bp}")

    if len(unobserved_flux)!=0:
        if printer==True:   print("\nUnobserved flux:")
        for flux in unobserved_flux:
            if printer==True:   print(f"  - {flux}")

    if printer==True:   print("\nFiltered Fields:")
    if printer==True:   print(f'{len(f_tc_match)} out of {len(tc_match)} observed')
    if len(f_tc_match)!=0:
        if printer==True:   print("    target         calibrator")
    for if_tc_match in f_tc_match: 
        for iscanfield in scanfields:
            if iscanfield[0]==if_tc_match[0]:
                s_count=iscanfield[1]
                s_scans_time=iscanfield[2]
                s_scans=iscanfield[3]
            if iscanfield[0]==if_tc_match[1]:
                c_count=iscanfield[1]
                c_scans_time=iscanfield[2]
                c_scans=iscanfield[3]
        if printer==True:   
            print(f"  - {if_tc_match[0]}: scanum:{s_count}, scanlength:{s_scans_time}, scan:{s_scans}")
            print(f"  - {if_tc_match[1]}: scanum:{c_count}, scanlength:{c_scans_time}, scan:{c_scans}\n")

    if printer==True:   print(f'\nBandpass: {len(f_bp_match)} out of {len(bpfields)} observed')
    if printer==True:   print("  Field     TotNum    Length      Num")
    for ibpfield in f_bp_match:
        for iscanfield in scanfields:
            if iscanfield[0]==ibpfield[1]:
                count=iscanfield[1]
                scans_time=iscanfield[2]
                if scans_time!='random number':
                    scans=iscanfield[3]
        if printer==True:   print(f"  - {ibpfield[1]}    {count}        {scans_time}        {scans}")

    if printer==True:   print(f'\nFlux: \n{len(f_flux_match)} out of {len(fluxfields)} observed')
    if printer==True:   print("  Field     TotNum    Length      Num")
    for ifluxfield in f_flux_match:
        for iscanfield in scanfields:
            if iscanfield[0]==ifluxfield[1]:
                count=iscanfield[1]
                scans_time=iscanfield[2]
                if len(iscanfield)>3:
                    scans=iscanfield[3]
                else:
                    scans='NA'
        if printer==True:   print(f"  - {ifluxfield[1]}    {count}        {scans_time}        {scans}")
    if write==True:
        sys.stdout = original_stdout
fc_prefix=f'{namer}/images/clean/flag{flagnum}_cal{calnum}/'
getfit_params=[]
bpgetfit_params=[]
fluxgetfit_params=[]

for itc_match in tc_match:
    science_field = itc_match[0]
    phase_field=itc_match[1]

    temploc=f'sci_temp_{science_field}'
    ptemploc=f'p_temp_{phase_field}'

    if os.path.exists(f'{temploc}.psf'):
        os.system(f'rm -rf {temploc}*')

    if os.path.exists(f'{ptemploc}.psf'):
        os.system(f'rm -rf {ptemploc}*')

j=0
for itc_match in tc_match:
    tcoord=tcoords[j]
    j=j+1

    science_field = itc_match[0]
    phase_field=itc_match[1]

    if not _is_selected_clean_detection(science_field):
        continue

    sciloc=f'{fc_prefix}fields/science/{science_field}_test.image'
    phaseloc=f'{fc_prefix}fields/phasecal/{phase_field}_test.image'

    breaker2=0
    for iut in unobserved_targets:
        if science_field==iut[0]:
            breaker2=1
    if breaker2==1:
        print(f'no observed target for {science_field}')
        continue

    if tcleaner==True:
        temploc=f'sci_temp_{science_field}'
        tclean_safe(vis=myvis, 
            field=science_field,
            imagename=temploc,
            imsize=imsizer, 
            cell='1arcsec', 
            specmode='mfs',
            niter=0,
            usemask='pb')

        ptemploc=f'p_temp_{phase_field}'
        tclean_safe(vis=myvis, 
            field=phase_field,
            imagename=ptemploc,
            imsize=imsizer, 
            cell='1arcsec', 
            specmode='mfs',
            niter=0,
            usemask='pb')

        getfit_param=[tcoord,f'{temploc}.image',sciloc,f'{ptemploc}.image',phaseloc]
        getfit_params.append(getfit_param)
####################

for ibp in bpfields:
    bp_field = ibp[1]
    bploc=f'{fc_prefix}fields/bpcal/{bp_field}_test.image'
    if tcleaner==True:
        bptemploc=f'bp_temp_{bp_field}'
        tclean_safe(vis=myvis, 
            field=bp_field,
            imagename=bptemploc,
            imsize=imsizer, 
            cell='1arcsec', 
            specmode='mfs',
            niter=0,
            usemask='pb')

        bpgetfit_param=[f'{bptemploc}.image',bploc]
        bpgetfit_params.append(bpgetfit_param)
####################
#flux
for iflux in fluxfields:
    flux_field = iflux[1]
    fluxloc=f'{fc_prefix}fields/fluxcal/{flux_field}_test.image'
    if tcleaner==True:
        fluxtemploc=f'flux_temp_{flux_field}'
        tclean_safe(vis=myvis, 
            field=flux_field,
            imagename=fluxtemploc,
            imsize=imsizer, 
            cell='1arcsec', 
            specmode='mfs',
            niter=0,
            usemask='pb')

        fluxgetfit_param=[f'{fluxtemploc}.image',fluxloc]
        fluxgetfit_params.append(fluxgetfit_param)

####################
if tcleaner==True:
    regions=[]
    pregions=[]
    bpregions=[]
    fluxregions=[]

    j=0
    for iparam in getfit_params:
        isum_vals=[]
        coord=iparam[0]
        imfile=iparam[1]
        majunit=imhead(imagename=imfile)['restoringbeam']['major']['unit']
        majval=imhead(imagename=imfile)['restoringbeam']['major']['value']

        minunit=imhead(imagename=imfile)['restoringbeam']['minor']['unit']
        minval=imhead(imagename=imfile)['restoringbeam']['minor']['value']

        rotunit=imhead(imagename=imfile)['restoringbeam']['positionangle']['unit']
        rotval=imhead(imagename=imfile)['restoringbeam']['positionangle']['value']

        if majunit!='arcsec' or minunit!='arcsec':
            input(f'Beam Axis does not use arcsec for units.  check imhead for:\n{imfile}')

        if rotunit!='deg':
            input(f'Beam Rotation does not use deg for units.  check imhead for:\n{imfile}')

        icoord1=coord[0].split(':')
        coord1=f'{icoord1[0]}h{icoord1[1]}m{icoord1[2]}s'
        coord2=coord[1].replace(':','.')
        #majval_m=majval*2
        majval_m=majval*multiplier
        minval_m=minval*multiplier

        region=f'rotbox[[{coord1},{coord2}],[{majval_m}{majunit},{minval_m}{minunit}],{rotval}{rotunit}]'
        regions.append(region)

        ##################
        #phase calibrator center region
        pimfile=iparam[3]
        blc=imstat(imagename=pimfile)['blc']
        trc=imstat(imagename=pimfile)['trc']
        pcoord1=(blc[0]+trc[0])/2
        pcoord2=(blc[1]+trc[1])/2

        #phase calibrator region format
        pisum_vals=[]
        pmajunit=imhead(imagename=pimfile)['restoringbeam']['major']['unit']
        pmajval=imhead(imagename=pimfile)['restoringbeam']['major']['value']

        pminunit=imhead(imagename=pimfile)['restoringbeam']['minor']['unit']
        pminval=imhead(imagename=pimfile)['restoringbeam']['minor']['value']

        protunit=imhead(imagename=pimfile)['restoringbeam']['positionangle']['unit']
        protval=imhead(imagename=pimfile)['restoringbeam']['positionangle']['value']

        if pmajunit!='arcsec' or pminunit!='arcsec':
            input(f'Beam Axis does not use arcsec for units.  check imhead for:\n{pimfile}')

        if protunit!='deg':
            input(f'Beam Rotation does not use deg for units.  check imhead for:\n{pimfile}')

        pmajval_m=pmajval*2
        pminval_m=pminval*2

        pregion=f'rotbox[[{pcoord1}pix,{pcoord2}pix],[{pmajval_m}{pmajunit},{pminval_m}{pminunit}],{protval}{protunit}]'
        pregions.append(pregion)
        #########################
        ##################
    #bandpass calibrator center region
    for iparam in bpgetfit_params:
        imfile=iparam[0]
        blc=imstat(imagename=imfile)['blc']
        trc=imstat(imagename=imfile)['trc']
        coord1=(blc[0]+trc[0])/2
        coord2=(blc[1]+trc[1])/2

        #calibrator region format
        isum_vals=[]
        majunit=imhead(imagename=imfile)['restoringbeam']['major']['unit']
        majval=imhead(imagename=imfile)['restoringbeam']['major']['value']

        minunit=imhead(imagename=imfile)['restoringbeam']['minor']['unit']
        minval=imhead(imagename=imfile)['restoringbeam']['minor']['value']

        rotunit=imhead(imagename=imfile)['restoringbeam']['positionangle']['unit']
        rotval=imhead(imagename=imfile)['restoringbeam']['positionangle']['value']

        if majunit!='arcsec' or minunit!='arcsec':
            input(f'Beam Axis does not use arcsec for units.  check imhead for:\n{imfile}')

        if rotunit!='deg':
            input(f'Beam Rotation does not use deg for units.  check imhead for:\n{imfile}')

        majval_m=majval*2
        minval_m=minval*2

        bpregion=f'rotbox[[{coord1}pix,{coord2}pix],[{majval_m}{majunit},{minval_m}{minunit}],{rotval}{rotunit}]'
        bpregions.append(bpregion)
        #########################
        #bandpass calibrator center region
    for iparam in fluxgetfit_params:
        imfile=iparam[0]
        blc=imstat(imagename=imfile)['blc']
        trc=imstat(imagename=imfile)['trc']
        coord1=(blc[0]+trc[0])/2
        coord2=(blc[1]+trc[1])/2
        #calibrator region format
        isum_vals=[]
        majunit=imhead(imagename=imfile)['restoringbeam']['major']['unit']
        majval=imhead(imagename=imfile)['restoringbeam']['major']['value']

        minunit=imhead(imagename=imfile)['restoringbeam']['minor']['unit']
        minval=imhead(imagename=imfile)['restoringbeam']['minor']['value']

        rotunit=imhead(imagename=imfile)['restoringbeam']['positionangle']['unit']
        rotval=imhead(imagename=imfile)['restoringbeam']['positionangle']['value']

        if majunit!='arcsec' or minunit!='arcsec':
            input(f'Beam Axis does not use arcsec for units.  check imhead for:\n{imfile}')

        if rotunit!='deg':
            input(f'Beam Rotation does not use deg for units.  check imhead for:\n{imfile}')

        majval_m=majval*2
        minval_m=minval*2

        fluxregion=f'rotbox[[{coord1}pix,{coord2}pix],[{majval_m}{majunit},{minval_m}{minunit}],{rotval}{rotunit}]'
        fluxregions.append(fluxregion)
        #########################

    for itc_match in tc_match:
        science_field = itc_match[0]
        phase_field=itc_match[1]
        temploc=f'sci_temp_{science_field}'
        ptemploc=f'p_temp_{phase_field}'
        if os.path.exists(f'{temploc}.psf'):
            os.system(f'rm -rf {temploc}*')
        if os.path.exists(f'{ptemploc}.psf'):
            os.system(f'rm -rf {ptemploc}*')

    for ibp in bpfields:
        bp_field = ibp[0]
        bptemploc=f'bp_temp_{ibp[1]}'
        if os.path.exists(f'{bptemploc}.psf'):
            os.system(f'rm -rf {bptemploc}*')   

    for iflux in fluxfields:
        flux_field = iflux[0]
        fluxtemploc=f'flux_temp_{iflux[1]}'
        if os.path.exists(f'{fluxtemploc}.psf'):
            os.system(f'rm -rf {fluxtemploc}*')   


# In[6]:


#print SPW and channels
from casatools import table
import numpy as np

myvis = f'{fundir2}/{namer}/uncal_data/{namer}_{flagnum}.ms'

#store back up 
src_ms = os.path.join(fundir2, "sma_data", f"{namer}_master.ms")
dst_ms = os.path.join(f'{fundir2}/{namer}/uncal_data', f"{namer}_1.ms")
if not os.path.exists(dst_ms):
    shutil.copytree(src_ms, dst_ms)

if flagnum!=1:
        if os.path.exists(myvis):
          shutil.rmtree(myvis)
        shutil.copytree(f'{fundir2}/{namer}/uncal_data/{namer}_1.ms',myvis)



ms = myvis
tb = table()
tb.open(f'{ms}/SPECTRAL_WINDOW')

all_freqs = []   # we’ll extend this list with edges instead of just centers

print('=== SPWs 0–5 ===')
for spw in range(0, 6):
        freqs = tb.getcell('CHAN_FREQ', spw).astype(float)
        widths = tb.getcell('CHAN_WIDTH', spw)
        # ensure widths are arrays of |w|
        if hasattr(widths, '__len__'):
                w = np.abs(np.asarray(widths, dtype=float))
        else:
                w = np.full_like(freqs, abs(float(widths)), dtype=float)
        # compute edges per channel and extend the list
        all_freqs.extend((freqs - 0.5 * w).tolist())
        all_freqs.extend((freqs + 0.5 * w).tolist())

        spw_low = np.min(freqs - 0.5 * w)
        spw_high = np.max(freqs + 0.5 * w)
        bw = (spw_high - spw_low) / 1e9
        print(f'SPW {spw}: {spw_low/1e9:.3f}–{spw_high/1e9:.3f} GHz  (width = {bw:.3f} GHz)')

print('\n=== SPWs 6–12 ===')
for spw in range(6, min(13, tb.nrows())):
        freqs = tb.getcell('CHAN_FREQ', spw).astype(float)
        widths = tb.getcell('CHAN_WIDTH', spw)
        if hasattr(widths, '__len__'):
                w = np.abs(np.asarray(widths, dtype=float))
        else:
                w = np.full_like(freqs, abs(float(widths)), dtype=float)
        all_freqs.extend((freqs - 0.5 * w).tolist())
        all_freqs.extend((freqs + 0.5 * w).tolist())

        spw_low = np.min(freqs - 0.5 * w)
        spw_high = np.max(freqs + 0.5 * w)
        bw = (spw_high - spw_low) / 1e9
        print(f'SPW {spw}: {spw_low/1e9:.3f}–{spw_high/1e9:.3f} GHz  (width = {bw:.3f} GHz)')

# Now overall range uses the extended list of edges
lowest = min(all_freqs) / 1e9
highest = max(all_freqs) / 1e9
print(f'\nOverall frequency range: {lowest:.3f}–{highest:.3f} GHz')

tb.close()



# In[7]:


#molecular transitions
with open(f"{fundir2}/molecular_transitions.txt", "r") as file:
    for line in file:
        if 'http' in line:
            print(line)
        words = line.strip().split()  # split by spaces
        segments = [words[i:i+3] for i in range(0, len(words), 3)]

thresh=1000
threshvals=[]
for i  in segments:
    if float(i[2])>thresh:
        threshvals.append(i)

# sort by frequency (i[1])
threshvals.sort(key=lambda x: float(x[1]))

# print in ascending order
for i in threshvals:
    print(i)


# In[8]:


#Get flagging region from distances for transitions
#distances for Pesce, 2018 ordered to Burridge, 2025
import numpy as np

totshifts=[]

val=466.87
valerr=0.49
dist1=val+valerr
dist2=val-valerr
totshifts.append(["4258",f"{dist1:.2f}",f"{dist2:.2f}"])

val=4088.8
valerr=5.3
dist1=val+valerr
dist2=val-valerr
totshifts.append(["1194",f"{dist1:.2f}",f"{dist2:.2f}"])

val=8322.22
valerr=1.13
dist1=val+valerr
dist2=val-valerr
totshifts.append(["5765b",f"{dist1:.2f}",f"{dist2:.2f}"])

val=1850.8
valerr1=13.5
valerr2=13.9
dist1=val+valerr2
dist2=val-valerr2
totshifts.append(["2273",f"{dist1:.2f}",f"{dist2:.2f}"])

val=10189.26
valerr=1.20
dist1=val+valerr
dist2=val-valerr
totshifts.append(["6264",f"{dist1:.2f}",f"{dist2:.2f}"])

val=3259.75
valerr=1.00
dist1=val+valerr
dist2=val-valerr
totshifts.append(["3789",f"{dist1:.2f}",f"{dist2:.2f}"])

val=4954.5
valerr=15
dist1=val+valerr
dist2=val-valerr
totshifts.append(["2960",f"{dist1:.2f}",f"{dist2:.2f}"])

val=7618.2
valerr=14.0
dist1=val+valerr
dist2=val-valerr
totshifts.append(["558",f"{dist1:.2f}",f"{dist2:.2f}"])

val=7834.28
valerr1=2.1
valerr2=2.2
dist1=val+valerr2
dist2=val-valerr2
totshifts.append(["6323",f"{dist1:.2f}",f"{dist2:.2f}"])

val=4818
valerr=10.5
dist1=val+valerr
dist2=val-valerr
totshifts.append(["437",f"{dist1:.2f}",f"{dist2:.2f}"])

#distance for other sources 

#NGC 1068 
#Baer-Way, 2024
#2024ApJ...964..172B
val=1137
valerr=3
dist1=val+valerr
dist2=val-valerr
totshifts.append(["1068",f"{dist1:.2f}",f"{dist2:.2f}"])

#NGC 3393 
#Koss, 2022
#2022ApJS..261....6K
val=3833
valerr=5
dist1=val+valerr
dist2=val-valerr
totshifts.append(["3393",f"{dist1:.2f}",f"{dist2:.2f}"])

#Circinus
#	Koribalski, 2004
#2004AJ....128...16K
val=434
valerr=3
dist1=val+valerr
dist2=val-valerr
totshifts.append(["Circinus",f"{dist1:.2f}",f"{dist2:.2f}"])

#CGCG 074-06
#Albareti, 2017
#2016SDSSD.C...0000:
val=6915
valerr=3
dist1=val+valerr
dist2=val-valerr
totshifts.append(["074",f"{dist1:.2f}",f"{dist2:.2f}"])

#NGC 4388
#Lu, 1993
#1993ApJS...88..383L
val=2524
valerr=1
dist1=val+valerr
dist2=val-valerr
totshifts.append(["4388",f"{dist1:.2f}",f"{dist2:.2f}"])

#NGC 4945
#Kanehisa, 2023
#2023MNRAS.519.6184K
val=563
valerr=3
dist1=val+valerr
dist2=val-valerr
totshifts.append(["4945",f"{dist1:.2f}",f"{dist2:.2f}"])

#UGC 6093
#Albareti, 2017
#2016SDSSD.C...0000:
val=10803
valerr=3
dist1=val+valerr
dist2=val-valerr
totshifts.append(["6093",f"{dist1:.2f}",f"{dist2:.2f}"])

#NGC 3079
#Springbob, 2005
#2005ApJS..160..149S
val=1106
valerr=1
dist1=val+valerr
dist2=val-valerr
totshifts.append(["3079",f"{dist1:.2f}",f"{dist2:.2f}"])

#NGC 1320
#Theureau, 2017
#2017A&A...599A.104T
val=2783
valerr=10
dist1=val+valerr
dist2=val-valerr
totshifts.append(["1320",f"{dist1:.2f}",f"{dist2:.2f}"])

#IC 2560
#Pesce, 2015
#2015ApJ...810...65P
val=2925
valerr=2
dist1=val+valerr
dist2=val-valerr
totshifts.append(["2560",f"{dist1:.2f}",f"{dist2:.2f}"])

#Mrk 1029
#Huchra, 1999
#1999ApJS..121..287H
val=9076
valerr=32
dist1=val+valerr
dist2=val-valerr
totshifts.append(["1029",f"{dist1:.2f}",f"{dist2:.2f}"])

from astropy.constants import c
import astropy.units as u

c_kms = c.to(u.km / u.s).value    # convert c to km/s

totflags=[]


for i in totshifts:
    #beta1 is smaller than beta2
    val1=float(i[1])
    beta1 = val1 / c_kms
    #z1 is smaller than z2
    z1 = ((1 + beta1) / (1 - beta1))**0.5 - 1
    val2=float(i[2])
    beta2 = val2 / c_kms
    z2 = ((1 + beta2) / (1 - beta2))**0.5 - 1
    itotflags=[]
    source=i[0]
    itotflags.append(source)
    iitotflags=[]
    for j in threshvals:
      mol=j[0]
      molval=float(j[1])
      #z1 is smaller than z2
      molval1 = molval/(1+z1)
      molval2 = molval/(1+z2)

      #width1 is smaller than width2
      #just used width2
      width1=150/c_kms*molval1
      width2=150/c_kms*molval2
      #smaller
      newmolval1=molval1-width2*2
      #bigger
      newmolval2=molval2+width2*2
      #if the smallest is bigger than highest or if bigger is smaller than the lowest
      #if newmolval1=>highest or newmolval2=<lowest: continue
      if newmolval1>=highest or newmolval2<=lowest: 
         continue
      else:
            iitotflags.append([mol,newmolval1,newmolval2])
    itotflags.append(iitotflags)
    totflags.append(itotflags)

for i in totflags:
    j=-1
    print()
    for ii in i:
        j=j+1
        if j!=0:
          for iii in ii:
              print(iii[0])
              diff=iii[2]-iii[1]
              print(f'{diff:.2f}')
          print()
        else:
            print(ii)


# In[9]:


#getting flag regions for each source
import numpy as np
from casatools import table

myvis = f'{namer}/uncal_data/{namer}_{flagnum}.ms'
ms = myvis
tb = table()
tb.open(f'{ms}/SPECTRAL_WINDOW')

bands = [(209.367e9, 221.667e9), (229.367e9, 241.667e9)]
# two observing sidebands (in Hz)
bands = [(209.367e9, 221.667e9), (229.367e9, 241.667e9)]

# compute total and per-band widths
for i, (lo, hi) in enumerate(bands, start=1):
    bw = hi - lo
    print(f"Sideband {i}: {lo/1e9:.3f}–{hi/1e9:.3f} GHz  (width = {bw/1e9:.3f} GHz)")

# total continuous bandwidth (sum of both sidebands)
total_bw = sum(hi - lo for lo, hi in bands)
print(f"\nTotal continuous bandwidth before flagging: {total_bw/1e9:.3f} GHz")


totflags2=[]
for i in totflags:
    print()
    src, lines = i[0], i[1]
    print(src)
    allflagvals=[]
    for iii in lines:
        lowerval = iii[1] * 1e9   # GHz -> Hz
        upperval = iii[2] * 1e9   # GHz -> Hz
        for spw in range(tb.nrows()):
            freqs = tb.getcell('CHAN_FREQ', spw)         # Hz (centers)
            widths = tb.getcell('CHAN_WIDTH', spw)       # Hz (scalar or array)

            # make per-channel |width| array
            if hasattr(widths, '__len__'):
                w_arr = np.abs(np.asarray(widths, dtype=float))
            else:
                w_arr = np.full_like(freqs, abs(float(widths)), dtype=float)

            flagged = []
            for ch, fc in enumerate(freqs):
                w = w_arr[ch]
                ch_lo = fc - 0.5 * w
                ch_hi = fc + 0.5 * w
                if max(ch_lo, lowerval) <= min(ch_hi, upperval):  # overlap
                    flagged.append((ch, fc, w, ch_lo, ch_hi))

            if flagged:
                lowest_ch = min(ch for ch, _, _, _, _ in flagged)
                highest_ch = max(ch for ch, _, _, _, _ in flagged)

                # find corresponding frequencies (edges)
                lowest_freq = min(ch_lo for _, _, _, ch_lo, _ in flagged)
                highest_freq = max(ch_hi for _, _, _, _, ch_hi in flagged)

                print(f'{spw}:{lowest_ch}~{highest_ch}')
                totflags2.append([src,int(spw),int(lowest_ch),int(highest_ch)])                       # channel range
                flagvals=[lowest_freq,highest_freq]
                allflagvals.append(flagvals)

                import numpy as np

    # ---------------------------------------------------------------------
    # DEFINE YOUR SIDE BAND COVERAGE (in Hz)
    # ---------------------------------------------------------------------
    # Each tuple represents (lowest frequency, highest frequency) for a sideband
    # These are the total usable frequency ranges of your data.
    bands = [
        (209.367e9, 221.667e9),   # Lower sideband
        (229.367e9, 241.667e9)    # Upper sideband
    ]

    # ---------------------------------------------------------------------
    # FUNCTION 1: Merge overlapping or touching frequency intervals
    # ---------------------------------------------------------------------
    def merge_intervals(intervals):
        """
        Combine any overlapping or adjacent [lo, hi] intervals into a single one.
        Example:
            [(1, 3), (2, 4), (5, 6)] -> [(1, 4), (5, 6)]
        This prevents double-counting when intervals overlap.
        """
        if not intervals:
            return []

        # Clean up: make sure each interval is ordered correctly (lo < hi)
        # and sort them by their starting frequency.
        ints = [(min(a, b), max(a, b)) for a, b in intervals if a < b]
        ints.sort(key=lambda x: x[0])

        merged = [ints[0]]
        for lo, hi in ints[1:]:
            # If the next interval starts before the previous one ends, merge them
            if lo <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
            else:
                # Otherwise, it’s a separate interval
                merged.append((lo, hi))
        return merged


    # ---------------------------------------------------------------------
    # STEP 1: Gather all flagged frequency intervals
    # ---------------------------------------------------------------------
    # Your variable allflagvals is expected to look like this:
    # allflagvals = [[low1, high1], [low2, high2], ...] in Hz
    flag_ints = [(float(a), float(b)) for a, b in allflagvals]

    # ---------------------------------------------------------------------
    # STEP 2: Restrict flagged intervals to your two sidebands
    # ---------------------------------------------------------------------
    flag_in_band = flag_ints

    # ---------------------------------------------------------------------
    # STEP 3: Merge any overlapping flagged regions
    # ---------------------------------------------------------------------
    flag_merged = merge_intervals(flag_in_band)

    # ---------------------------------------------------------------------
    # STEP 4: Calculate how much frequency is flagged (in Hz)
    # ---------------------------------------------------------------------
    flagged_bw = sum(hi - lo for lo, hi in flag_merged)

    # ---------------------------------------------------------------------
    # STEP 5: Calculate total usable bandwidth (both sidebands)
    # ---------------------------------------------------------------------
    total_bw = sum(bhi - blo for blo, bhi in bands)

    # ---------------------------------------------------------------------
    # STEP 6: Compute percentage of total bandwidth flagged
    # ---------------------------------------------------------------------
    pct_total = 100.0 * flagged_bw / total_bw if total_bw > 0 else 0.0

    # ---------------------------------------------------------------------
    # STEP 7: (Optional) Compute per-sideband flagging percentages
    # ---------------------------------------------------------------------
    per_band = []
    for blo, bhi in bands:
        # Keep only portions of flagged intervals that fall within this band
        seg = merge_intervals([
            (max(lo, blo), min(hi, bhi))
            for lo, hi in flag_merged
            if min(hi, bhi) > max(lo, blo)
        ])
        bw = sum(hi - lo for lo, hi in seg)
        percent = 100.0 * bw / (bhi - blo)
        per_band.append(percent)

    # ---------------------------------------------------------------------
    # STEP 8: Print results in a clear format
    # ---------------------------------------------------------------------
    print(f"Total frequency coverage: {(total_bw / 1e9):.3f} GHz")
    print(f"Unique flagged frequency width: {(flagged_bw / 1e9):.3f} GHz")
    print(f"Total percentage flagged: {pct_total:.2f}%")
    print(f"Flagged per sideband:  LSB = {per_band[0]:.2f}%,  USB = {per_band[1]:.2f}%")




tb.close()


# In[10]:


#print flags
for i in totflags2:
    print(i)


# In[11]:


#calibrate

#science calibration phasecal is scan not int
#science calibration uses calflag strict

########################
import numpy as np
########################
from casatasks import flagmanager
from casatasks import flagdata
from casatasks import gaincal
from casatasks import fluxscale
from casatasks import setjy
from casatasks import bandpass
from casatasks import applycal



fc_prefix_cal=f'{namer}/images/clean/flag{flagnum}_cal{calnum}/caltables/'
fc_prefix=f'{namer}/images/clean/flag{flagnum}_cal{calnum}/'

for iflux in fluxfields:
    flux = iflux[1]
    for ibp in bpfields:
        bpcal= ibp[1]

        j=-1
        for itc_match in tc_match:
            pcal = itc_match[1]
            science_field = itc_match[0]
            found=0
            targets=[]
            for i in totflags:
                if i[0] in science_field:
                    k=-1
                    for ii in i:
                      k=k+1
                      if k==0:
                          continue
                      for iii in ii:
                        found=1
                        target= [iii[1] * 10**9,  iii[2] * 10**9]
                        targets.append(target)

            if found==0 and '5495' not in science_field:
                print(science_field)
                input('could not find')

            if not _is_selected_clean_detection(science_field):
                print(f'skipping clean target for date={date} flag={flagnum}: {science_field} does not match the manual chooser')
                continue

            breaker2=0
            for iut in unobserved_targets:
                if science_field==iut[0]:
                    breaker2=1
            if breaker2==1:
                print(f'no observed target for {science_field}')
                continue
            #only until after the checks can j increase.  regions are not evaluated for 4 of the 19 fields
            j=j+1

            #print(f'flux: {flux}\nbandpass:{bpcal}\nphase:{pcal}\nscience:{science_field}')
            calfields= ",".join([bpcal, pcal, flux])

            if flagger==True:        flagmanager(vis=myvis, mode='save', versionname='original')

            #edge flags
            rechunk = binsize
            chan_num = int(16384 / rechunk)
            edgechan_frac = 0.025
            edge_chan_low = int(np.floor(chan_num * edgechan_frac))
            edge_chan_high = int(np.floor(chan_num * (1. - edgechan_frac)))
            edgechan = "*:0~{0};{1}~{2}".format(edge_chan_low, edge_chan_high, chan_num-1)

            if flagger==True:        flagdata(vis=myvis, mode='manual', spw=edgechan, flagbackup=False)
            if flagger==True:        flagmanager(vis=myvis, mode='save', versionname='edge_flagging')


            if flagnum!=1:
                  for i in totflags2:
                      if i[0] in science_field:
                        flagdata(vis=myvis,mode='manual',spw=f'{i[1]}:{i[2]}~{i[3]}',flagbackup=False)
            if flagnum==3:
              flagdata(vis=myvis,mode='manual',spw='0',flagbackup=False)
              flagdata(vis=myvis,mode='manual',spw='1',flagbackup=False)
              flagdata(vis=myvis,mode='manual',spw='2',flagbackup=False)
              flagdata(vis=myvis,mode='manual',spw='3',flagbackup=False)
              flagdata(vis=myvis,mode='manual',spw='4',flagbackup=False)
              flagdata(vis=myvis,mode='manual',spw='5',flagbackup=False)
            if flagnum==4:
              flagdata(vis=myvis,mode='manual',spw='6',flagbackup=False)
              flagdata(vis=myvis,mode='manual',spw='7',flagbackup=False)
              flagdata(vis=myvis,mode='manual',spw='8',flagbackup=False)
              flagdata(vis=myvis,mode='manual',spw='9',flagbackup=False)
              flagdata(vis=myvis,mode='manual',spw='10',flagbackup=False)
              flagdata(vis=myvis,mode='manual',spw='11',flagbackup=False)


            #for i in totflags2:
            #    if i[0] in science_field:
            #      flagdata(vis=myvis,mode='manual',spw=f'{i[1]}:{i[2]}~{i[3]}',flagbackup=False)

            flagmanager(vis=myvis, mode='save', versionname='all_flagging')

            #options for standard:
            #standard = 'Perley-Butler 2017'
            #standard = 'Perley-Butler 2013'
            #standard = 'Perley-Butler 2010'
            #standard = Butler-JPL-Horizons 2012
            if namer =='230316_128':
              standard = 'Butler-JPL-Horizons 2012'
              standval=standard
            if namer == '220131_128':
              pass
              #standard = 'Perley-Butler 2017'
            #standard = manual
            #standard = fluxscale

            if namer=='230316_128':
                setjy(vis=myvis, 
                    field=flux, 
                    scalebychan=True,#scale by channel
                    standard=standval,
                    usescratch=True#save model data
                )

            if 0==1:#namer=='220131_128':
                setjy(vis=myvis, 
                    field=flux, 
                    scalebychan=True,#scale by channel
                    standard=standval,
                    usescratch=False#save model data
                )

            if namer=='220131_128':
                standval='manual'
                setjy(vis=myvis, 
                    field=flux, 
                    spw='',
                    fluxdensity=[1.98, 0, 0, 0],
                    scalebychan=True,#scale by channel
                    standard=standval,
                    usescratch=False#save model data
                )
            if inputs==True:        input(   'setjy done')


            refant = 'Ant8'
            min_solint = 'int'
            long_solint = 'inf'

            phaseshortgaincal_table = f'{fc_prefix_cal}bpself.gcal'
            if os.path.exists(phaseshortgaincal_table):
                os.system(f'rm -rf {phaseshortgaincal_table}')

            os.makedirs( phaseshortgaincal_table)
            gaincal(vis=myvis,
                    caltable=phaseshortgaincal_table,
                    field=bpcal,
                    spw="",
                    refant=refant, 
                    scan="",
                    calmode='p',
                    solint=min_solint,
                    minsnr=2.0, 
                    minblperant=3,
                    gaintable=[])
            if inputs==True:        input(   'phaseshortgaincal done')

            bpsolint = 'inf'
            bandpass_table = f'{fc_prefix_cal}bandpass.solnorm_true.bcal'
            if os.path.exists(bandpass_table):
                os.system(f'rm -rf {bandpass_table}')
            bandpass(vis=myvis,
                    caltable=bandpass_table,
                    bandtype='B', 
                    scan="",
                    field=bpcal, 
                    spw="",
                    refant=refant,
                    solint=bpsolint,
                    solnorm=False,
                    minblperant=3,
                    gaintable=[phaseshortgaincal_table])
            if inputs==True:        input(   'bandpass done')

            ########################
            #phase cal is for atmospheric and systematic changes in time
            #to boost S/N for phase gain we combine all of SPW per side-band. 
            #We can do this because adjacent SPW have similair phase variation

            #we use per-int phase and combine spws across sidebands with append='True'
            #Use short soltime as well as long soltime for each table

            gain_phase_int_table = f'{fc_prefix_cal}intphase_combinespw.gcal'
            if os.path.exists(gain_phase_int_table):
                os.system(f'rm -rf {gain_phase_int_table}')
            # LSB
            if flagnum!=3 and flagnum!=4:
                gaincal(vis=myvis,
                        caltable=gain_phase_int_table,
                        field=calfields,
                        refant=refant,
                        combine='spw',
                        spw='0~5',
                        calmode='p',
                        solint=min_solint,
                        minsnr=2.0,
                        minblperant=3,
                        gaintable=[bandpass_table])
                if inputs==True:        input(   'gain_phase_int_table lsb done')
                # USB
                gaincal(vis=myvis,
                        caltable=gain_phase_int_table, 
                        append=True,
                        field=calfields, 
                        refant=refant,
                        combine='spw',
                        spw='6~11',
                        calmode='p',
                        solint=min_solint,
                        minsnr=2.0,
                        minblperant=3,
                        gaintable=[bandpass_table])

            if flagnum==3:
                #USB
                gaincal(vis=myvis,
                        caltable=gain_phase_int_table, 
                        field=calfields, 
                        refant=refant,
                        combine='spw',
                        spw='6~8',
                        calmode='p',
                        solint=min_solint,
                        minsnr=2.0,
                        minblperant=3,
                        gaintable=[bandpass_table])           
                gaincal(vis=myvis,
                        caltable=gain_phase_int_table, 
                        append=True,
                        field=calfields, 
                        refant=refant,
                        combine='spw',
                        spw='9~11',
                        calmode='p',
                        solint=min_solint,
                        minsnr=2.0,
                        minblperant=3,
                        gaintable=[bandpass_table])     
            if flagnum==4:
                #LSB
                gaincal(vis=myvis,
                        caltable=gain_phase_int_table, 
                        field=calfields, 
                        refant=refant,
                        combine='spw',
                        spw='0~2',
                        calmode='p',
                        solint=min_solint,
                        minsnr=2.0,
                        minblperant=3,
                        gaintable=[bandpass_table])      
                gaincal(vis=myvis,
                        caltable=gain_phase_int_table, 
                        append=True,
                        field=calfields, 
                        refant=refant,
                        combine='spw',
                        spw='3~5',
                        calmode='p',
                        solint=min_solint,
                        minsnr=2.0,
                        minblperant=3,
                        gaintable=[bandpass_table])      


            if inputs==True:        input(   'gain_phase_int_table usb done')

            gain_phase_scan_table = f'{fc_prefix_cal}scanphase_combinespw.gcal'
            if os.path.exists(gain_phase_scan_table):
                os.system(f'rm -rf {gain_phase_scan_table}')
            if flagnum!=3 and flagnum!=4:
                # LSB
                gaincal(vis=myvis,
                        caltable=gain_phase_scan_table,
                        field=calfields,
                        refant=refant,
                        combine='spw',
                        spw='0~5',
                        calmode='p',
                        solint=long_solint,
                        minsnr=2.0,
                        minblperant=3,
                        gaintable=[bandpass_table])
                if inputs==True:        input(   'gain_phase_scan_table lsb done')
                # USB
                gaincal(vis=myvis,
                        caltable=gain_phase_scan_table, 
                        append=True,
                        field=calfields, 
                        refant=refant,
                        combine='spw',
                        spw='6~11',
                        calmode='p',
                        solint=long_solint,
                        minsnr=2.0,
                        minblperant=3,
                        gaintable=[bandpass_table])

            if flagnum==3:
                # USB
                gaincal(vis=myvis,
                        caltable=gain_phase_scan_table,
                        field=calfields,
                        refant=refant,
                        combine='spw',
                        spw='6~8',
                        calmode='p',
                        solint=long_solint,
                        minsnr=2.0,
                        minblperant=3,
                        gaintable=[bandpass_table])    

                gaincal(vis=myvis,
                        caltable=gain_phase_scan_table,
                        append=True,
                        field=calfields,
                        refant=refant,
                        combine='spw',
                        spw='9~11',
                        calmode='p',
                        solint=long_solint,
                        minsnr=2.0,
                        minblperant=3,
                        gaintable=[bandpass_table])   

            if flagnum==4:  
                #LSB  
                gaincal(vis=myvis,
                        caltable=gain_phase_scan_table,
                        field=calfields,
                        refant=refant,
                        combine='spw',
                        spw='0~2',
                        calmode='p',
                        solint=long_solint,
                        minsnr=2.0,
                        minblperant=3,
                        gaintable=[bandpass_table])     
                gaincal(vis=myvis,
                        caltable=gain_phase_scan_table,
                        append=True,
                        field=calfields,
                        refant=refant,
                        combine='spw',
                        spw='3~5',
                        calmode='p',
                        solint=long_solint,
                        minsnr=2.0,
                        minblperant=3,
                        gaintable=[bandpass_table])   



            if inputs==True:        input(   'gain_phase_scan_table usb done')

            amp_solint=long_solint
            if flagnum==1 or flagnum==2:
                  this_spwmap = [0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 6 ,6]
            elif flagnum==3:
                  this_spwmap = [6, 6, 6, 9, 9, 9]
            elif flagnum==4:
                  this_spwmap = [0, 0, 0, 4, 4, 4]
            gain_amp_scan_table = f'{fc_prefix_cal}amp.gcal'
            #long amplitude and short phase corrections
            if os.path.exists(gain_amp_scan_table):
                os.system(f'rm -rf {gain_amp_scan_table}')
            gaincal(vis=myvis,
                    caltable=gain_amp_scan_table,
                    field=calfields, 
                    spw="",
                    refant=refant,
                    combine='',
                    spwmap=[[],this_spwmap],
                    calmode='a', 
                    solint=amp_solint, 
                    minsnr=3.0, 
                    minblperant=3,
                    gaintable=[bandpass_table, gain_phase_int_table])# amp scan table has not yet had short_int applied until now
            if inputs==True:        input(   'gain_amp_scan_table done')
            ##########################

            # Apply to all calibrators except the flux calibrator itself.
            transfer_fields = [bpcal, pcal]
            fluxboot_table = f'{fc_prefix_cal}flux.cal'
            if os.path.exists(fluxboot_table):
                os.system(f'rm -rf {fluxboot_table}')
            fluxresults = fluxscale(vis=myvis,
                                    caltable=gain_amp_scan_table,
                                    refspwmap=[-1],
                                    transfer=transfer_fields,
                                    fluxtable=fluxboot_table,
                                    reference=flux,
                                    fitorder=1,
                                    incremental=False
                                    )
            if inputs==True:        input(   'fluxscale done')

            if flagger==True:        flagmanager(vis=myvis, mode='save', versionname='beforeapplycal')
            # The following shows how to restore flags
            # flagmanager(vis=myvis, mode='restore', versionname='beforeapplycal')

            clearcal(vis=myvis)

            applycal(vis=myvis,field=science_field,
                    spw="",
                    gaintable=[bandpass_table,gain_phase_scan_table,fluxboot_table],#science table for phasecal is scan not int
                    interp=['nearest','linear','nearest'],
                    spwmap=[[], this_spwmap, []],
                    gainfield=[bpcal,pcal, pcal],
                    flagbackup=False, calwt=False,
                    applymode='calflagstrict')#Another feature of for science calibration
            sciloc=f'{fc_prefix}fields/science/{science_field}_test'
            if os.path.exists(sciloc):
                os.system(f'rm -rf {sciloc}')
            os.makedirs(sciloc)
            tclean_safe(vis=myvis, 
                field=science_field, 
                imagename=f'{sciloc}',
                imsize=imsizer, 
                cell='1arcsec', 
                specmode='mfs',
                #weighting='briggs', 
                #robust=2., 
                niter=sciniter,
                mask=regions[j])



            applycal(vis=myvis,field=pcal,
                    spw="",
                    gaintable=[bandpass_table,gain_phase_int_table,fluxboot_table],
                    interp=['nearest','linear','nearest'],
                    spwmap=[[], this_spwmap, []],
                    gainfield=[bpcal,pcal, pcal],
                    flagbackup=False, calwt=False)

            phaseloc=f'{fc_prefix}fields/phasecal/{pcal}_test'
            if os.path.exists(phaseloc):
                os.system(f'rm -rf {phaseloc}')
            os.makedirs(phaseloc)
            tclean_safe(vis=myvis, 
                field=pcal, 
                imagename=phaseloc,
                imsize=imsizer, 
                cell='1arcsec', 
                specmode='mfs',
                #weighting='briggs', 
                robust=0., 
                niter=sciniter,
                usemask='pb')
            #imview(f'{phaseloc}.image')

        ###############################
        applycal(vis=myvis,field=bpcal,
                spw="",
                gaintable=[bandpass_table,
                            gain_phase_int_table,
                            fluxboot_table], # This order for tables is reflected in `interp`, `spwmap`, and `gainfield`
                interp=['nearest',
                        'linear',
                        'nearest'], # for each field, this is how to average/interpolate the solutions in "time,freq"
                spwmap=[[], 
                        this_spwmap, 
                        []],
                gainfield=[bpcal, 
                            bpcal, 
                            bpcal],
                flagbackup=False,
                calwt=False)
        bploc=f'{fc_prefix}fields/bpcal/{bpcal}_test'
        if os.path.exists(bploc):
            os.system(f'rm -rf {bploc}')
        os.makedirs(bploc)

        tclean_safe(vis=myvis, field=bpcal, imagename=bploc,
            imsize=imsizer, cell='1arcsec', specmode='mfs',
            #weighting='briggs', 
            #robust=0., 
            niter=sciniter)
        #imview(f'{bploc}.image')

        ###################
        applycal(vis=myvis,field=flux,
                spw="",
                gaintable=[bandpass_table,
                            gain_phase_int_table,
                            fluxboot_table], # This order for tables is reflected in `interp`, `spwmap`, and `gainfield`
                interp=['nearest',
                        'linear',
                        'nearest'], # for each field, this is how to average/interpolate the solutions in "time,freq"
                spwmap=[[], 
                        this_spwmap, 
                        []],
                gainfield=[bpcal, 
                            flux, 
                            flux],
                flagbackup=False,
                calwt=False)

        fluxloc=f'{fc_prefix}fields/fluxcal/{flux}_test'
        if os.path.exists(fluxloc):
            os.system(f'rm -rf {fluxloc}')
        os.makedirs(fluxloc)

        tclean_safe(vis=myvis, field=flux, imagename=fluxloc,
            imsize=imsizer, cell='1arcsec', specmode='mfs',
            #weighting='briggs', 
            #robust=0., 
            niter=sciniter)
        #imview(f'{fluxloc}.image')

        ###################




# In[12]:


print(fluxfields)
print(bpfields)


# In[13]:


for i in totflags2:
    print(f'manual ,spw={i[1]}:{i[2]}~{i[3]},flagbackup=False')
print(chooser)


# In[14]:


print(flagnum)
print(chooser)

# move this run's products into an iteration-specific folder
import os
import shutil

iter_dir = f"{namer}/images/clean/nitercal_{sciniter}"
src_dir = f"{namer}/images/clean/flag{flagnum}_cal{calnum}"

os.makedirs(iter_dir, exist_ok=True)

dst_dir = os.path.join(iter_dir, f"flag{flagnum}_cal{calnum}")

if os.path.exists(src_dir):
            if os.path.exists(dst_dir):
                        shutil.rmtree(dst_dir)
            shutil.move(src_dir, dst_dir)
