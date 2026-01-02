#in the command line the last / must be ommitted
#the output file should look like _master.ms

#!/opt/anaconda3/envs/pyuvdata/bin/python
import os
import argparse
from pyuvdata import UVData

def process_folder(folder):
    folder=folder.replace("/","")
    cwd = '/home/peacesea/Burridge2025_WM_AGN/Public/analysis/sma_calibration/sma_data'
    UV = UVData()

    UV.read_mir(f"{folder}")

    UV.write_ms(cwd + f"/{folder}_master.ms")
    UV.write_uvfits(cwd + f"/{folder}_master.uvfits")

    uvfits_file = f"{folder}_master.uvfits"
    os.remove(uvfits_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process a folder with pyuvdata.')
    parser.add_argument('folder', type=str, help='Folder name to process with pyuvdata.')

    args = parser.parse_args()

    process_folder(args.folder)

