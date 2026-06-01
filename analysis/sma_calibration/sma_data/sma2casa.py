#!/usr/bin/env python3

# In the command line the last / must be omitted.
# The output file should look like _master.ms.
import os
import argparse
from pathlib import Path
from pyuvdata import UVData

def process_folder(folder):
    folder=folder.replace("/","")
    cwd = Path(__file__).resolve().parent
    UV = UVData()

    UV.read_mir(f"{folder}")

    UV.write_ms(str(cwd / f"{folder}_master.ms"))
    UV.write_uvfits(str(cwd / f"{folder}_master.uvfits"))

    uvfits_file = cwd / f"{folder}_master.uvfits"
    os.remove(uvfits_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process a folder with pyuvdata.')
    parser.add_argument('folder', type=str, help='Folder name to process with pyuvdata.')

    args = parser.parse_args()

    process_folder(args.folder)

