#!/usr/bin/env bash

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python_bin="${PYTHON_BIN:-python3}"

for date in 220131 230316
do
   for flagnum in 1 2 3 4
   do
      for sciniter in 0
      do
         "${python_bin}" "${script_dir}/calibrate_sma_dirty.py" --flagnum "$flagnum" --sciniter "$sciniter" --date "$date"
      done
   done
done
