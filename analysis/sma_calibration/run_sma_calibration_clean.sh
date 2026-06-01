#!/usr/bin/env bash

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python_bin="${PYTHON_BIN:-python3}"

for date in 230316 220131
do
   for sciniter in {1..100} #10 100 300
   do
      # Visit every requested flag state for a given clean iteration before
      # advancing to the next niter value.
      for flagnum in 2 3 4
      do
         "${python_bin}" "${script_dir}/calibrate_sma_clean.py" --flagnum "$flagnum" --sciniter "$sciniter" --date "$date"
      done
   done
done
