#!/usr/bin/env bash

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
fundir2="${repo_root}/analysis/sma_calibration"
python_bin="${PYTHON_BIN:-python3}"

if ! command -v "${python_bin}" >/dev/null 2>&1; then
      echo "ERROR: ${python_bin} not found. Set PYTHON_BIN to a CASA-capable python interpreter."
      exit 1
fi

if ! "${python_bin}" - <<'PY' >/dev/null 2>&1
import casatasks
import casatools
PY
then
      echo "ERROR: ${python_bin} cannot import casatasks/casatools."
      echo "Activate the CASA environment you used for calibration, or rerun with:"
      echo "  PYTHON_BIN=/path/to/casa-python bash run_continuum_fits_clean.sh"
      exit 1
fi

nitermover=1
niterarray=(0)
datearray=(220131 230316)
flagarray=(2)

for date in "${datearray[@]}"
do
      for flagnum in "${flagarray[@]}"
      do
            for nitercal in "${niterarray[@]}"
            do
                  smapath="${fundir2}/${date}_128/images/nitercal_${nitercal}/flag${flagnum}_cal1/fields/science"
                  if [[ ! -d "${smapath}" ]]; then
                        echo "Skipping date=${date} flag=${flagnum} nitercal=${nitercal}: missing ${smapath}"
                        continue
                  fi

                  "${python_bin}" "${script_dir}/get_continuum_params.py" \
                        --flagnum "$flagnum" \
                        --calnum 1 \
                        --mode fixtochosenfreq \
                        --datestr "$date" \
                        --chooser 0 \
                        --choose_science NGC4258 \
                        --printimages 0 \
                        --writesummary 1 \
                        --nofixmultval 10 \
                        --printbadfix 0 \
                        --telechooser 0 \
                        --telechoice SMA \
                        --moveniter "$nitermover" \
                        --niternum "$nitercal" \
                        --smapath "${smapath}"
            done
      done
done
