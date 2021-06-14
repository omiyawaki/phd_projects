#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python
source activate netcdf

rea="merra2"
yr_begin="1980"
mn_begin="01"
dy_begin="01"
yr_end="2005"
mn_end="12"
dy_end="31"

cwd=$(pwd) # save current working directory

#if ls $cwd/${rea}_sfcWind_${yr_begin}_${yr_end}.nc 1> /dev/null 2>&1; then # check if data is already there
#    echo "tend was already created. Skipping..."
#else
    full_uas=${cwd}/${rea}_uvas_${yr_begin}_${yr_end}.nc
    full_vas=${cwd}/${rea}_uvas_${yr_begin}_${yr_end}.nc
    full_sfcWind=${cwd}/${rea}_sfcWind_${yr_begin}_${yr_end}.nc

    # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_tend_mon.py ${full_ps} ${full_ta} ${full_zg} ${full_hus} ${full_tend} ${full_tas} ${full_orog} ${full_huss}
    python ${cwd}/make_sfcWind.py ${full_uas} ${full_vas} ${full_sfcWind}
#fi
