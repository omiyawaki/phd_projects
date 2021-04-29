#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python
source activate netcdf

rea="merra2c"
yr_begin="1980"
mn_begin="01"
dy_begin="01"
yr_end="2005"
mn_end="12"
dy_end="31"

cwd=$(pwd) # save current working directory

mkdir -p $cwd/tendmon # make model directory in processed data folder if it doesn't exist yet

if ls $cwd/tendmon/${rea}_tend_${yr_begin}_${yr_end}.nc 1> /dev/null 2>&1; then # check if data is already there
    echo "tend was already created. Skipping..."
else
    full_ps=${cwd}/srfc/${rea}_srfc_${yr_begin}_${yr_end}.nc
    full_tas=${cwd}/srfc/${rea}_srfc_${yr_begin}_${yr_end}.nc
    full_huss=${cwd}/srfc/${rea}_srfc_${yr_begin}_${yr_end}.nc
    full_orog=${cwd}/srfc/${rea}_orog_${yr_begin}_${yr_end}.nc
    full_ta=${cwd}/temp/${rea}_temp_${yr_begin}_${yr_end}.nc
    full_zg=${cwd}/zg/${rea}_zg_${yr_begin}_${yr_end}.nc
    full_hus=${cwd}/hus/${rea}_hus_${yr_begin}_${yr_end}.nc
    full_tend=${cwd}/tendmon/${rea}_tendmon_${yr_begin}_${yr_end}.nc

    srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_tend_mon.py ${full_ps} ${full_ta} ${full_zg} ${full_hus} ${full_tend} ${full_tas} ${full_orog} ${full_huss}
    # python ${cwd}/make_tend_mon.py ${full_ps} ${full_ta} ${full_zg} ${full_hus} ${full_tend} ${full_tas} ${full_orog} ${full_huss}
fi
