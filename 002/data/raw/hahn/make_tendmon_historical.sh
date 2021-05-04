#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python
source activate netcdf

clim="ControlSOM1850"

cwd=$(pwd) # save current working directory

#if ls $cwd/${clim}.tendmon.nc 1> /dev/null 2>&1; then # check if data is already there
#    echo "tend was already created. Skipping..."
#else
    full_ps=${cwd}/lapserateclima/${clim}.PS.nc
    full_tas=${cwd}/lapserateclima/${clim}.TREFHT.nc
    full_huss=${cwd}/lapserateclima/${clim}.Q.nc
    full_orog=${cwd}/lapserateclima/${clim}.PS.nc
    full_ta=${cwd}/lapserateclima/${clim}.T.nc
    full_zg=${cwd}/lapserateclima/${clim}.Z3.nc
    full_hus=${cwd}/lapserateclima/${clim}.Q.nc
    full_tend=${cwd}/${clim}.tendmon.nc

    # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_tend_mon.py ${full_ps} ${full_ta} ${full_zg} ${full_hus} ${full_tend} ${full_tas} ${full_orog} ${full_huss}
    python ${cwd}/make_tend_mon.py ${full_ps} ${full_ta} ${full_zg} ${full_hus} ${full_tend} ${full_tas} ${full_orog} ${full_huss}
#fi
