#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python
source activate netcdf

rea="jra55"
yr_begin="1980"
mn_begin="01"
dy_begin="01"
yr_end="2005"
mn_end="12"
dy_end="31"

cwd=$(pwd) # save current working directory
arch_path="/project2/tas1/abacus/data1/tas/archive/Reanalysis/JRA55"

if ls $cwd/tendmon/${rea}_tend_${yr_begin}_${yr_end}.nc 1> /dev/null 2>&1; then # check if data is already there
    echo "tend was already created. Skipping..."
else
    for yr in $(seq $yr_begin $yr_end); do
        full_uas=${arch_path}/uv10m/*ugrd.${yr}*.nc
        full_vas=${arch_path}/uv10m/*vgrd.${yr}*.nc
        full_sfcWind=${cwd}/${rea}_sfcWind_${yr}.nc

        echo "Computing surface wind for ${yr}..."
        # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_sfcWind.py ${full_uas} ${full_vas} ${full_sfcWind}
        python ${cwd}/make_sfcWind.py ${full_uas} ${full_vas} ${full_sfcWind}
        echo -e "Done.\n"
    done

fi
