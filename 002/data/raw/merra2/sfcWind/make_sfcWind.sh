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
arch_path="/project2/tas1/abacus/data1/tas/archive/Reanalysis/MERRA2"

#if ls $cwd/${rea}_sfcWind_${yr_begin}_${yr_end}.nc 1> /dev/null 2>&1; then # check if data is already there
#    echo "tend was already created. Skipping..."
#else
for yr in $(seq ${yr_begin} 1 ${yr_end}); do
    for mn in $(seq -w 01 1 12); do
        for dy in $(seq -w 01 1 31); do

            if [ ! -f ${arch_path}/uv10m/*${yr}${mn}${dy}*.nc ]; then
                echo -e "File doesn't exist for ${yr}-${mn}-${dy}, skipping...\n"
                continue
            fi

            full_uas=${arch_path}/uv10m/*${yr}${mn}${dy}*.nc
            full_vas=${arch_path}/uv10m/*${yr}${mn}${dy}*.nc
            full_sfcWind=${cwd}/${rea}_sfcWind_${yr}${mn}${dy}.nc

            echo "Computing surface wind for ${yr}-${mn}-${dy}..."
            # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_tend_mon.py ${full_ps} ${full_ta} ${full_zg} ${full_hus} ${full_tend} ${full_tas} ${full_orog} ${full_huss}
            python ${cwd}/make_sfcWind.py ${full_uas} ${full_vas} ${full_sfcWind}
            echo -e "Done.\n"

        done
    done
done
#fi
