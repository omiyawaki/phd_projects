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
yr_end="1980"
mn_end="01"
dy_end="31"

cwd=$(pwd) # save current working directory

#if ls $cwd/${rea}_sfcWind_${yr_begin}_${yr_end}.nc 1> /dev/null 2>&1; then # check if data is already there
#    echo "tend was already created. Skipping..."
#else
for yr in $(seq ${yr_begin} 1 ${yr_end}); do
    for mn in $(seq -w 01 1 01); do
        for dy in $(seq -w 01 1 31); do

            if [ ! -f *${yr}${mn}${dy}*.nc4 ]; then
                echo "File doesn't exist for ${yr}${mn}${dy}, skipping..."
                continue
            fi

            full_uas=${cwd}/*${yr}${mn}${dy}*.nc4
            full_vas=${cwd}/*${yr}${mn}${dy}*.nc4
            full_sfcWind=${cwd}/${rea}_sfcWind_${yr}${mn}${dy}.nc

            # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_tend_mon.py ${full_ps} ${full_ta} ${full_zg} ${full_hus} ${full_tend} ${full_tas} ${full_orog} ${full_huss}
            python ${cwd}/make_sfcWind.py ${full_uas} ${full_vas} ${full_sfcWind}

        done
    done
done
#fi
