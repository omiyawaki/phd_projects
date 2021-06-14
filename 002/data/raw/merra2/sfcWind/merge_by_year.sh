#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python
source activate netcdf

rea="merra2"
yr_begin="1981"
mn_begin="01"
dy_begin="01"
yr_end="2005"
mn_end="12"
dy_end="31"

cwd=$(pwd) # save current working directory

#if ls $cwd/${rea}_sfcWind_${yr_begin}_${yr_end}.nc 1> /dev/null 2>&1; then # check if data is already there
#    echo "tend was already created. Skipping..."
#else
for yr in $(seq ${yr_begin} 1 ${yr_end}); do
    infiles=$(ls *${yr}*.nc)
    echo "Merging data for ${yr}..."
    cdo mergetime ${infiles} ${rea}_sfcWind_${yr}.nc
    echo -e "Done.\n"

    rm ${infiles}
done

echo "Merging all data..."
cdo mergetime *${1980..2005}.nc ${rea}_sfcWind_1980_2005.nc
cdo ymonmean ${rea}_sfcWind_1980_2005.nc ${rea}_sfcWind_1980_2005.ymonmean.nc
echo -e "Done.\n"

#fi
