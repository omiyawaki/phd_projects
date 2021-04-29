#!/bin/sh

declare -a vtypes=("rad" "hydro" "stf" "temp" "srfc" "zg" "lfrac" "orog") # variable type
#declare -a vtypes=("temp_ml") # variable type
yr_start="1980"
mon_start="1"
day_start="1"
yr_end="2005"
mon_end="12"
day_end="30"

cwd=$(pwd)

for vtype in ${vtypes[@]}; do
    cd $vtype

    filename=era5_${vtype}_${yr_start}_${yr_end}

    if test -f "$filename.ymonmean.nc"; then
        echo "$filename.ymonmean.nc exists. Skipping..."
    else
        cdo seldate,$yr_start-$mon_start-$day_start,$yr_end-$mon_end-$day_end era5_${vtype}_1979_2005.nc $filename.nc

        if [[ $vtype != "tend" ]]; then
            cdo ymonmean $filename.nc $filename.ymonmean.nc
        fi
    fi

    cd $cwd

done
