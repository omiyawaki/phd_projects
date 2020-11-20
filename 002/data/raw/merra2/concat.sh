#!/bin/sh

declare -a vtypes=("rad" "hydro" "stf" "temp" "srfc" "zg") # variable type
yr_start="2000"
mon_start="3"
day_start="1"
yr_end="2018"
mon_end="2"
day_end="1"

cwd=$(pwd)

for vtype in ${vtypes[@]}; do
    cd $vtype

    filename=merra2_${vtype}_${yr_start}_${yr_end}

    if test -f "$filename.ymonmean.nc"; then
        echo "$filename.ymonmean.nc exists. Skipping..."
    else
        cdo -b F64 mergetime MERRA2_*.nc4 $filename.merge.nc
        cdo seldate,$yr_start-$mon_start-$day_start,$yr_end-$mon_end-$day_end $filename.merge.nc $filename.nc
        rm $filename.merge.nc

        if [[ $vtype != "tend" ]]; then
            cdo ymonmean $filename.nc $filename.ymonmean.nc
        fi
    fi

    cd $cwd

done
