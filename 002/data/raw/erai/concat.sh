#!/bin/sh

# $1 is vtype
# $2 is yr_start
# $3 is yr_end
# $4 is mon_start
# $5 is mon_end
# $6 is day_start
# $7 is day_end

filename=erai_$1_$2_$3

if test -f "$filename.ymonmean.nc"; then
    echo "$filename.ymonmean.nc exists. Skipping..."
else
    if test -f "$filename.nc"; then
        echo "$filename.nc exists. Skipping mergetime..."
    else
        cdo -b F64 mergetime interim_$1_*.nc $filename.merge.nc
        cdo seldate,$2-$4-$6,$3-$5-$7 $filename.merge.nc $filename.nc
        rm $filename.merge.nc
    fi

    if [[ $1 != "tend" ]]; then
        cdo ymonmean $filename.nc $filename.ymonmean.nc
        rm $filename.nc
    fi
fi
