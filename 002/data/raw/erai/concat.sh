#!/bin/sh

filename=erai_$1_$2_$3

if test -f "$filename.nc"; then
    echo "$filename.nc exists. Skipping..."
else
    cdo -b F64 mergetime interim_$1_*.nc $filename.nc
fi

if test -f "$filename.ymonmean.nc"; then
    echo "$filename.ymonmean.nc exists. Skipping..."
else
    cdo ymonmean $filename.nc $filename.ymonmean.nc
    rm $filename.nc
fi
