#!/bin/sh

co2_vmrs=("300" "400" "500" "600" "700" "800" "900" "1000" "1100" "1200")
# co2_vmrs=("300")
exp="co2vmr"
var="T"

for co2 in ${co2_vmrs[@]}; do
    file="./proc/$var/${exp}_${co2}.$var"
    cdo -O -seasmean -selseas,DJF $file.nc ${file}.djfmean.nc.tmp
    cdo seldate,0002-12-31,0003-01-31 $file.djfmean.nc.tmp $file.djfmean.nc
    rm $file.djfmean.nc.tmp
done
