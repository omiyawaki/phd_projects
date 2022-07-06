#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

co2_vmrs=("300" "400" "500" "600" "700" "800" "900" "1000" "1100" "1200")
# co2_vmrs=("300")
exp="qco2vmr"
var="TMQ"

mkdir ./proc/$var

for co2 in ${co2_vmrs[@]}; do
    ifile="./proc/merged/${exp}_${co2}.nc"
    ofile="./proc/$var/${exp}_${co2}.$var.nc"

    cdo selvar,$var ${ifile} ${ofile}
done
