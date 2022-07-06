#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

co2_vmrs=("300" "400" "500" "600" "700" "800" "900" "1000" "1100" "1200")
# co2_vmrs=("300")
exp="qco2vmr"

mkdir ./proc/ra

for co2 in ${co2_vmrs[@]}; do
    ifile="./proc/merged/${exp}_${co2}.nc"
    ofile="./proc/ra/${exp}_${co2}.ra.nc"

    python make_ra.py ${ifile} ${ofile}
done
