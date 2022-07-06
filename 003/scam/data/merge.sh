#!/bin/sh

co2_vmrs=("300" "400" "500" "600" "700" "800" "900" "1000" "1100" "1200")
exp="qco2vmr"

mkdir ./proc/merged

for co2 in ${co2_vmrs[@]}; do
    ifile=$(ls ./raw/${exp}_$co2.FSCAM.padv/atm/hist/*h0*.nc)
    ofile="./proc/merged/${exp}_${co2}.nc"
    cdo -O mergetime $ifile ${ofile}.tmp
    cdo shifttime,-15day ${ofile}.tmp ${ofile}
    rm ${ofile}.tmp
done
