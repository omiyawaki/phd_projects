#!/bin/sh

varnames=("ftoa")

# clims=("rp000130" "rp000131" "rp000134" "rp000135") # climate name
clims=("rp000135") # climate name

mean=""
yr_span="0001_0039"

# save path to current directory
cwd=$(pwd)

for clim in ${clims[@]}; do

    echo ${clim}

    cd ${cwd}/${clim}

    for varname in ${varnames[@]}; do

        echo ${varname}

        filename="${varname}_${clim}_${yr_span}${mean}"

        # create ymonmean mean file if it doesn't exist yet
        # if [ -f "${filename}.ymonmean-20.nc" ]; then
        #     echo "ymonmean mean already taken, skipping..."
        cdo -yearmean ${filename}.nc ${filename}.yearmean.nc 
        cdo -zonmean ${filename}.yearmean.nc ${filename}.yearmean.zonmean.nc 
        cdo -mermean ${filename}.yearmean.zonmean.nc ${filename}.gmean.nc 
        rm ${filename}.yearmean.zonmean.nc

    done # varnames
done # clims
