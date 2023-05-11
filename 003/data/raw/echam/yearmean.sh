#!/bin/sh

varnames=("tradsu" "sradsu" "fsfc" "tsurf" "siced" "albedo" "srad0d" "trads" "srads" "srad0" "trad0" "ahfl" "ahfs" "seaice" "temp2" "precip" "aprl" "aprc" "aps" "tend" "ra" "stf" "ra" "stgadv" "adv" "r1") # list of GCM variables that we want to process
# varnames=("tsurf" "trads" "srads" "ahfl" "ahfs" "fsfc") # list of GCM variables that we want to process
# varnames=("tsurf" "temp2")
# varnames=("srad0d" "albedo" "sradsu" "seaice" "siced" "temp2" "tsurf" "trads" "srads" "ahfl" "ahfs" "fsfc")

clims=("rp000201") # climate name

mean=""
yr_span="0040_0083"

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
        cdo -shifttime,-15days ${filename}.nc ${filename}.shifttime.nc 
        cdo -yearmean ${filename}.shifttime.nc ${filename}.yearmean.nc 

        rm ${filename}.shifttime.nc

    done # varnames
done # clims
