#!/bin/sh

# varnames=("srafs" "srad0u" "sradsu" "fsfc" "tsurf" "siced" "albedo" "srad0d" "trads" "srads" "srad0" "trad0" "ahfl" "ahfs" "seaice" "temp2" "precip" "aprl" "aprc" "aps" "tend" "ra" "stf" "ra" "stgadv" "adv" "r1") # list of GCM variables that we want to process
# varnames=("tsurf" "trads" "srads" "ahfl" "ahfs" "fsfc") # list of GCM variables that we want to process
varnames=("srafs")
# varnames=("srad0d" "albedo" "sradsu" "seaice" "siced" "temp2" "tsurf" "trads" "srads" "ahfl" "ahfs" "fsfc")

clims=("rp000190" "rp000190a" "rp000190b" "rp000190c") # climate name
# clims=("rp000190a") # climate name

mean=""
yr_span="0001_0009"

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
        cdo seldate,0005-01-01,0009-12-31 ${filename}.nc ${filename}.sel.nc
        cdo -shifttime,-15days ${filename}.sel.nc ${filename}.shifttime.nc 
        cdo -ymonmean ${filename}.shifttime.nc ${filename}.ymonmean-5.nc 

        rm ${filename}.sel.nc
        rm ${filename}.shifttime.nc

    done # varnames
done # clims
