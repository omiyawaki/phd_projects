#!/bin/sh

varnames=("trads" "srads" "srad0" "trad0" "ahfl" "ahfs" "seaice" "temp2" "precip" "aprl" "aprc" "aps" "tend" "ra" "stf" "ra" "stgadv" "adv" "r1" "trafs" "srafs" "sraf0" "traf0" "siced") # list of GCM variables that we want to process
# varnames=("trafs" "srafs" "sraf0" "traf0") # list of GCM variables that we want to process
# varnames=("trads" "srads" "srad0" "trad0" "ahfl" "ahfs" "seaice" "temp2" "precip" "aprl" "aprc" "tend") # list of GCM variables that we want to process
# varnames=("seaice" "siced")
# varnames=("ahfres")

# clims=("rp000190") # climate name
# yr_span="0001_0009"
mean=""

# clims=("rp000134") # climate name
# yr_span="0001_0049"
# mean=""

# clims=("rp000130") # climate name
# yr_span="0001_0079"
# mean=""

# clims=("rp000184") # climate name
# clims=("rp000184" "rp000185" "rp000188" "rp000189") # climate name
# yr_span="0040_0139"
# mean=""

clims=("rp000134") # climate name
# clims=("rp000130" "rp000131" "rp000134" "rp000135") # climate name
yr_span="0020_0039"
# mean=".ymonmean-20"

# save path to current directory
cwd=$(pwd)

for clim in ${clims[@]}; do

    echo ${clim}

    cd ${cwd}/${clim}

    for varname in ${varnames[@]}; do

        echo ${varname}

        filename="${varname}_${clim}_${yr_span}${mean}"

        # create DJF mean file if it doesn't exist yet
        # if [ -f "${filename}.djfmean.nc" ]; then
        #     echo "DJF mean already taken, skipping..."
        if [[ $mean == *"ymonmean-20"* ]]; then
            cdo -selseas,DJF ${filename}.nc ${filename}.djfsel.nc 
            cdo -timmean ${filename}.djfsel.nc ${filename}.djfmean.nc 
            rm ${filename}.djfsel.nc 
        else
            cdo -shifttime,-15days ${filename}.nc ${filename}.shifttime.nc 
            cdo -seasmean -selseas,DJF ${filename}.shifttime.nc ${filename}.djfmean.nc 
            rm ${filename}.shifttime.nc
        fi

    done # varnames
done # clims
