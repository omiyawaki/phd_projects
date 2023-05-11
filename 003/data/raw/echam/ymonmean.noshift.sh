#!/bin/sh

# varnames=("sni" "xlvi" "xivi" "aclcov" "prfrac" "ftoa" "twtend" "titend" "ttend" "stend" "sn" "ahfres" "tsi" "tsw" "ahfcon" "trafs" "sraf0" "traf0" "srafs" "srad0d" "srad0u" "fsfc" "tsurf" "sradsu" "albedo" "trads" "srads" "srad0" "trad0" "ahfl" "ahfs" "siced" "seaice" "temp2" "precip" "aprl" "aprc" "aps" "tend" "ra" "stf" "ra" "stgadv" "adv" "r1") # list of GCM variables that we want to process
# varnames=("trafs" "srafs" "sraf0" "traf0") # list of GCM variables that we want to process
# varnames=("smelt" "snmelt")
varnames=("cond")

# clims=("rp000130" "rp000131" "rp000134" "rp000135") # climate name
# clims=("rp000134" "rp000135" "rp000190b" "rp000190" "rp000190a") # climate name
clims=("rp000134") # climate name

mean=""
yr_span="0020_0039"

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
        cdo -ymonmean ${filename}.nc ${filename}.ymonmean-20.nc 

    done # varnames
done # clims
