#!/bin/sh

# varnames=("ftoacs" "swtau" "swref" "swabs" "alb" "tradsu" "aclcov" "sni" "xlvi" "xivi" "aclcov" "prfrac" "ftoa" "twtend" "titend" "ttend" "stend" "sn" "ahfres" "tsi" "tsw" "ahfcon" "trafs" "sraf0" "traf0" "srafs" "srad0d" "srad0u" "fsfc" "tsurf" "sradsu" "albedo" "trads" "srads" "srad0" "trad0" "ahfl" "ahfs" "siced" "seaice" "temp2" "precip" "aprl" "aprc" "aps" "tend" "ra" "stf" "ra" "stgadv" "adv" "r1") # list of GCM variables that we want to process
# varnames=("trafs" "srafs" "sraf0" "traf0") # list of GCM variables that we want to process
# varnames=("t" "rhumidity" "q" "tendv")
# varnames=("swref" "swabs" "alb")
varnames=("t")

# clims=("rp000190a" "rp000190b" "rp000190d" "rp000130" "rp000131" "rp000134" "rp000135") # climate name
# clims=("rp000134" "rp000190f") # climate name
clims=("rp000188" "rp000191f") # climate name

mean=""
# yr_span="0020_0039"
yr_span="0040_0252"

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
        cdo -ymonmean ${filename}.shifttime.nc ${filename}.ymonmean-20.nc 
        rm ${filename}.shifttime.nc

    done # varnames
done # clims
