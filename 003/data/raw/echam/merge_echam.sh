#!/usr/bin/env bash
set -euo pipefail

# declare -a vars_gcm=("trads" "srads" "srad0" "trad0" "ahfl" "ahfs") # list of GCM variables that we want to process
declare -a vars_gcm=("friac") # list of GCM variables that we want to process
declare -a clims=("rp000140") # climate name
declare -a skip_files=("200601-209912.nc" "_eady.nc")

cwd=$(pwd) # save current working directory
cd /project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1 # switch to directory with data
rwd=$(pwd) # save raw data directory

# for dirs in */; do # loop through all the models
for dirs in ${clims[@]}; do # loop through models
    echo $dirs
    mkdir -p $cwd/$dirs # make model directory in processed data folder if it doesn't exist yet
    cd $rwd/$dirs # go in the model directory

    for vars in ${vars_gcm[@]}; do
        echo $vars
        # if ls $cwd/${dirs}${vars}_*historical*${out_yr_begin}${out_mn_begin}-${out_yr_end}${out_mn_end}*.nc 1> /dev/null 2>&1; then # check if data is already there
        #     echo "${vars} was already converted. Skipping..."
        # else
            pattern=$(ls BOT_${dirs}_{0001..0039}.nc)
            files=( $pattern )

            if [[ ! -f "$cwd/$dirs/BOT_${dirs}_0001_0039.nc" ]]; then # check if merged file already exists
                cdo -O mergetime ${files[@]} $cwd/$dirs/BOT_${dirs}_0001_0039.nc # combine multiple files into one
            fi

            cdo -O selvar,$vars $cwd/$dirs/BOT_${dirs}_0001_0039.nc $cwd/$dirs/${vars}_${dirs}_0001_0039.nc 
            cdo -O yearmean $cwd/$dirs/${vars}_${dirs}_0001_0039.nc $cwd/$dirs/${vars}_${dirs}_0001_0039.yearmean.nc  # combine multiple files into one
        # fi
        cd $rwd/$dirs # go in the model directory
    done
    
done

cd $cwd
