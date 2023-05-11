#!/usr/bin/env bash
set -euo pipefail

# dim="BOT_dm"
# # declare -a vars_gcm=("siced" "sni") # list of GCM variables that we want to process
# declare -a vars_gcm=("tsurf" "tsi") # list of GCM variables that we want to process

dim="ATM_dm_ml"
# declare -a vars_gcm=("t" "q" "geopoth") # list of GCM variables that we want to process
declare -a vars_gcm=("t") # list of GCM variables that we want to process

# declare -a clims=("rp000190" "rp000190a" "rp000190b" "rp000190c") # climate name
# yr0=0001
# yr1=0009

# declare -a clims=("rp000190a") # climate name
# yr0=0001
# yr1=0009

# declare -a clims=("rp000190") # climate name
# yr0=0001
# yr1=0039

# declare -a clims=("rp000191") # climate name
# yr0=0040
# yr1=0129

# declare -a clims=("rp000130") # climate name
# yr0=0001
# yr1=0079

# declare -a clims=("rp000134" "rp000135" "rp000190" "rp000190a" "rp000190b") # climate name
declare -a clims=("rp000195") # climate name
# declare -a clims=("rp000130" "rp000131" "rp000134" "rp000135") # climate name
yr0=0020
yr1=0039

# declare -a clims=("rp000191") # climate name
# # declare -a clims=("rp000188" "rp000191b") # climate name
# yr0=0040
# yr1=0339

# declare -a clims=("rp000191b") # climate name
# # declare -a clims=("rp000184" "rp000185" "rp000188" "rp000189") # climate name
# yr0=0040
# yr1=0139

declare -a skip_files=("200601-209912.nc" "_eady.nc")

###################################
# MAKE SURE TO CHANGE YEAR RANGE
# MANUALLY BELOW AS WELL
###################################


cwd=$(pwd) # save current working directory
cd /project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1 # switch to directory with data
rwd=$(pwd) # save raw data directory

for dirs in ${clims[@]}; do # loop through models

    echo $dirs
    mkdir -p $cwd/$dirs # make model directory in processed data folder if it doesn't exist yet
    cd $rwd/$dirs # go in the model directory

    for vars in ${vars_gcm[@]}; do
        echo $vars
        # if ls $cwd/${dirs}${vars}_*historical*${out_yr_begin}${out_mn_begin}-${out_yr_end}${out_mn_end}*.nc 1> /dev/null 2>&1; then # check if data is already there
        #     echo "${vars} was already converted. Skipping..."
        # else

        ##########################################################
        # !!! UPDATE THE YEAR RANGE BELOW MANUALLY
        ##########################################################

            pattern=$(ls ${dim}_${dirs}_{0020..0039}.nc)
            # pattern=$(ls ${dim}_${dirs}_{0001..0009}.nc)
            # pattern=$(ls ${dim}_${dirs}_{0001..0039}.nc)
            # pattern=$(ls ${dim}_${dirs}_{0040..0100}.nc)
            # pattern=$(ls ${dim}_${dirs}_{0040..0239}.nc)
            # pattern=$(ls ${dim}_${dirs}_{0040..0339}.nc)
            # pattern=$(ls ${dim}_${dirs}_{0001..0079}.nc)
            files=( $pattern )
            echo $files

            if [[ ! -f "$cwd/$dirs/${dim}_${dirs}_${yr0}_${yr1}.nc" ]]; then # check if merged file already exists
                cdo -O mergetime ${files[@]} $cwd/$dirs/${dim}_${dirs}_${yr0}_${yr1}.nc # combine multiple files into one
            fi

            cdo -O selvar,$vars $cwd/$dirs/${dim}_${dirs}_${yr0}_${yr1}.nc $cwd/$dirs/${vars}_${dirs}_${yr0}_${yr1}.nc 
        # fi
        cd $rwd/$dirs # go in the model directory
    done
    
done

cd $cwd
