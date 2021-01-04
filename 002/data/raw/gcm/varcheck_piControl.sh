#!/usr/bin/env bash
set -euo pipefail

# declare -a vars_gcm=("va" "zg" "wap" "ta" "hur" "ps" "hurs" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "rldscs" "rlutcs" "rsdscs" "rsuscs" "rsutcs" "hfls" "hfss" "pr" "prc" "evspsbl" "vas") # list of GCM variables that we want to process
declare -a vars_gcm=("zg" "ta" "hur" "ps" "hurs" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "prc" "evspsbl") # list of GCM variables that we want to process
# declare -a models=("MPI-ESM-LR/" "IPSL-CM5A-LR/") # list of GCM models to process
declare -a models=$(cd /project2/tas1/CMIP5_piControl/ && ls -d */) # list of GCM models to process

cwd=$(pwd) # save current working directory
cd ../piControl_raw # switch to directory with raw data
rwd=$(pwd) # save raw data directory

# for dirs in */; do # loop through all the models
for dirs in ${models[@]}; do # loop through models
    cd $rwd/$dirs # go in the model directory
    for vars in ${vars_gcm[@]}; do
        if ls $cwd/${dirs}${vars}_*piControl*.nc 1> /dev/null 2>&1; then # check if data is already there
            :
        else
            echo "${dirs}: ${vars} is missing."
        fi
    done
done

cd $cwd
