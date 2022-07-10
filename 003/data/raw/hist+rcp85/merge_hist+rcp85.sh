#!/usr/bin/env bash
set -euo pipefail

table="Amon"
# declare -a vars_gcm=("pr" "prc" "evspsbl" "ftoa" "fsfc" "rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "tend" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "tas" "ts" "r1" "stgadv" "adv" "ra" "stf")
# declare -a vars_gcm=("zg" "ta" "hur" "ps" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "prc" "evspsbl") # list of GCM variables that we want to process
# declare -a vars_gcm=("rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss") # list of GCM variables that we want to process
declare -a vars_gcm=("tend") # list of GCM variables that we want to process

# table="OImon"
# declare -a vars_gcm=("sic") # list of GCM variables that we want to process

# mean=".zonmean.djfmean"
mean=".djfmean"

declare -a ens="r1i1p1" # ensemble specification 
# declare -a ens="r0i0p0" # ensemble specification 
# declare -a models=$(cd /project2/tas1/ockham/data9/tas/CMIP5_RAW && ls -d */) # list of GCM models to process
declare -a models=("HadGEM2-ES/" "bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# declare -a models=("bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# declare -a models=("bcc-csm1-1/") # extended RCP runs

# out_yr_begin=1850
# out_mn_begin=01
# out_yr_end=2005
# out_mn_end=12

cwd=$(pwd) # save current working directory

# for dirs in */; do # loop through all the models
for dirs in ${models[@]}; do # loop through models

    echo $dirs
    mkdir -p $cwd/$dirs # make model directory in processed data folder if it doesn't exist yet
    # rm -f $cwd/$dirs*12.nc # remove nc files that are not ymonmean.nc
    # rm -f $cwd/$dirs*1231.nc
    # rm -f $cwd/$dirs*30.nc 

    echo $(pwd)
    for vars in ${vars_gcm[@]}; do
        echo $vars

        common=${dirs}${vars}_${table}_${dirs%/}
        histnc=$cwd/../historical/${common}_historical_${ens}_186001-200512${mean}.nc 
        rcp85nc=$cwd/../rcp85/${common}_rcp85_${ens}_200601-229912${mean}.nc 
        mergenc=$cwd/../hist+rcp85/${common}_hist+rcp85_${ens}_186001-229912${mean}.nc 
        cdo -O mergetime $histnc $rcp85nc $mergenc

    done
        
done

cd $cwd
