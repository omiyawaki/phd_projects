#!/usr/bin/env bash
set -euo pipefail

# declare -a vars_gcm=("zg" "ta" "hur" "ps" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "prc" "evspsbl") # list of GCM variables that we want to process
# declare -a vars_gcm=("rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss") # list of GCM variables that we want to process
# declare -a vars_gcm=("rlutcs" "rsutcs" "rldscs" "rsdscs" "rsuscs") # list of GCM variables that we want to process
# declare -a vars_gcm=("ps" "tas" "ta" "zg" "hus") # list of GCM variables that we want to process
# declare -a vars_gcm=("aht" "vmmmc" "vmse" "vmte") # list of GCM variables that we want to process
declare -a vars_gcm=("qaht" "saht" "vqmmc" "vsmmc" "vqse" "vsse" "vqte" "vste") # list of GCM variables that we want to process
# declare -a vars_gcm=("ta" "hus" "hur") # list of GCM variables that we want to process
# declare -a vars_gcm=("ra" "stf" "stgadv" "adv" "r1") # list of GCM variables that we want to process
declare -a realm=("atmos")
declare -a freq="mon" # data output frequency (e.g. fx for fixed, mon for monthly, day for daily)
table="Amon"
declare -a ens="r1i1p1" # ensemble specification 
# declare -a ens="r0i0p0" # ensemble specification 
# declare -a models=$(cd /project2/tas1/ockham/data9/tas/CMIP5_RAW && ls -d */) # list of GCM models to process
# declare -a models=("HadGEM2-ES/" "bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "GISS-E2-H/" "GISS-E2-R/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# declare -a models=("bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
declare -a models=("CCSM4/") # extended RCP runs

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
        histnc=$cwd/../historical/${common}_historical_${ens}_186001-200512.nc 
        rcp85nc=$cwd/../rcp85/${common}_rcp85_${ens}_200601-229912.nc 
        mergenc=$cwd/../hist+rcp85/${common}_hist+rcp85_${ens}_186001-229912.nc 
        cdo -O mergetime $histnc $rcp85nc $mergenc

    done
        
done

cd $cwd
