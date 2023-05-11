#!/usr/bin/env bash
set -euo pipefail

table="Amon"
# declare -a vars_gcm=("zg" "evspsbl" "prfrac" "prc" "pr" "tend" "stgadv" "r1" "ps" "ta" "hur" "hus" "tendv" "rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "adv" "ftoa" "ftoacs" "fsfc" "tas" "ts" "ra" "racs") # list of GCM variables that we want to process
# declare -a vars_gcm=("evspsbl" "prfrac" "prc" "pr" "tend" "stgadv" "r1" "ps" "ta" "hur" "hus" "tendv" "adv" "ftoa" "ftoacs" "fsfc" "tas" "ts") # list of GCM variables that we want to process
# declare -a vars_gcm=("evspsbl")
# declare -a vars_gcm=("rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "ra" "racs") # list of GCM variables that we want to process
declare -a vars_gcm=("zg" "ta" "hur" "hus" "mse" "tendv") # list of GCM variables that we want to process
# declare -a vars_gcm=("rlutcs" "rsutcs" "rldscs" "rsdscs" "rsuscs") # list of GCM variables that we want to process
# declare -a vars_gcm=("zg") # list of GCM variables that we want to process

# table="SImon"
# declare -a vars_gcm=("siconc") # list of GCM variables that we want to process

declare -a freq="mon" # data output frequency (e.g. fx for fixed, mon for monthly, day for daily)
mean=""
# mean=".zonmean"
declare -a ens="r1i1p1f1" # ensemble specification 
# declare -a ens="r0i0p0" # ensemble specification 
declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# declare -a models=("IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# declare -a models=("CanESM5/") # extended RCP runs

# declare -a ens="r1i1p1f2" # ensemble specification 
# declare -a models=("MIROC-ES2L/") # extended RCP runs
# declare -a ens="r3i1p1f2" # ensemble specification 
# declare -a models=("GISS-E2-1-G/" "GISS-E2-1-H/") # extended RCP runs
# declare -a ens="r4i1p1f2" # ensemble specification 
# declare -a models=("UKESM1-0-LL/") # extended RCP runs

# out_yr_begin=1850
# out_mn_begin=01
# out_yr_end=2014
# out_mn_end=12

cwd=$(pwd) # save current working directory

# for dirs in */; do # loop through all the models
for dirs in ${models[@]}; do # loop through models

    model=${dirs%/}

    echo $dirs
    mkdir -p $cwd/$dirs # make model directory in processed data folder if it doesn't exist yet
    # rm -f $cwd/$dirs*12.nc # remove nc files that are not ymonmean.nc
    # rm -f $cwd/$dirs*1231.nc
    # rm -f $cwd/$dirs*30.nc 

    echo $(pwd)
    for vars in ${vars_gcm[@]}; do
        echo $vars

        common=${dirs}${vars}_${table}_${model}

        grid=$(basename $cwd/../historical/${model}/rlut_${table}_${model}_historical_${ens}_*_186001-201412.nc)
        echo $grid
        grid=${grid#rlut_${table}_${model}_historical_${ens}_}
        grid=${grid%_186001-201412.nc}

        histnc=$cwd/../historical/${common}_historical_${ens}_${grid}_186001-201412${mean}.nc 
        ssp585nc=$cwd/../ssp585/${common}_ssp585_${ens}_${grid}_201501-230012${mean}.nc 
        mergenc=$cwd/../hist+ssp585/${common}_hist+ssp585_${ens}_${grid}_186001-230012${mean}.nc 
        cdo -O mergetime $histnc $ssp585nc $mergenc

    done
        
done

cd $cwd
