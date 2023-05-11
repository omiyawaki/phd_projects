#!/usr/bin/env bash
set -euo pipefail

table="fx"
declare -a vars_gcm=("orog") # list of GCM variables that we want to process

declare -a freq="fx" # data output frequency (e.g. fx for fixed, mon for monthly, day for daily)
declare -a ens="r1i1p1f1" # ensemble specification 
# declare -a ens="r0i0p0" # ensemble specification 
declare -a models=("CESM2-WACCM/") # extended RCP runs
# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# declare -a models=("ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs

# declare -a ens="r1i1p1f2" # ensemble specification 
# declare -a models=("MIROC-ES2L/") # extended RCP runs
# declare -a ens="r3i1p1f2" # ensemble specification 
# declare -a models=("GISS-E2-1-G/" "GISS-E2-1-H/") # extended RCP runs
# declare -a ens="r4i1p1f2" # ensemble specification 
# declare -a models=("UKESM1-0-LL/") # extended RCP runs


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

        grid=$(basename $cwd/../historical/${model}/orog_${table}_${model}_historical_${ens}_*.nc)
        echo $grid
        grid=${grid#orog_${table}_${model}_historical_${ens}_}
        grid=${grid%.nc}

        histnc=$cwd/../historical/${common}_historical_${ens}_${grid}.nc 
        ssp585nc=$cwd/../ssp585/${common}_ssp585_${ens}_${grid}.nc 
        mergenc=$cwd/../hist+ssp585/${common}_hist+ssp585_${ens}_${grid}.nc 
        # cp $ssp585nc $mergenc
        cp $histnc $mergenc

    done
        
done

cd $cwd
