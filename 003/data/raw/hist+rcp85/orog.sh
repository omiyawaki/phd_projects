#!/usr/bin/env bash
set -euo pipefail

table="fx"
declare -a vars_gcm=("orog") # list of GCM variables that we want to process

declare -a ens="r0i0p0" # ensemble specification 
# declare -a models=$(cd /project2/tas1/ockham/data9/tas/CMIP5_RAW && ls -d */) # list of GCM models to process
# declare -a models=("HadGEM2-ES/" "bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# declare -a models=("bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
declare -a models=("GISS-E2-H/" "GISS-E2-R/") # extended RCP runs

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
        histnc=$cwd/../historical/${common}_historical_${ens}.nc 
        rcp85nc=$cwd/../rcp85/${common}_rcp85_${ens}.nc 
        mergenc=$cwd/../hist+rcp85/${common}_hist+rcp85_${ens}.nc 

        # cp $rcp85nc $mergenc
        cp $histnc $mergenc

    done
        
done

cd $cwd
