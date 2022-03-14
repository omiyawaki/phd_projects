#!/usr/bin/env bash
set -euo pipefail

# declare -a vars_gcm=("zg" "ta" "hur" "ps" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "prc" "evspsbl") # list of GCM variables that we want to process
# declare -a vars_gcm=("rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "prc" "sic") # list of GCM variables that we want to process
# declare -a vars_gcm=("ta" "hus" "zg")
# declare -a vars_gcm=("rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs")
# declare -a vars_gcm=("ps" "va" "tas") # list of GCM variables that we want to process
# declare -a vars_gcm=("ps" "tas" "ts" "ta" "hus" "zg") # list of GCM variables that we want to process
declare -a vars_gcm=("evspsbl") # list of GCM variables that we want to process
# declare -a vars_gcm=("sic") # list of GCM variables that we want to process
# declare -a vars_gcm=("clt" "clwvi") # list of GCM variables that we want to process
declare -a clim="rcp85" # climate name
declare -a freq="Amon" # data output frequency (e.g. fx for fixed, mon for monthly, day for daily)
declare -a ens="r1i1p1" # ensemble specification 
declare -a models="HadGEM2-ES/" # list of GCM models to process

in_yr_out=2300

out_yr_begin=2005
out_mn_begin=11
# out_yr_end=2299
out_mn_end=12

model="HadGEM2-ES"

common=${freq}_${model}_historical_${ens}
common_rcp=${freq}_${model}_rcp85_${ens}

cwd=$(pwd) # save current working directory
rcpdir=${cwd}/../rcp85 # rcp dir

cd ${cwd}/${model}

for varname in ${vars_gcm[@]}; do
    if ls $cwd/${model}/mse_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "MSE was already created. Skipping..."
    else
        hist_0511=${cwd}/${model}/${varname}_${common}_186001-200511.nc
        hist_0512=${cwd}/${model}/${varname}_${common}_186001-200512.nc
        rcp_0512=${rcpdir}/${model}/${varname}_${common_rcp}_200511-200512.nc
        # hist_0512=${cwd}/${model}/${varname}_${common}_186001-200512.nc

        ntime_hist=$(cdo -ntime $hist_0512)

        if [ $ntime_hist == 1751 ]; then
            echo "Appending 2005-12 data for ${varname} file..."
            mv $hist_0512 $hist_0511
            cdo mergetime $hist_0511 $rcp_0512 $hist_0512
        fi

    fi
done
