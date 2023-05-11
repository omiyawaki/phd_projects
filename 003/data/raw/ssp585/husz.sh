#!/usr/bin/env bash
# set -euo pipefail

# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "CESM2-WACCM/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# declare -a models=("CESM2-WACCM/") # extended RCP runs
# ens="r1i1p1f1"
# declare -a ens="r1i1p1f2" # ensemble specification 
# declare -a models=("MIROC-ES2L/") # extended RCP runs
declare -a ens="r3i1p1f2" # ensemble specification 
declare -a models=("GISS-E2-1-G/" "GISS-E2-1-H/") # extended RCP runs
# declare -a ens="r4i1p1f2" # ensemble specification 
# declare -a models=("UKESM1-0-LL/") # extended RCP runs
yr_span="201501-229912"

sim="ssp585"
freq="Amon"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    grid=$(basename ${cwd}/${model}/rlut_${freq}_${model}_${sim}_${ens}_*_${yr_span}.nc)
    grid=${grid#rlut_${freq}_${model}_${sim}_${ens}_}
    grid=${grid%_${yr_span}.nc}

    common=${freq}_${model}_${sim}_${ens}_${grid}_${yr_span}

    cd ${cwd}/${model}

    # if ls $cwd/${model}/hus_${common} 1> /dev/null 2>&1; then # check if data is already there
    #     echo "hus was already created. Skipping..."
    # else
        full_hus=${cwd}/${model}/hus_${common}.nc
        full_ps=${cwd}/${model}/ps_${common}.timmean.nc
        full_husz=${cwd}/${model}/hus_${common}.zonmean.nc

        python ${cwd}/make_husz.py ${full_hus} ${full_ps} ${full_husz}
    # fi

done
