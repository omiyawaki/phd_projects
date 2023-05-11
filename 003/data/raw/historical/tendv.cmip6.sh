#!/usr/bin/env bash
# set -euo pipefail

# load python environment
# module unload python
# module load python

# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# declare -a models=("CESM2-WACCM/") # extended RCP runs
# ens="r1i1p1f1"
# declare -a ens="r1i1p1f2" # ensemble specification 
# declare -a models=("MIROC-ES2L/") # extended RCP runs
declare -a ens="r3i1p1f2" # ensemble specification 
declare -a models=("GISS-E2-1-G/" "GISS-E2-1-H/") # extended RCP runs
# declare -a ens="r4i1p1f2" # ensemble specification 
# declare -a models=("UKESM1-0-LL/") # extended RCP runs
yr_span="186001-201412"

sim="historical"
freq="Amon"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    grid=$(basename ${cwd}/${model}/rlut_${freq}_${model}_${sim}_${ens}*${yr_span}.nc)
    grid=${grid#rlut_${freq}_${model}_${sim}_${ens}_}
    grid=${grid%_${yr_span}.nc}

    common=${freq}_${model}_${sim}_${ens}_${grid}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/tendv_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "tendv was already created. Skipping..."
    else
        full_ps=${cwd}/${model}/ps_${common}.timmean.nc
        full_mse=${cwd}/${model}/mse_${common}.nc
        full_beta=${cwd}/${model}/beta_${common}.timmean.nc
        full_tendv=${cwd}/${model}/tendv_${common}.nc

        python ${cwd}/make_tendv.py ${full_ps} ${full_mse} ${full_beta} ${full_tendv}
    fi

done
