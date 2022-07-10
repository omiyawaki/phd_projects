#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

declare -a ens="r1i1p1f1" # ensemble specification 
declare -a models=("MRI-ESM2-0/") # extended RCP runs
# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs

sim="ssp585"
freq="Amon"
yr_begin=2015
yr_end=2300
yr_span="${yr_begin}01-${yr_end}12"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=$(basename ${cwd}/${model}/rlut_${freq}_${model}_${sim}_${ens}_*${yr_span}.nc)
    common=${common#rlut_}
    common=${common%.nc}


    cd ${cwd}/${model}

    if ls $cwd/${model}/tend_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "tend was already created. Skipping..."
    else
        full_ps=${cwd}/${model}/ps_${common}.timmean.nc
        full_mse=${cwd}/${model}/mse_${common}.nc
        full_beta=${cwd}/${model}/beta_${common}.timmean.nc
        full_tend=${cwd}/${model}/tend_${common}.nc

        python ${cwd}/make_tend.py ${full_ps} ${full_mse} ${full_beta} ${full_tend}
    fi

done
