#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

declare -a ens="r1i1p1f1" # ensemble specification 
declare -a models=("MRI-ESM2-0/" "CESM2-WACCM/") # extended RCP runs
# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/" "CESM2-WACCM/") # extended RCP runs
# declare -a ens="r1i1p1f2" # ensemble specification 
# declare -a models=("MIROC-ES2L/") # extended RCP runs
# declare -a ens="r3i1p1f2" # ensemble specification 
# declare -a models=("GISS-E2-1-G/" "GISS-E2-1-H/") # extended RCP runs
# declare -a ens="r4i1p1f2" # ensemble specification 
# declare -a models=("UKESM1-0-LL/") # extended RCP runs

sim="ssp585"
freq="Amon"
yr_begin=2200
yr_end=2299
yr_span="${yr_begin}01-${yr_end}12"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=$(basename ${cwd}/${model}/mse_${freq}_${model}_${sim}_${ens}_*${yr_span}.nc)
    common=${common#mse_}
    common=${common%.nc}


    cd ${cwd}/${model}

    if ls $cwd/${model}/tend_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "tend was already created. Skipping..."
    else
        full_ps=${cwd}/${model}/ps_${common%220001-229912}201501-229912.timmean.nc
        full_mse=${cwd}/${model}/mse_${common}.nc
        full_beta=${cwd}/${model}/beta_${common%220001-229912}201501-229912.timmean.nc
        full_tend=${cwd}/${model}/tend_${common}.nc

        python ${cwd}/make_tend.py ${full_ps} ${full_mse} ${full_beta} ${full_tend}
    fi

done
