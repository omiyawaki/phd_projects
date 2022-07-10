#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

declare -a ens="r1i1p1f1" # ensemble specification 
declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs

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

    if ls $cwd/${model}/r1_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else

        full_hfls=${cwd}/${model}/hfls_${common}.nc
        full_hfss=${cwd}/${model}/hfss_${common}.nc
        full_stf=${cwd}/${model}/stf_${common}.nc

        python ${cwd}/make_stf.py ${full_hfls} ${full_hfss} ${full_stf}
    fi

done
