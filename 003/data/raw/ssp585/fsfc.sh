#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

declare -a ens="r1i1p1f1" # ensemble specification 
declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "CESM2-WACCM/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs

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

    if ls $cwd/${model}/fsfc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else

        full_rlus=${cwd}/${model}/rlus_${common}.nc
        full_rlds=${cwd}/${model}/rlds_${common}.nc
        full_rsds=${cwd}/${model}/rsds_${common}.nc
        full_rsus=${cwd}/${model}/rsus_${common}.nc
        full_hfls=${cwd}/${model}/hfls_${common}.nc
        full_hfss=${cwd}/${model}/hfss_${common}.nc
        full_fsfc=${cwd}/${model}/fsfc_${common}.nc

        python ${cwd}/make_fsfc.py ${full_rlus} ${full_rlds} ${full_rsds} ${full_rsus} ${full_hfls} ${full_hfss} ${full_fsfc}
    fi

done
