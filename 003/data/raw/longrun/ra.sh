#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

declare -a models=("MPIESM12_abrupt32x") # climate name
freq="mon"
yr_span="103"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/r1_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else

        full_rsdt=${cwd}/${model}/rsdt_${common}.nc
        full_rsut=${cwd}/${model}/rsut_${common}.nc
        full_rlut=${cwd}/${model}/rlut_${common}.nc
        full_rsds=${cwd}/${model}/rsds_${common}.nc
        full_rsus=${cwd}/${model}/rsus_${common}.nc
        full_rlus=${cwd}/${model}/rlus_${common}.nc
        full_rlds=${cwd}/${model}/rlds_${common}.nc
        full_ra=${cwd}/${model}/ra_${common}.nc

        python ${cwd}/make_ra.py ${full_rsdt} ${full_rsut} ${full_rlut} ${full_rsds} ${full_rsus} ${full_rlds} ${full_rlus} ${full_ra}
    fi

done
