#!/usr/bin/env bash
# set -euo pipefail

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# models=("MPI-ESM-LR/")
models=("IPSL-CM5A-LR/")
sim="historical"
freq="Amon"
ens="r1i1p1"
yr_span="185001-200512"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/mse_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "MSE was already created. Skipping..."
    else
        full_ta=${cwd}/${model}/ta_${common}.nc
        full_zg=${cwd}/${model}/zg_${common}.nc
        full_hus=${cwd}/${model}/hus_${common}.nc
        full_mse=${cwd}/${model}/mse_${common}.nc

        # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_mse.py ${full_ta} ${full_zg} ${full_hus} ${full_mse}
        python ${cwd}/make_mse.py ${full_ta} ${full_zg} ${full_hus} ${full_mse}
    fi

done
