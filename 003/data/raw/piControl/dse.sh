#!/usr/bin/env bash
# set -euo pipefail

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# models=("ACCESS1-0/")
# yr_span="075001-079912"
# models=("CCSM4/")
# yr_span="125001-129912"
models=("MPI-ESM-LR/")
yr_span="280001-284912"
sim="piControl"
freq="Amon"
ens="r1i1p1"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/dse_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "dse was already created. Skipping..."
    else
        full_ta=${cwd}/${model}/ta_${common}.nc
        full_zg=${cwd}/${model}/zg_${common}.nc
        full_dse=${cwd}/${model}/dse_${common}.nc

        # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_dse.py ${full_ta} ${full_zg} ${full_hus} ${full_dse}
        python ${cwd}/make_dse.py ${full_ta} ${full_zg} ${full_dse}
    fi

done
