#!/usr/bin/env bash
# set -euo pipefail

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
models=("MPI-ESM-LR/")
sim="historical"
freq="Amon"
ens="r1i1p1"
yr_span="186001-200512"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/tasi_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "tasi was already created. Skipping..."
    else
        full_ps=${cwd}/${model}/ps_${common}.nc
        full_tas=${cwd}/${model}/tas_${common}.nc
        full_ta=${cwd}/${model}/ta_${common}.nc
        full_tasi=${cwd}/${model}/tasi_${common}.nc

        # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_sigma.py ${full_ps} ${full_tas} ${full_ta} ${full_tasi}
        python ${cwd}/make_sigma.py ${full_ps} ${full_tas} ${full_ta} ${full_tasi}
    fi

done
