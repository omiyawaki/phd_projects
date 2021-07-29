#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
models=("MPI-ESM-LR/")
sim="rcp85"
freq="Amon"
ens="r1i1p1"
yr_span="200601-200612"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/vmmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_vmmmc=${cwd}/${model}/vmmmc_${common}.nc
        full_vmse=${cwd}/${model}/vmse_${common}.nc
        full_aht=${cwd}/${model}/aht_${common}.nc
        full_vmte=${cwd}/${model}/vmte_${common}.nc

        # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_te.py ${full_vmmmc} ${full_vmse} ${full_aht} ${full_vmte}
        python ${cwd}/make_te.py ${full_vmmmc} ${full_vmse} ${full_aht} ${full_vmte}
    fi

done
