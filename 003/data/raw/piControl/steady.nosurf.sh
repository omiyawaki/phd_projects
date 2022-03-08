#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
models=("CCSM4/")
sim="piControl"
freq="Amon"
ens="r1i1p1"
yr_span="125001-129912"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/vmmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_ps=${cwd}/${model}/ps_${common}.nc
        full_va=${cwd}/${model}/va_${common}.nc
        full_mse=${cwd}/${model}/mse_${common}.nc
        # full_divmmc=${cwd}/${model}/divmmc_${common}.nc
        # full_divse=${cwd}/${model}/divse_${common}.nc
        full_vmmmc=${cwd}/${model}/vmmmc_${common}.nc
        full_vmse=${cwd}/${model}/vmse_${common}.nc

        python ${cwd}/make_steady.nosurf.py ${full_ps} ${full_va} ${full_mse} ${full_vmmmc} ${full_vmse}
    fi

done
