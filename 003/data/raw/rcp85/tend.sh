#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
models=("IPSL-CM5A-LR/")
sim="rcp85"
freq="Amon"
ens="r1i1p1"
yr_span="200601-230012"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/tend_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "tend was already created. Skipping..."
    else
        full_ps=${cwd}/${model}/ps_${common}.nc
        full_mse=${cwd}/${model}/mse_${common}.nc
        full_tend=${cwd}/${model}/tend_${common}.nc

        # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_tend.py ${full_ps} ${full_mse} ${full_tend}
        python ${cwd}/make_tend.py ${full_ps} ${full_mse} ${full_tend}
    fi

done
