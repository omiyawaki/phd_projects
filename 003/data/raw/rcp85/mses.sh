#!/usr/bin/env bash
# set -euo pipefail

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

    if ls $cwd/${model}/mses_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mses was already created. Skipping..."
    else
        full_ts=${cwd}/${model}/ts_${common}.nc
        full_orog=/project2/tas1/ockham/data9/tas/CMIP5_RAW/MPI-ESM-LR/rcp85/atmos/fx/orog/r0i0p0/orog_fx_${model}_${sim}_r0i0p0.nc
        full_huss=${cwd}/${model}/huss_${common}.nc
        full_mses=${cwd}/${model}/mses_${common}.nc

        # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_mses.py ${full_ts} ${full_orog} ${full_huss} ${full_mses}
        python ${cwd}/make_mses.py ${full_ts} ${full_orog} ${full_huss} ${full_mses}
    fi

done
