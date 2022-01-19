#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
models=("MPI-ESM-LR/")
# models=("IPSL-CM5A-LR/")
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

    # if ls $cwd/${model}/vmmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
    #     echo "mmc and se was already created. Skipping..."
    # else
        full_aht=${cwd}/${model}/aht_${common}.nc
        full_vmte=${cwd}/${model}/vmte_${common}.nc
        full_divaht=${cwd}/${model}/divaht_${common}.nc
        full_divte=${cwd}/${model}/divte_${common}.nc

        # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_fluxdiv.py ${full_aht} ${full_vmmmc} ${full_vmse} ${full_vmte} ${full_divaht} ${full_divmmc} ${full_divse} ${full_divte}
        python ${cwd}/make_fluxdiv.py ${full_aht} ${full_vmte} ${full_divaht} ${full_divte}
    # fi

done
