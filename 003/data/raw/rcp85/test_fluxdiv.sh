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
yr_span="20060101-20061231"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/divte_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "divte was already created. Skipping..."
    else
        full_vmte=${cwd}/${model}/te/vmte_${common}.nc
        full_divte=${cwd}/${model}/te/divte_${common}.nc

        # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_fluxdiv.py ${full_aht} ${full_vmmmc} ${full_vmse} ${full_vmte} ${full_divaht} ${full_divmmc} ${full_divse} ${full_divte}
        python ${cwd}/make_test_fluxdiv.py ${full_vmte} ${full_divte}
    fi

done
