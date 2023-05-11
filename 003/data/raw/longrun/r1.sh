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
        full_ra=${cwd}/${model}/ra_${common}.nc
        full_stf=${cwd}/${model}/stf_${common}.nc
        full_r1=${cwd}/${model}/r1_${common}.nc

        python ${cwd}/make_r1.py ${full_ra} ${full_stf} ${full_r1}
    fi

done
