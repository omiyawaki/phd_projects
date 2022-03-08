#!/usr/bin/env bash
# set -euo pipefail

# load python environment
# module unload python
# module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
models=("ACCESS1-0/")
yr_begin="0750"
yr_end="0799"
# models=("CCSM4/")
# yr_begin="1250"
# yr_end="1299"
# models=("MPI-ESM-LR/")
# yr_begin="2800"
# yr_end="2849"
sim="piControl"
freq="Amon"
ens="r1i1p1"
# mean=".timmean"
mean=""
yr_span="${yr_begin}01-${yr_end}12"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/beta_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "beta was already created. Skipping..."
    else
        full_ps=${cwd}/${model}/ps_${common}${mean}.nc
        full_va=${cwd}/${model}/va_${common}.nc
        full_beta=${cwd}/${model}/beta_${common}${mean}.nc

        python ${cwd}/make_beta.py ${full_ps} ${full_va} ${full_beta}

    fi

done
