#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

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

    if ls $cwd/${model}/qaht_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_evap=${cwd}/${model}/evspsbl_${common}.nc
        full_prec=${cwd}/${model}/pr_${common}.nc
        full_qtend=${cwd}/${model}/qtend_${common}.nc
        full_qaht=${cwd}/${model}/qaht_${common}.nc

        python ${cwd}/make_qaht.py ${full_evap} ${full_prec} ${full_qtend} ${full_qaht} 
    fi

done
