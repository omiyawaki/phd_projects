#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# models=("IPSL-CM5A-LR/")
models=("MPI-ESM-LR/")
sim="historical"
freq="Amon"
ens="r1i1p1"
yr_span="195001-200512"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    full_eke=${cwd}/${model}/ke/eke_*.nc
    full_merge=${cwd}/${model}/eke_${common}.nc

    cdo mergetime $full_eke $full_merge

    full_veke=${cwd}/${model}/ke/veke_*.nc
    full_merge=${cwd}/${model}/veke_${common}.nc

    cdo mergetime $full_veke $full_merge

done
