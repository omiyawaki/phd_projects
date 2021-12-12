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
yr_span="200601-230012"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    full_vmte=${cwd}/${model}/te/vmte_*.nc
    full_merge=${cwd}/${model}/vmte_${common}.nc

    cdo mergetime $full_vmte $full_merge

done
