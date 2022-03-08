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

    full_vmse=${cwd}/${model}/steady/vmse_*.nc
    full_merge=${cwd}/${model}/vmse_${common}.nc

    cdo -O mergetime $full_vmse $full_merge

    full_vmmmc=${cwd}/${model}/steady/vmmmc_*.nc
    full_merge=${cwd}/${model}/vmmmc_${common}.nc

    cdo -O mergetime $full_vmmmc $full_merge

done
