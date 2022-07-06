#!/usr/bin/env bash
# set -euo pipefail

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# models=("HadGEM2-ES/")
models=("CCSM4/")
# declare -a models=("bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
sim="rcp85"
freq="Amon"
ens="r1i1p1"
yr_span="200601-229912"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    # if ls $cwd/${model}/hus_${common} 1> /dev/null 2>&1; then # check if data is already there
    #     echo "hus was already created. Skipping..."
    # else
        full_hus=${cwd}/${model}/hus_${common}.nc
        full_ps=${cwd}/${model}/ps_${common}.nc
        full_husz=${cwd}/${model}/hus_${common}.zonmean.nc
        full_huss=${cwd}/${model}/huss_${common}.nc

        python ${cwd}/make_husz_mask.py ${full_hus} ${full_ps} ${full_husz} ${full_huss}
    # fi

done
