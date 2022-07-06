#!/usr/bin/env bash
# set -euo pipefail

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# models=("HadGEM2-ES/")
# models=("bcc-csm1-1/")
declare -a models=("bcc-csm1-1/" "HadGEM2-ES/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
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

    # if ls $cwd/${model}/tendv_${common} 1> /dev/null 2>&1; then # check if data is already there
    #     echo "tendv was already created. Skipping..."
    # else
        full_tendv=${cwd}/${model}/tendv_${common}.nc
        full_tendvz=${cwd}/${model}/tendv_${common}.zonmean.nc

        python ${cwd}/make_tendvz.py ${full_tendv} ${full_tendvz}
    # fi

done
