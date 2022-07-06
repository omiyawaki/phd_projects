#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
declare -a models=("HadGEM2-ES/" "bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# models=("MPI-ESM-LR/")
# models=("IPSL-CM5A-LR/")
sim="hist+rcp85"
freq="Amon"
ens="r1i1p1"
yr_span="186001-229912"
lev="92500"
mean=".djfmean"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    zonmean=${freq}_${model}_${sim}_${ens}_${yr_span}.zonmean${mean}
    common=${freq}_${model}_${sim}_${ens}_${yr_span}${mean}

    cd ${cwd}/${model}

    if ls $cwd/${model}/diffv${lev}_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "diffv was already created. Skipping..."
    else
        full_vmte=${cwd}/${model}/vmte_sm_${common}.nc
        full_gmse=${cwd}/${model}/gmse${lev}_sm_${zonmean}.nc
        full_diffv=${cwd}/${model}/diffv${lev}_sm_${common}.nc

        python ${cwd}/make_diffv.py ${full_vmte} ${full_gmse} ${full_diffv}
    fi

done
