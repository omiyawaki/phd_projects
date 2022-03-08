#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# declare -a models=("bcc-csm1-1") # extended RCP runs
sim="rcp85"
freq="Amon"
ens="r1i1p1"
yr_begin=2006
yr_end=2299
yr_span="${yr_begin}01-${yr_end}12"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/vqmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_vqmmc=${cwd}/${model}/vqmmc_${common}.nc
        full_vqse=${cwd}/${model}/vqse_${common}.nc
        full_qaht=${cwd}/${model}/qaht_${common}.nc
        full_vqte=${cwd}/${model}/vqte_${common}.nc

        python ${cwd}/infer_qte.py ${full_vqmmc} ${full_vqse} ${full_qaht} ${full_vqte}
    fi

done
