#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
declare -a models=("HadGEM2-ES") # extended RCP runs
sim="historical"
freq="Amon"
ens="r1i1p1"
yr_begin=1860
yr_end=2005
yr_span="${yr_begin}01-${yr_end}12"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/saht_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_aht=${cwd}/${model}/aht_${common}.nc
        full_vmmmc=${cwd}/${model}/vmmmc_${common}.nc
        full_vmse=${cwd}/${model}/vmse_${common}.nc
        full_vmte=${cwd}/${model}/vmte_${common}.nc

        full_qaht=${cwd}/${model}/qaht_${common}.nc
        full_vqmmc=${cwd}/${model}/vqmmc_${common}.nc
        full_vqse=${cwd}/${model}/vqse_${common}.nc
        full_vqte=${cwd}/${model}/vqte_${common}.nc

        full_saht=${cwd}/${model}/saht_${common}.nc
        full_vsmmc=${cwd}/${model}/vsmmc_${common}.nc
        full_vsse=${cwd}/${model}/vsse_${common}.nc
        full_vste=${cwd}/${model}/vste_${common}.nc

        python ${cwd}/make_dry.py ${full_aht} ${full_vmmmc} ${full_vmse} ${full_vmte} ${full_qaht} ${full_vqmmc} ${full_vqse} ${full_vqte} ${full_saht} ${full_vsmmc} ${full_vsse} ${full_vste}
    fi

done
