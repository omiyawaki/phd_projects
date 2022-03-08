#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# declare -a models=("bcc-csm1-1") # extended RCP runs
sim="hist+rcp85"
clim=".djfmean"
freq="Amon"
ens="r1i1p1"
yr_begin=1860
yr_end=2299
yr_span="${yr_begin}01-${yr_end}12"

ref_sim="historical"
ref_clim=".ymonmean-30"
ref_yr_begin=1860
ref_yr_end=2005
ref_yr_span="${ref_yr_begin}01-${ref_yr_end}12"

cwd=$(pwd) # save current working directory
rfd=${cwd}/../${ref_sim} # reference data directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}${clim}
    ref_common=${freq}_${model}_${ref_sim}_${ens}_${ref_yr_span}${ref_clim}${clim}

    cd ${cwd}/${model}

    if ls $cwd/${model}/vqmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_qaht=${cwd}/${model}/qaht_${common}.nc
        full_vqmmc=${cwd}/${model}/vqmmc_${common}.nc
        full_vqse=${cwd}/${model}/vqse_${common}.nc
        full_vqte=${cwd}/${model}/vqte_${common}.nc
        ref_qaht=${rfd}/${model}/qaht_${ref_common}.nc
        ref_vqmmc=${rfd}/${model}/vqmmc_${ref_common}.nc
        ref_vqse=${rfd}/${model}/vqse_${ref_common}.nc
        ref_vqte=${rfd}/${model}/vqte_${ref_common}.nc
        diff_qaht=${cwd}/${model}/dqaht_${common}.nc
        diff_vqmmc=${cwd}/${model}/dvqmmc_${common}.nc
        diff_vqse=${cwd}/${model}/dvqse_${common}.nc
        diff_vqte=${cwd}/${model}/dvqte_${common}.nc

        python ${cwd}/diff_qdyn.py ${full_qaht} ${full_vqmmc} ${full_vqse} ${full_vqte} ${ref_qaht} ${ref_vqmmc} ${ref_vqse} ${ref_vqte} ${diff_qaht} ${diff_vqmmc} ${diff_vqse} ${diff_vqte}


    fi

done
