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

    if ls $cwd/${model}/vmmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_aht=${cwd}/${model}/aht_${common}.nc
        full_vmmmc=${cwd}/${model}/vmmmc_${common}.nc
        full_vmse=${cwd}/${model}/vmse_${common}.nc
        full_vmte=${cwd}/${model}/vmte_${common}.nc
        ref_aht=${rfd}/${model}/aht_${ref_common}.nc
        ref_vmmmc=${rfd}/${model}/vmmmc_${ref_common}.nc
        ref_vmse=${rfd}/${model}/vmse_${ref_common}.nc
        ref_vmte=${rfd}/${model}/vmte_${ref_common}.nc
        diff_aht=${cwd}/${model}/daht_${common}.nc
        diff_vmmmc=${cwd}/${model}/dvmmmc_${common}.nc
        diff_vmse=${cwd}/${model}/dvmse_${common}.nc
        diff_vmte=${cwd}/${model}/dvmte_${common}.nc

        python ${cwd}/diff_dyn.py ${full_aht} ${full_vmmmc} ${full_vmse} ${full_vmte} ${ref_aht} ${ref_vmmmc} ${ref_vmse} ${ref_vmte} ${diff_aht} ${diff_vmmmc} ${diff_vmse} ${diff_vmte}


    fi

done
