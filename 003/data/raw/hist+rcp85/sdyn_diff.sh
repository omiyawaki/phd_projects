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

    if ls $cwd/${model}/vsmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_saht=${cwd}/${model}/saht_${common}.nc
        full_vsmmc=${cwd}/${model}/vsmmc_${common}.nc
        full_vsse=${cwd}/${model}/vsse_${common}.nc
        full_vste=${cwd}/${model}/vste_${common}.nc
        ref_saht=${rfd}/${model}/saht_${ref_common}.nc
        ref_vsmmc=${rfd}/${model}/vsmmc_${ref_common}.nc
        ref_vsse=${rfd}/${model}/vsse_${ref_common}.nc
        ref_vste=${rfd}/${model}/vste_${ref_common}.nc
        diff_saht=${cwd}/${model}/dsaht_${common}.nc
        diff_vsmmc=${cwd}/${model}/dvsmmc_${common}.nc
        diff_vsse=${cwd}/${model}/dvsse_${common}.nc
        diff_vste=${cwd}/${model}/dvste_${common}.nc

        python ${cwd}/diff_sdyn.py ${full_saht} ${full_vsmmc} ${full_vsse} ${full_vste} ${ref_saht} ${ref_vsmmc} ${ref_vsse} ${ref_vste} ${diff_saht} ${diff_vsmmc} ${diff_vsse} ${diff_vste}


    fi

done
