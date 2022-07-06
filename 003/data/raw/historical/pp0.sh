#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process

# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
declare -a models=("CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# declare -a models=("bcc-csm1-1") # extended RCP runs
sim="historical"
freq="Amon"
ens="r1i1p1"
yr_span="186001-200512"
# clim=".ymonmean-30"
clim=""

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}${clim}

    cd ${cwd}/${model}

    if ls $cwd/${model}/pp0_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "pp0 was already created. Skipping..."
    else
        full_psic=${cwd}/${model}/psic_${common}.nc
        full_pp0=${cwd}/${model}/pp0_${common}.nc

        python ${cwd}/make_pp0.py ${full_psic} ${full_pp0}
    fi

done
