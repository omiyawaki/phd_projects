#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extmaxed RCP runs
# declare -a models=("HadGEM2-ES") # extmaxed RCP runs
# ens="r1i1p1"

sim="historical"
freq="Amon"
yr_span="186001-200512"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/tmax_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "tmax was already created. Skipping..."
    else
        full_ps=${cwd}/${model}/ps_${common}.nc
        full_ta=${cwd}/${model}/ta_${common}.nc
        full_tmax=${cwd}/${model}/tmax_${common}.nc

        python ${cwd}/make_tmax.py ${full_ps} ${full_ta} ${full_tmax}
    fi

done
