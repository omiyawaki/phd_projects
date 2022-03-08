#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
sim="historical"
freq="Amon"
ens="r1i1p1"
yr_span="186001-200512"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/qtend_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "qtend was already created. Skipping..."
    else
        full_ps=${cwd}/${model}/ps_${common}.timmean.nc
        full_hus=${cwd}/${model}/hus_${common}.nc
        full_beta=${cwd}/${model}/beta_${common}.timmean.nc
        full_qtend=${cwd}/${model}/qtend_${common}.nc

        python ${cwd}/make_qtend.py ${full_ps} ${full_hus} ${full_beta} ${full_qtend}
    fi

done
