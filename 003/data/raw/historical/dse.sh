#!/usr/bin/env bash
# set -euo pipefail

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# models=("HadGEM2-ES/")
# models=("IPSL-CM5A-LR/")
# declare -a models=("bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "GISS-E2-H/" "GISS-E2-R/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
declare -a models=("bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
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

    if ls $cwd/${model}/dse_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "dse was already created. Skipping..."
    else
        full_ta=${cwd}/${model}/ta_${common}.nc
        full_zg=${cwd}/${model}/zg_${common}.nc
        full_dse=${cwd}/${model}/dse_${common}.nc

        # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_dse.py ${full_ta} ${full_zg} ${full_hus} ${full_dse}
        python ${cwd}/make_dse.py ${full_ta} ${full_zg} ${full_dse}
    fi

done
