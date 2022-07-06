#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# declare -a models=("CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
declare -a models=("CNRM-CM5") # extended RCP runs
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

    if ls $cwd/${model}/r1_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else

        full_rlutcs=${cwd}/${model}/rlutcs_${common}.nc
        full_rsutcs=${cwd}/${model}/rsutcs_${common}.nc
        full_rsdt=${cwd}/${model}/rsdt_${common}.nc
        full_rlus=${cwd}/${model}/rlus_${common}.nc
        full_rldscs=${cwd}/${model}/rldscs_${common}.nc
        full_rsdscs=${cwd}/${model}/rsdscs_${common}.nc
        full_rsuscs=${cwd}/${model}/rsuscs_${common}.nc
        full_racs=${cwd}/${model}/racs_${common}.nc

        python ${cwd}/make_racs.py ${full_rlutcs} ${full_rsutcs} ${full_rsdt} ${full_rlus} ${full_rldscs} ${full_rsdscs} ${full_rsuscs} ${full_racs}
    fi

done
