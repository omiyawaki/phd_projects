#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# declare -a models=("CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
declare -a models=("GISS-E2-H/" "GISS-E2-R/") # extended RCP runs
ens="r1i1p1"
yr_end=2005

# # declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# declare -a models=("CESM2-WACCM/") # extended RCP runs
# ens="r1i1p1f1"
# # declare -a ens="r1i1p1f2" # ensemble specification 
# # declare -a models=("MIROC-ES2L/") # extended RCP runs
# # declare -a ens="r3i1p1f2" # ensemble specification 
# # declare -a models=("GISS-E2-1-G/" "GISS-E2-1-H/") # extended RCP runs
# # declare -a ens="r4i1p1f2" # ensemble specification 
# # declare -a models=("UKESM1-0-LL/") # extended RCP runs
# yr_end=2014

sim="historical"
freq="Amon"
yr_begin=1860
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
