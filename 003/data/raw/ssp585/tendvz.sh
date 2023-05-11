#!/usr/bin/env bash
# set -euo pipefail

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# models=("HadGEM2-ES/")
# models=("bcc-csm1-1/")
# declare -a models=("bcc-csm1-1/" "HadGEM2-ES/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# # declare -a models=("bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# yr_span="186001-200512"
# ens="r1i1p1"

ens="r1i1p1f1"
# declare -a models=("CESM2-WACCM/") # extended RCP runs
# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "CESM2-WACCM/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
declare -a ens="r1i1p1f2" # ensemble specification 
declare -a models=("MIROC-ES2L/") # extended RCP runs
# declare -a ens="r3i1p1f2" # ensemble specification 
# declare -a models=("GISS-E2-1-G/" "GISS-E2-1-H/") # extended RCP runs
# declare -a ens="r4i1p1f2" # ensemble specification 
# declare -a models=("UKESM1-0-LL/") # extended RCP runs
yr_span="201501-229912"

sim="ssp585"
freq="Amon"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    grid=$(basename ${cwd}/${model}/tendv_${freq}_${model}_${sim}_${ens}_*_${yr_span}.nc)
    grid=${grid#tendv_${freq}_${model}_${sim}_${ens}_}
    grid=${grid%_${yr_span}.nc}

    common=${freq}_${model}_${sim}_${ens}_${grid}_${yr_span}

    cd ${cwd}/${model}

    # if ls $cwd/${model}/tendv_${common} 1> /dev/null 2>&1; then # check if data is already there
    #     echo "tendv was already created. Skipping..."
    # else
        full_tendv=${cwd}/${model}/tendv_${common}.nc
        full_tendvz=${cwd}/${model}/tendv_${common}.zonmean.nc

        python ${cwd}/make_tendvz.py ${full_tendv} ${full_tendvz}
    # fi

done
