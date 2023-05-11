#!/usr/bin/env bash
# set -euo pipefail

# # declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# # models=("HadGEM2-ES/")
# models=("CCSM4/")
# # declare -a models=("bcc-csm1-1/" "HadGEM2-ES/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# # declare -a models=("bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# ens="r1i1p1"
# yr_span="186001-200512"

# declare -a models=("CESM2-WACCM/" "ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# ens="r1i1p1f1"
# declare -a ens="r1i1p1f2" # ensemble specification 
# declare -a models=("MIROC-ES2L/") # extended RCP runs
declare -a ens="r3i1p1f2" # ensemble specification 
declare -a models=("GISS-E2-1-G/" "GISS-E2-1-H/") # extended RCP runs
# declare -a ens="r4i1p1f2" # ensemble specification 
# declare -a models=("UKESM1-0-LL/") # extended RCP runs
yr_span="186001-201412"

sim="historical"
freq="Amon"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    grid=$(basename ${cwd}/${model}/rlut_${freq}_${model}_${sim}_${ens}_*_${yr_span}.nc)
    grid=${grid#rlut_${freq}_${model}_${sim}_${ens}_}
    grid=${grid%_${yr_span}.nc}

    common=${freq}_${model}_${sim}_${ens}_${grid}_${yr_span}

    cd ${cwd}/${model}

    # if ls $cwd/${model}/hus_${common} 1> /dev/null 2>&1; then # check if data is already there
    #     echo "hus was already created. Skipping..."
    # else
        full_hus=${cwd}/${model}/hus_${common}.nc
        full_ps=${cwd}/${model}/ps_${common}.timmean.nc
        full_husz=${cwd}/${model}/hus_${common}.zonmean.nc

        python ${cwd}/make_husz.py ${full_hus} ${full_ps} ${full_husz}
    # fi

done
