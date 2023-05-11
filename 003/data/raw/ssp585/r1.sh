#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a ens="r1i1p1f1" # ensemble specification 
# declare -a models=("CESM2-WACCM/") # extended RCP runs
# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/") # extended RCP runs
# declare -a models=("MRI-ESM2-0/") # extended RCP runs
# declare -a ens="r1i1p1f2" # ensemble specification 
# declare -a models=("MIROC-ES2L/") # extended RCP runs
declare -a ens="r3i1p1f2" # ensemble specification 
declare -a models=("GISS-E2-1-G/" "GISS-E2-1-H/") # extended RCP runs
# declare -a ens="r4i1p1f2" # ensemble specification 
# declare -a models=("UKESM1-0-LL/") # extended RCP runs

sim="ssp585"
freq="Amon"
yr_begin=2015
yr_end=2299
yr_span="${yr_begin}01-${yr_end}12"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=$(basename ${cwd}/${model}/rlut_${freq}_${model}_${sim}_${ens}_*${yr_span}.nc)
    common=${common#rlut_}
    common=${common%.nc}

    cd ${cwd}/${model}

    if ls $cwd/${model}/r1_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_ra=${cwd}/${model}/ra_${common}.nc
        full_stf=${cwd}/${model}/stf_${common}.nc
        full_tend=${cwd}/${model}/tend_${common}.nc
        full_stgadv=${cwd}/${model}/stgadv_${common}.nc
        full_adv=${cwd}/${model}/adv_${common}.nc
        full_r1=${cwd}/${model}/r1_${common}.nc

        python ${cwd}/make_r1.py ${full_ra} ${full_stf} ${full_tend} ${full_stgadv} ${full_adv} ${full_r1}
    fi

done
