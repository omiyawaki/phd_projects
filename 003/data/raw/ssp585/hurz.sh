#!/usr/bin/env bash
# set -euo pipefail

# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "CESM2-WACCM/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
declare -a models=("CESM2-WACCM/") # extended RCP runs
ens="r1i1p1f1"
# declare -a ens="r1i1p1f2" # ensemble specification 
# declare -a models=("MIROC-ES2L/") # extended RCP runs
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

    grid=$(basename ${cwd}/${model}/rlut_${freq}_${model}_${sim}_${ens}_*_${yr_span}.nc)
    grid=${grid#rlut_${freq}_${model}_${sim}_${ens}_}
    grid=${grid%_${yr_span}.nc}

    common=${freq}_${model}_${sim}_${ens}_${grid}_${yr_span}

    cd ${cwd}/${model}

    # if ls $cwd/${model}/hur_${common} 1> /dev/null 2>&1; then # check if data is already there
    #     echo "hur was already created. Skipping..."
    # else
        full_hur=${cwd}/${model}/hur_${common}.nc
        full_ps=${cwd}/${model}/ps_${common}.timmean.nc
        full_hurz=${cwd}/${model}/hur_${common}.zonmean.nc

        python ${cwd}/make_hurz.py ${full_hur} ${full_ps} ${full_hurz}
    # fi

done
