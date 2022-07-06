#!/usr/bin/env bash
# set -euo pipefail

models=("CCSM4/")
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

    # if ls $cwd/${model}/ta_${common} 1> /dev/null 2>&1; then # check if data is already there
    #     echo "ta was already created. Skipping..."
    # else
        full_ta=${cwd}/${model}/ta_${common}.nc
        full_ps=${cwd}/${model}/ps_${common}.timmean.nc
        full_taz=${cwd}/${model}/ta_${common}.zonmean.nc

        python ${cwd}/make_taz_ccsm.py ${full_ta} ${full_ps} ${full_taz}
    # fi

done
