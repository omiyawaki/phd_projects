#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

declare -a ens="r1i1p1f1" # ensemble specification 
declare -a models=("MRI-ESM2-0/") # extended RCP runs

sim="ssp585"
freq="Amon"

yr_spans=("201501-206412" "206501-210012" "210101-215012" "215101-220012" "220101-225012" "225101-230012")

cwd=$(pwd) # save current working directory

for yr_span in ${yr_spans[@]}; do
    echo $yr_span

    for model in ${models[@]}; do
        
        model=${model%/}
        echo ${model}

        common=$(basename ${cwd}/${model}/mse_${freq}_${model}_${sim}_${ens}_*${yr_span}.nc)
        common=${common#mse_}
        common=${common%.nc}

        ref_common=$(basename ${cwd}/${model}/mse_${freq}_${model}_${sim}_${ens}_*201501-230012.nc)
        ref_common=${ref_common#mse_}
        ref_common=${ref_common%.nc}

        cd ${cwd}/${model}

        if ls $cwd/${model}/tend_${common} 1> /dev/null 2>&1; then # check if data is already there
            echo "tend was already created. Skipping..."
        else
            full_ps=${cwd}/${model}/ps_${ref_common}.timmean.nc
            full_mse=${cwd}/${model}/mse_${common}.nc
            full_beta=${cwd}/${model}/beta_${ref_common}.timmean.nc
            full_tend=${cwd}/${model}/tend_${common}.nc

            python ${cwd}/make_tend.py ${full_ps} ${full_mse} ${full_beta} ${full_tend}
        fi

    done

done
