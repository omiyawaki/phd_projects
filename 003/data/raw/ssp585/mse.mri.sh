#!/usr/bin/env bash
# set -euo pipefail

declare -a ens="r1i1p1f1" # ensemble specification 
declare -a models=("MRI-ESM2-0/") # extended RCP runs

# yr_spans=("201501-206412" "206501-210012" "210101-215012" "215101-220012" "220101-225012" "225101-230012")
yr_spans=("210101-215012" "215101-220012" "220101-225012" "225101-230012")

sim="ssp585"
freq="Amon"
# yr_begin=2015
# yr_end=2064
# yr_begin=2065
# yr_end=2100
# yr_begin=2101
# yr_end=2150
# yr_begin=2151
# yr_end=2200
# yr_begin=2201
# yr_end=2250
# yr_begin=2251
# yr_end=2300
# yr_span="${yr_begin}01-${yr_end}12"

cwd=$(pwd) # save current working directory

for yr_span in ${yr_spans[@]}; do
    echo $yr_span

    for model in ${models[@]}; do
        
        model=${model%/}
        echo ${model}

        common=$(basename ${cwd}/${model}/ta_${freq}_${model}_${sim}_${ens}_*${yr_span}.nc)
        common=${common#ta_}
        common=${common%.nc}

        cd ${cwd}/${model}

        if ls $cwd/${model}/mse_${common} 1> /dev/null 2>&1; then # check if data is already there
            echo "MSE was already created. Skipping..."
        else
            full_ta=${cwd}/${model}/ta_${common}.nc
            full_zg=${cwd}/${model}/zg_${common}.nc
            full_hus=${cwd}/${model}/hus_${common}.nc
            full_mse=${cwd}/${model}/mse_${common}.nc

            # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_mse.py ${full_ta} ${full_zg} ${full_hus} ${full_mse}
            python ${cwd}/make_mse.py ${full_ta} ${full_zg} ${full_hus} ${full_mse}
        fi

    done

done
