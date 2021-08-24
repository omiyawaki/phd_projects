#!/usr/bin/env bash
# set -euo pipefail

sims=("echr0001/")
freq="6h"
yr_begin="1030"
mn_begin="01"
yr_end="1030"
mn_end="12"

cwd=$(pwd) # save current working directory

for sim in ${sims[@]}; do

    sim=${sim%/}
    echo ${sim}

    for yr in $(seq $yr_begin $yr_end); do
        for mon in $(seq -f "%02g" ${mn_begin} ${mn_end}); do
            echo ${yr}-${mon}

            common=${freq}_ml_${sim}_${yr}${mon}

            cd ${cwd}/${sim}

            if ls $cwd/${sim}/mses_${common} 1> /dev/null 2>&1; then # check if data is already there
                echo "mses was already created. Skipping..."
            else
                full_mse=${cwd}/${sim}/mse_${common}.nc
                full_mses=${cwd}/${sim}/mses_${common}.nc

                # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_mses.py ${full_mse} ${full_mses}
                python ${cwd}/make_mses.py ${full_mse} ${full_mses}
            fi

        done

    done        

done
