#!/usr/bin/env bash
# set -euo pipefail

sims=("rp000189/")
# freq="_6h"
# lev="_pl"
freq=""
lev=""
yr_begin="0040"
mn_begin="01"
yr_end="0041"
mn_end="12"

cwd=$(pwd) # save current working directory

for sim in ${sims[@]}; do

    sim=${sim%/}
    echo ${sim}

    for yr in $(seq $yr_begin $yr_end); do
        for mon in $(seq -f "%02g" ${mn_begin} ${mn_end}); do
            echo ${yr}-${mon}
    
            common=${freq}${lev}_${sim}_${yr}${mon}

            cd ${cwd}/${sim}

            if ls $cwd/${sim}/mse_${common} 1> /dev/null 2>&1; then # check if data is already there
                echo "MSE was already created. Skipping..."
            else
                full_atm=/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/${sim}/ATM_${common}.nc
                full_mse=${cwd}/${sim}/mse_${common}.nc

                # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_mse.py ${full_atm} ${full_mse}
                python ${cwd}/make_mse.py ${full_atm} ${full_mse}
            fi

        done
    done        

done
