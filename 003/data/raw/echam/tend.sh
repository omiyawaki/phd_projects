#!/usr/bin/env bash
# set -euo pipefail

sims=("echr0001/")
freq="6h"
lev="ml"
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
    
            common=${freq}_${lev}_${sim}_${yr}${mon}

            cd ${cwd}/${sim}

            if ls $cwd/${sim}/tend_${common} 1> /dev/null 2>&1; then # check if data is already there
                echo "tend was already created. Skipping..."
            else
                full_atm=/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/${sim}/ATM_${common}.nc
                full_mse=${cwd}/${sim}/mse_${common}.nc
                full_tend=${cwd}/${sim}/tend_${common}.nc

                # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_tend.py ${full_atm} ${full_mse} ${full_tend}
                python ${cwd}/make_tend.py ${full_atm} ${full_mse} ${full_tend}
            fi

        done
    done        

done
