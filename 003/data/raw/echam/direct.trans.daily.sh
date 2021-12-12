#!/usr/bin/env bash
# set -euo pipefail

sims=("echr0001/")
freq="6h"
lev="pl"
yr_begin="1030"
mn_begin="01"
yr_end="1030"
mn_end="01"

cwd=$(pwd) # save current working directory

for sim in ${sims[@]}; do

    sim=${sim%/}
    echo ${sim}

    for yr in $(seq $yr_begin $yr_end); do
        for mon in $(seq -f "%02g" ${mn_begin} ${mn_end}); do
            echo ${yr}-${mon}
    
            common=${freq}_${lev}_${sim}_${yr}${mon}
            common_ml=${freq}_ml_${sim}_${yr}${mon}
            common_bot=${freq}_${sim}_${yr}${mon}

            # cdo seldate,${yr}-${mon}-01,${yr}-${mon}-31 /project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/${sim}/BOT_dm_${sim}_${yr}.nc ${cwd}/${sim}/BOT_dm_${sim}_${yr}${mon}.nc

            cd ${cwd}/${sim}

            if ls $cwd/${sim}/vmte_${common} 1> /dev/null 2>&1; then # check if data is already there
                echo "vmte was already created. Skipping..."
            else
                full_atm=${cwd}/${sim}/ATM_${common}.nc
                full_mse=${cwd}/${sim}/mse_${common}.nc
                full_vmte=${cwd}/${sim}/vmte_${common}.nc
                if [ $freq == 'dm' ]; then
                    full_ps=${cwd}/${sim}/BOT_${common_bot}.nc
                else
                    full_ps=${cwd}/${sim}/ATM_${common_ml}.nc
                fi

                python ${cwd}/make_direct.trans.${lev}.dm.py ${full_atm} ${full_mse} ${full_vmte} ${full_ps}

            fi

        done
    done        

done
