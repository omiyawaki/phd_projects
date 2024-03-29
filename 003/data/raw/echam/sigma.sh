#!/usr/bin/env bash
# set -euo pipefail

sims=("echr0001/")
freq="6h"
lev="pl"
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
            common_ml=${freq}_ml_${sim}_${yr}${mon}
            common_si=${freq}_si_${sim}_${yr}${mon}

            # cdo seldate,${yr}-${mon}-01,${yr}-${mon}-31 /project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/${sim}/BOT_dm_${sim}_${yr}.nc ${cwd}/${sim}/BOT_dm_${sim}_${yr}${mon}.nc

            cd ${cwd}/${sim}

            if ls $cwd/${sim}/vmte_${common} 1> /dev/null 2>&1; then # check if data is already there
                echo "vmte was already created. Skipping..."
            else
                full_atm=/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/${sim}/ATM_${common}.nc
                full_mses=${cwd}/${sim}/mses_${common_ml}.nc
                full_mse=${cwd}/${sim}/mse_${common}.nc
                full_vas=${cwd}/${sim}/vas_${common_ml}.nc
                full_ps=/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/${sim}/ATM_${common_ml}.nc
                full_vasi=${cwd}/${sim}/va_${common_si}.nc
                full_msesi=${cwd}/${sim}/mse_${common_si}.nc

                # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_sigma.py ${full_atm} ${full_mses} ${full_mse} ${full_vas} ${full_ps} ${full_vasi} ${full_msesi}
                python ${cwd}/make_sigma.py ${full_atm} ${full_mses} ${full_mse} ${full_vas} ${full_ps} ${full_vasi} ${full_msesi}

            fi

        done
    done        

done
