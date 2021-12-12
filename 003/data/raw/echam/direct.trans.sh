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
            common_ann=${freq}_${lev}_${sim}_${yr}
            common_atm=${sim}_${yr}
            common_ml=${freq}_ml_${sim}_${yr}${mon}
            common_ml_a=${freq}_ml_${sim}_${yr}
            common_aht=${freq}_${lev}_${sim}_${yr}

            # cdo seldate,${yr}-${mon}-01,${yr}-${mon}-31 /project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/${sim}/BOT_dm_${sim}_${yr}.nc ${cwd}/${sim}/BOT_dm_${sim}_${yr}${mon}.nc

            cd ${cwd}/${sim}

            if ls $cwd/${sim}/vmte_${common} 1> /dev/null 2>&1; then # check if data is already there
                echo "vmte was already created. Skipping..."
            else
                if [[ ${lev} == "si" ]]; then
                    full_atm=/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/${sim}/ATM_${common_ml}.nc
                    full_va=${cwd}/${sim}/va_${common}.nc
                else
                    full_atm=/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/${sim}/ATM_${common}.nc
                    # full_atm=${cwd}/${sim}/ATM_${common}.nc
                    full_atm_ann=${cwd}/${sim}/ATM_${common_atm}.yearmean.nc
                fi
                full_mses=${cwd}/${sim}/mses_${common_ml}.nc
                full_mses_a=${cwd}/${sim}/mses_${common_ml_a}.yearmean.nc
                full_mse=${cwd}/${sim}/mse_${common}.nc
                full_mse_ann=${cwd}/${sim}/mse_${common_ann}.yearmean.nc
                full_vmmmc=${cwd}/${sim}/vmmmc_${common}.nc
                full_vmse=${cwd}/${sim}/vmse_${common}.nc
                full_vmte=${cwd}/${sim}/vmte_${common}.nc
                full_vas=${cwd}/${sim}/vas_${common_ml}.nc
                full_vas_a=${cwd}/${sim}/vas_${common_ml_a}.yearmean.nc
                full_ps=/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/${sim}/ATM_${common_ml}.nc
                full_ps_a=${cwd}/${sim}/BOT_${sim}_${yr}.yearmean.nc
                full_aht=${cwd}/${sim}/aht_${common_aht}.monmean.nc

                # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_direct.trans.${lev}.py ${full_atm} ${full_mses} ${full_mse} ${full_vmmmc} ${full_vmse} ${full_vmte} ${full_vas} ${full_ps} ${full_aht}
                python ${cwd}/make_direct.trans.${lev}.py ${full_atm} ${full_mses} ${full_mse} ${full_vmmmc} ${full_vmse} ${full_vmte} ${full_vas} ${full_ps} ${full_aht} ${full_mse_ann} ${full_atm_ann} ${full_ps_a} ${full_mses_a} ${full_vas_a}

            fi

        done
    done        

done
