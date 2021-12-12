#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

cwd=$(pwd) # save current working directory

# for model in ${models[@]}; do
    
#     model=${model%/}
#     echo ${model}

#     common=${freq}_${model}_${sim}_${ens}_${yr_span}

#     cd ${cwd}/${model}

#     if ls $cwd/${model}/vmmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
#         echo "mmc and se was already created. Skipping..."
#     else
#         full_rlut=${cwd}/${model}/rlut_${common}.nc
#         full_rsut=${cwd}/${model}/rsut_${common}.nc
#         full_rsdt=${cwd}/${model}/rsdt_${common}.nc
#         full_rlus=${cwd}/${model}/rlus_${common}.nc
#         full_rlds=${cwd}/${model}/rlds_${common}.nc
#         full_rsds=${cwd}/${model}/rsds_${common}.nc
#         full_rsus=${cwd}/${model}/rsus_${common}.nc
#         full_hfls=${cwd}/${model}/hfls_${common}.nc
#         full_hfss=${cwd}/${model}/hfss_${common}.nc
#         full_tend=${cwd}/${model}/tend_${common}.nc
#         full_aht=${cwd}/${model}/aht_${common}.nc

#         # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_aht.py ${full_rlut} ${full_rsut} ${full_rsdt} ${full_rlus} ${full_rlds} ${full_rsds} ${full_rsus} ${full_hfls} ${full_hfss} ${full_tend} ${full_aht}
#         python ${cwd}/make_aht.py ${full_rlut} ${full_rsut} ${full_rsdt} ${full_rlus} ${full_rlds} ${full_rsds} ${full_rsus} ${full_hfls} ${full_hfss} ${full_tend} ${full_aht}
#     fi

# done

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
        echo ${yr}

        common=${freq}_${lev}_${sim}_${yr}

        cd ${cwd}/${sim}

        if ls $cwd/${sim}/aht_${common} 1> /dev/null 2>&1; then # check if data is already there
            echo "aht was already created. Skipping..."
        else
            full_bot=/project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/${sim}/BOT_${sim}_${yr}.nc
            full_tend=${cwd}/${sim}/tend_${common}.monmean.nc
            full_aht=${cwd}/${sim}/aht_${common}.monmean.nc

            # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_aht.py ${full_bot} ${full_tend} ${full_aht}
            python ${cwd}/make_aht.py ${full_bot} ${full_tend} ${full_aht}
        fi

    done        

done
