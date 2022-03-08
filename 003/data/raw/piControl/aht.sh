#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# models=("ACCESS1-0/")
# yr_span="075001-079912"
# models=("CCSM4/")
# yr_span="125001-129912"
models=("MPI-ESM-LR/")
yr_span="280001-284912"
sim="piControl"
freq="Amon"
ens="r1i1p1"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/vmmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_rlut=${cwd}/${model}/rlut_${common}.nc
        full_rsut=${cwd}/${model}/rsut_${common}.nc
        full_rsdt=${cwd}/${model}/rsdt_${common}.nc
        full_rlus=${cwd}/${model}/rlus_${common}.nc
        full_rlds=${cwd}/${model}/rlds_${common}.nc
        full_rsds=${cwd}/${model}/rsds_${common}.nc
        full_rsus=${cwd}/${model}/rsus_${common}.nc
        full_hfls=${cwd}/${model}/hfls_${common}.nc
        full_hfss=${cwd}/${model}/hfss_${common}.nc
        full_tend=${cwd}/${model}/tend_${common}.nc
        full_aht=${cwd}/${model}/aht_${common}.nc

        # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_aht.py ${full_rlut} ${full_rsut} ${full_rsdt} ${full_rlus} ${full_rlds} ${full_rsds} ${full_rsus} ${full_hfls} ${full_hfss} ${full_tend} ${full_aht}
        python ${cwd}/make_aht.py ${full_rlut} ${full_rsut} ${full_rsdt} ${full_rlus} ${full_rlds} ${full_rsds} ${full_rsus} ${full_hfls} ${full_hfss} ${full_tend} ${full_aht}
    fi

done
