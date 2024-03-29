#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
declare -a models=("bcc-csm1-1") # extended RCP runs
sim="rcp85"
freq="Amon"
ens="r1i1p1"
yr_span="200601-229912"

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
        cp $full_aht ${full_aht%.nc}.zonmean.nc
        cdo zonmean $full_aht $full_aht.tmp
        mv $full_aht.tmp $full_aht
    fi

done
