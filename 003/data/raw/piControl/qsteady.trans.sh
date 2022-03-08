#!/usr/bin/env bash
# set -euo pipefail

# load python environment
# module unload python
# module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# models=("ACCESS1-0/")
# yr_begin="0750"
# yr_end="0799"
# models=("CCSM4/")
# yr_begin="1250"
# yr_end="1299"
models=("MPI-ESM-LR/")
yr_begin="2800"
yr_end="2849"
sim="piControl"
freq="Amon"
ens="r1i1p1"
yr_span="${yr_begin}01-${yr_end}12"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    mkdir ${cwd}/${model}/qsteady

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/vqmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_ps=${cwd}/${model}/ps_${common}.timmean.nc
        full_va=${cwd}/${model}/va_${common}.nc
        full_hus=${cwd}/${model}/hus_${common}.nc
        full_beta=${cwd}/${model}/beta_${common}.timmean.nc

        for yr in $(seq -f "%04g" ${yr_begin} ${yr_end}); do
            echo $yr
            sel_span=${yr}01-${yr}12
            sel_common=${freq}_${model}_${sim}_${ens}_${sel_span}

            sel_va=${cwd}/${model}/qsteady/va_${sel_common}.nc
            sel_hus=${cwd}/${model}/qsteady/hus_${sel_common}.nc
            sel_vqmmc=${cwd}/${model}/qsteady/vqmmc_${sel_common}.nc
            sel_vqse=${cwd}/${model}/qsteady/vqse_${sel_common}.nc

            # cdo seldate,${yr}-01-01,${yr}-12-31 $full_ps $sel_ps
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_va $sel_va
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_hus $sel_hus

            python ${cwd}/make_qsteady.trans.py ${full_ps} ${sel_va} ${sel_hus} ${full_beta} ${sel_vqmmc} ${sel_vqse}

            # rm $sel_ps $sel_va $sel_hus
            rm $sel_va $sel_hus

            # exit 1 

        done

    fi

done
