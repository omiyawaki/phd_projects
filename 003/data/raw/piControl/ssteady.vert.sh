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

    mkdir ${cwd}/${model}/ssteady

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/vsmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_ps=${cwd}/${model}/ps_${common}.timmean.nc
        full_va=${cwd}/${model}/va_${common}.nc
        full_dse=${cwd}/${model}/dse_${common}.nc
        full_beta=${cwd}/${model}/beta_${common}.timmean.nc

        for yr in $(seq -f "%04g" ${yr_begin} ${yr_end}); do
            echo $yr
            sel_span=${yr}01-${yr}12
            sel_common=${freq}_${model}_${sim}_${ens}_${sel_span}

            sel_va=${cwd}/${model}/ssteady/va_${sel_common}.nc
            sel_dse=${cwd}/${model}/ssteady/dse_${sel_common}.nc
            sel_vsmmc=${cwd}/${model}/ssteady/vvsmmc_${sel_common}.nc
            sel_vsse=${cwd}/${model}/ssteady/vvsse_${sel_common}.nc

            # cdo seldate,${yr}-01-01,${yr}-12-31 $full_ps $sel_ps
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_va $sel_va
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_dse $sel_dse

            python ${cwd}/make_ssteady.vert.py ${sel_va} ${sel_dse} ${full_beta} ${sel_vsmmc} ${sel_vsse}

            rm $sel_va $sel_dse

            exit 1 

        done

    fi

done
