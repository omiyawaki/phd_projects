#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
models=("MPI-ESM-LR/")
sim="historical"
freq="day"
ens="r1i1p1"
yr_span="186001-200512"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    mkdir ${cwd}/${model}/ke

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd /project2/tas1/ockham/data9/tas/CMIP5_RAW/${model}/${sim}/atmos/${freq}/
    basedir=$(pwd)
    
    files=$(ls ${basedir}/ua/r1i1p1/ua_*.nc)

    for file in $files; do

        full_eke=${cwd}/${model}/ke/$(basename ${file/ua_/eke_})
        full_eke=${full_eke/day/Amon}
        full_veke=${full_eke/eke/veke}
        full_ua=${file}
        full_va=${file//ua/va}

        full_ps=${cwd}/${model}/ps_*${yr_span}.nc

        # extract ps data for same time segment
        time_span=${file: -20}
        yr_begin=${time_span: 0:4}
        mn_begin=${time_span: 4:2}
        dy_begin=${time_span: 6:2}
        yr_end=${time_span: 9:4}
        mn_end=${time_span: 13:2}
        dy_end=${time_span: 15:2}

        sel_ps=${cwd}/${model}/te/ps_${yr_begin}_${yr_end}.nc

        if ls $sel_ps 1> /dev/null 2>&1; then # check if data is already there
            echo "select ps already created. Skipping..."
        else
            cdo seldate,$yr_begin-$mn_begin-$dy_begin,$yr_end-$mn_end-$dy_end $full_ps $sel_ps
        fi

        echo ${time_span}

        # if ls $full_vmte 1> /dev/null 2>&1; then # check if data is already there
        #     echo "te was already created. Skipping..."
        # else

            # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_te.py ${full_vmmmc} ${full_vmse} ${full_aht} ${full_vmte}
            python ${cwd}/make_eke.py ${sel_ps} ${full_ua} ${full_va} ${full_eke} ${full_veke} 
        # fi

    done

done
