#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
models=("IPSL-CM5A-LR/")
sim="historical"
freq="day"
ens="r1i1p1"
yr_span="185001-200512"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    mkdir ${cwd}/${model}/te

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd /project2/tas1/ockham/data9/tas/CMIP5_RAW/${model}/${sim}/atmos/${freq}/
    basedir=$(pwd)
    
    files=$(ls ${basedir}/hus/r1i1p1/hus_*.nc)

    for file in $files; do

        full_vmte=${cwd}/${model}/te/$(basename ${file/hus_/vmte_})
        full_vmte=${full_vmte/day/Amon}
        full_hus=${file}
        full_va=${file//hus/va}
        full_ta=${file//hus/ta}
        full_zg=${file//hus/zg}

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

        # if ls $full_vmte 1> /dev/null 2>&1; then # check if data is already there
        #     echo "te was already created. Skipping..."
        # else

            # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_te.py ${full_vmmmc} ${full_vmse} ${full_aht} ${full_vmte}
            python ${cwd}/make_te.py ${sel_ps} ${full_va} ${full_ta} ${full_hus} ${full_zg} ${full_vmte}

        # fi

    done

done
