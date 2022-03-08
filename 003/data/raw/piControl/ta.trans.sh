#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
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

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/vtammc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_ps=${cwd}/${model}/ps_${common}.nc
        full_va=${cwd}/${model}/va_${common}.nc
        full_ta=${cwd}/${model}/ta_${common}.nc
        full_vtammc=${cwd}/${model}/vtammc_${common}.nc
        full_vtase=${cwd}/${model}/vtase_${common}.nc

        for yr in $(seq -f "%04g" ${yr_begin} ${yr_end}); do
            echo $yr
            sel_span=${yr}01-${yr}12
            sel_common=${freq}_${model}_${sim}_${ens}_${sel_span}

            sel_ps=${cwd}/${model}/ps_${sel_common}.nc
            sel_va=${cwd}/${model}/va_${sel_common}.nc
            sel_ta=${cwd}/${model}/ta_${sel_common}.nc
            sel_vtammc=${cwd}/${model}/vtammc_${sel_common}.nc
            sel_vtase=${cwd}/${model}/vtase_${sel_common}.nc

            cdo seldate,${yr}-01-01,${yr}-12-31 $full_ps $sel_ps
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_va $sel_va
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_ta $sel_ta

            python ${cwd}/make_ta.trans.py ${sel_ps} ${sel_va} ${sel_ta} ${sel_vtammc} ${sel_vtase}

            rm $sel_ps $sel_va $sel_ta

            exit 1 

        done

        # combine all years into a single file
        rm $full_vtammc $full_vtase
        cdo mergetime vtammc_*12.nc $full_vtammc
        cdo mergetime vtase_*12.nc $full_vtase

    fi

done
