#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
models=("MPI-ESM-LR/")
sim="piControl"
freq="Amon"
ens="r1i1p1"
yr_begin="2800"
yr_end="2849"
yr_span="${yr_begin}01-${yr_end}12"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/vmmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_ps=${cwd}/${model}/ps_${common}.nc
        full_va=${cwd}/${model}/va_${common}.nc
        full_vas=${cwd}/${model}/vas_${common}.nc
        full_mse=${cwd}/${model}/mse_${common}.nc
        full_mses=${cwd}/${model}/mses_${common}.nc
        full_vmmmc=${cwd}/${model}/vmmmc_${common}.nc
        full_vmse=${cwd}/${model}/vmse_${common}.nc

        for yr in $(seq -f "%04g" ${yr_begin} ${yr_end}); do
            echo $yr
            sel_span=${yr}01-${yr}12
            sel_common=${freq}_${model}_${sim}_${ens}_${sel_span}

            sel_ps=${cwd}/${model}/ps_${sel_common}.nc
            sel_va=${cwd}/${model}/va_${sel_common}.nc
            sel_vas=${cwd}/${model}/vas_${sel_common}.nc
            sel_mse=${cwd}/${model}/mse_${sel_common}.nc
            sel_mses=${cwd}/${model}/mses_${sel_common}.nc
            sel_vmmmc=${cwd}/${model}/vmmmc_${sel_common}.nc
            sel_vmse=${cwd}/${model}/vmse_${sel_common}.nc

            cdo seldate,${yr}-01-01,${yr}-12-31 $full_ps $sel_ps
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_va $sel_va
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_vas $sel_vas
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_mse $sel_mse
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_mses $sel_mses

            # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_steady.trans.py ${sel_ps} ${sel_va} ${sel_vas} ${sel_mse} ${sel_mses} ${sel_vmmmc} ${sel_vmse}
            python ${cwd}/make_steady.trans.py ${sel_ps} ${sel_va} ${sel_vas} ${sel_mse} ${sel_mses} ${sel_vmmmc} ${sel_vmse}

            rm $sel_ps $sel_va $sel_vas $sel_mse $sel_mses

        done

        # combine all years into a single file
        rm $full_vmmmc $full_vmse
        cdo mergetime vmmmc_*12.nc $full_vmmmc
        cdo mergetime vmse_*12.nc $full_vmse

    fi

done
