#!/usr/bin/env bash
# set -euo pipefail

# load python environment
# module unload python
# module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# declare -a models=("CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
sim="historical"
freq="Amon"
ens="r1i1p1"
yr_begin=1860
yr_end=2005
yr_span="${yr_begin}01-${yr_end}12"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    mkdir ${cwd}/${model}/steady

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/vmmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_ps=${cwd}/${model}/ps_${common}.timmean.nc
        # full_ps=${cwd}/${model}/ps_${common}.nc
        full_va=${cwd}/${model}/va_${common}.nc
        full_mse=${cwd}/${model}/mse_${common}.nc
        full_beta=${cwd}/${model}/beta_${common}.timmean.nc

        for yr in $(seq -f "%04g" ${yr_begin} ${yr_end}); do
            echo $yr
            sel_span=${yr}01-${yr}12
            sel_common=${freq}_${model}_${sim}_${ens}_${sel_span}

            # sel_ps=${cwd}/${model}/steady/ps_${sel_common}.nc
            sel_va=${cwd}/${model}/steady/va_${sel_common}.nc
            sel_mse=${cwd}/${model}/steady/mse_${sel_common}.nc
            sel_vmmmc=${cwd}/${model}/steady/vmmmc_${sel_common}.nc
            sel_vmse=${cwd}/${model}/steady/vmse_${sel_common}.nc

            # cdo seldate,${yr}-01-01,${yr}-12-31 $full_ps $sel_ps
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_va $sel_va
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_mse $sel_mse

            python ${cwd}/make_steady.trans.py ${full_ps} ${sel_va} ${sel_mse} ${full_beta} ${sel_vmmmc} ${sel_vmse}

            # rm $sel_ps $sel_va $sel_mse
            rm $sel_va $sel_mse

            # exit 1 

        done

    fi

done
