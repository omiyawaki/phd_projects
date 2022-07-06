#!/usr/bin/env bash
# set -euo pipefail

# load python environment
# module unload python
# module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
declare -a models=("IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
sim="rcp85"
freq="Amon"
ens="r1i1p1"
yr_begin=2006
yr_end=2299
yr_stop=2299
yr_span="${yr_begin}01-${yr_end}12"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    mkdir ${cwd}/${model}/vmmmc_lu

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
        full_lbeta=${cwd}/${model}/lbeta_${common}.nc
        full_ubeta=${cwd}/${model}/ubeta_${common}.nc

        for yr in $(seq -f "%04g" ${yr_begin} ${yr_stop}); do
            echo $yr
            sel_span=${yr}01-${yr}12
            sel_common=${freq}_${model}_${sim}_${ens}_${sel_span}

            # sel_ps=${cwd}/${model}/vmmmc_lu/ps_${sel_common}.nc
            sel_va=${cwd}/${model}/vmmmc_lu/va_${sel_common}.nc
            sel_mse=${cwd}/${model}/vmmmc_lu/mse_${sel_common}.nc
            sel_lbeta=${cwd}/${model}/vmmmc_lu/lbeta_${sel_common}.nc
            sel_ubeta=${cwd}/${model}/vmmmc_lu/ubeta_${sel_common}.nc
            sel_lvmmmc=${cwd}/${model}/vmmmc_lu/lvmmmc_${sel_common}.nc
            sel_uvmmmc=${cwd}/${model}/vmmmc_lu/uvmmmc_${sel_common}.nc

            # cdo seldate,${yr}-01-01,${yr}-12-31 $full_ps $sel_ps
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_va $sel_va
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_mse $sel_mse
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_lbeta $sel_lbeta
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_ubeta $sel_ubeta

            python ${cwd}/make_vmmmc.lu.py ${full_ps} ${sel_va} ${sel_mse} ${full_beta} ${sel_lbeta} ${sel_ubeta} ${sel_lvmmmc} ${sel_uvmmmc}

            # rm $sel_ps $sel_va $sel_mse
            rm $sel_va $sel_mse $sel_lbeta $sel_ubeta

            # exit 1 

        done

    fi

done
