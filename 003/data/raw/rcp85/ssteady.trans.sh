#!/usr/bin/env bash
# set -euo pipefail

# load python environment
# module unload python
# module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
declare -a models=("CCSM4") # extended RCP runs
sim="rcp85"
freq="Amon"
ens="r1i1p1"
yr_begin=2006
yr_resume=2271
yr_end=2299
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
        full_ta=${cwd}/${model}/ta_${common}.nc
        full_zg=${cwd}/${model}/zg_${common}.nc
        full_beta=${cwd}/${model}/beta_${common}.timmean.nc

        for yr in $(seq -f "%04g" ${yr_resume} ${yr_end}); do
            echo $yr
            sel_span=${yr}01-${yr}12
            sel_common=${freq}_${model}_${sim}_${ens}_${sel_span}

            sel_va=${cwd}/${model}/ssteady/va_${sel_common}.nc
            sel_ta=${cwd}/${model}/ssteady/ta_${sel_common}.nc
            sel_zg=${cwd}/${model}/ssteady/zg_${sel_common}.nc
            sel_vsmmc=${cwd}/${model}/ssteady/vsmmc_${sel_common}.nc
            sel_vsse=${cwd}/${model}/ssteady/vsse_${sel_common}.nc

            # cdo seldate,${yr}-01-01,${yr}-12-31 $full_ps $sel_ps
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_va $sel_va
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_ta $sel_ta
            cdo seldate,${yr}-01-01,${yr}-12-31 $full_zg $sel_zg

            python ${cwd}/make_ssteady.trans.py ${full_ps} ${sel_va} ${sel_ta} ${sel_zg} ${full_beta} ${sel_vsmmc} ${sel_vsse}

            # rm $sel_ps $sel_va $sel_ta $sel_zg
            rm $sel_va $sel_ta $sel_zg

            # exit 1 

        done

    fi

done
