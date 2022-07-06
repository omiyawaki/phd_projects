#!/usr/bin/env bash
# set -euo pipefail

# load python environment
# module unload python
# module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
declare -a models=("bcc-csm1-1") # extended RCP runs
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

        for yr in $(seq -f "%04g" ${yr_begin} ${yr_stop}); do
            echo $yr
            sel_span=${yr}01-${yr}12
            sel_common=${freq}_${model}_${sim}_${ens}_${sel_span}

            sel_vmmmc=${cwd}/${model}/steady/vmmmc_${sel_common}.nc
            sel_lvmmmc=${cwd}/${model}/vmmmc_lu/lvmmmc_${sel_common}.nc
            sel_uvmmmc=${cwd}/${model}/vmmmc_lu/uvmmmc_${sel_common}.nc
            sel_res=${cwd}/${model}/vmmmc_lu/res_${sel_common}.nc

            python ${cwd}/check_vmmmc.lu.py ${sel_vmmmc} ${sel_lvmmmc} ${sel_uvmmmc} ${sel_res}

            # exit 1 

        done

    fi

done
