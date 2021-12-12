#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
sims=("echr0001/")
freq="6h"
lev="pl"
yr_begin="1030"
mn_begin="01"
yr_end="1030"
mn_end="01"


cwd=$(pwd) # save current working directory

for sim in ${sims[@]}; do
    
    sim=${sim%/}
    echo ${sim}

    for yr in $(seq $yr_begin $yr_end); do
        for mon in $(seq -f "%02g" ${mn_begin} ${mn_end}); do
            echo ${yr}-${mon}

            common=${freq}_${lev}_${sim}_${yr}${mon}
            common_aht=${freq}_${lev}_${sim}_${yr}

            cd ${cwd}/${sim}

            if ls $cwd/${sim}/vmmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
                echo "mmc and se was already created. Skipping..."
            else
                full_vmmmc=${cwd}/${sim}/vmmmc_${common}.nc
                full_vmse=${cwd}/${sim}/vmse_${common}.nc
                full_aht=${cwd}/${sim}/aht_${common_aht}.monmean.nc
                full_vmte=${cwd}/${sim}/vmte_${common}.nc

                # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_te.py ${full_vmmmc} ${full_vmse} ${full_aht} ${full_vmte}
                python ${cwd}/make_te.py ${full_vmmmc} ${full_vmse} ${full_aht} ${full_vmte}
            fi

        done

    done

done
