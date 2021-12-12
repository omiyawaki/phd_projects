#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

sims=("echr0001/")
freq="dm"
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

            full_vmte=${cwd}/${sim}/vmte_${common}.nc
            full_divte=${cwd}/${sim}/divte_${common}.nc

            python ${cwd}/make_fluxdiv.te.py ${full_vmte} ${full_divte}

        done
    done

done
