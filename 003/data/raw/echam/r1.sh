#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

sim="echam"

declare -a models=("rp000191") # climate name
yr_span="0040_0059"

# declare -a models=("rp000184") # extended RCP runs
# # declare -a models=("rp000185" "rp000189") # extended RCP runs
# yr_span="0040_0139"

# declare -a models=("rp000190") # extended RCP runs
# # declare -a models=("rp000130" "rp000131" "rp000134" "rp000135") # extended RCP runs
# yr_span="0020_0039"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${model}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/r1_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_ra=${cwd}/${model}/ra_${common}.nc
        full_stf=${cwd}/${model}/stf_${common}.nc
        full_tend=${cwd}/${model}/tend_${common}.nc
        full_stgadv=${cwd}/${model}/stgadv_${common}.nc
        full_adv=${cwd}/${model}/adv_${common}.nc
        full_r1=${cwd}/${model}/r1_${common}.nc

        python ${cwd}/make_r1.py ${full_ra} ${full_stf} ${full_tend} ${full_stgadv} ${full_adv} ${full_r1}
    fi

done
