#!/usr/bin/env bash
# set -euo pipefail

# declare -a models=("rp000188" "rp000191f") # extended RCP runs
# # declare -a models=("rp000185" "rp000189") # extended RCP runs
# yr_span="0040_0139"

# declare -a models=("rp000190" "rp000190b" "rp000134" "rp000135") # extended RCP runs
declare -a models=("rp000134" "rp000190f") # extended RCP runs
# declare -a models=("rp000130" "rp000131" "rp000134" "rp000135") # extended RCP runs
yr_span="0020_0039"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${model}_${yr_span}

    cd ${cwd}/${model}

    # if ls $cwd/${model}/hus_${common} 1> /dev/null 2>&1; then # check if data is already there
    #     echo "hus was already created. Skipping..."
    # else
        full_hus=${cwd}/${model}/q_${common}.nc
        full_ps=${cwd}/${model}/aps_${common}.timmean.nc
        full_husz=${cwd}/${model}/q_${common}.zonmean.nc

        python ${cwd}/make_husz.py ${full_hus} ${full_ps} ${full_husz}
    # fi

done
