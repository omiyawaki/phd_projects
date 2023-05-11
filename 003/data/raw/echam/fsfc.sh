#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

sim="echam"

# declare -a models=("rp000190b") # extended RCP runs
# # declare -a models=("rp000185" "rp000189") # extended RCP runs
# yr_span="0001_0009"

# declare -a models=("rp000188" "rp000191" "rp000191b") # extended RCP runs
# yr_span="0040_0339"

declare -a models=("rp000188" "rp000191f") # extended RCP runs
# declare -a models=("rp000185" "rp000189") # extended RCP runs
yr_span="0040_0252"

# declare -a models=("rp000190e" "rp000190f") # extended RCP runs
# # declare -a models=("rp000131" "rp000134") # extended RCP runs
# # declare -a models=("rp000130" "rp000131" "rp000134" "rp000135") # extended RCP runs
# yr_span="0020_0039"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${model}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/fsfc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else

        full_hfls=${cwd}/${model}/ahfl_${common}.nc
        full_hfss=${cwd}/${model}/ahfs_${common}.nc
        full_trads=${cwd}/${model}/trads_${common}.nc
        full_srads=${cwd}/${model}/srads_${common}.nc
        full_fsfc=${cwd}/${model}/fsfc_${common}.nc

        python ${cwd}/make_fsfc.py ${full_hfls} ${full_hfss} ${full_trads} ${full_srads} ${full_fsfc}
    fi

done
