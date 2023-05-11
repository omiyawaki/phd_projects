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

# declare -a models=("rp000191b") # extended RCP runs
# # declare -a models=("rp000185" "rp000189") # extended RCP runs
# yr_span="0040_0139"

declare -a models=("rp000190f") # extended RCP runs
# declare -a models=("rp000131" "rp000134") # extended RCP runs
# declare -a models=("rp000190a" "rp000190b" "rp000190d" "rp000130" "rp000131" "rp000134" "rp000135") # extended RCP runs
yr_span="0020_0039"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${model}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/swrefu_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else

        full_sradsd=${cwd}/${model}/sradsd_${common}.nc
        full_srafsd=${cwd}/${model}/srafsd_${common}.nc
        full_swrefu=${cwd}/${model}/swrefu_${common}.nc

        python ${cwd}/make_swrefu.py ${full_sradsd} ${full_srafsd} ${full_swrefu} 
    fi

done
