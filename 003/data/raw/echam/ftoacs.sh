#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

sim="echam"

# declare -a models=("rp000190") # extended RCP runs
# # declare -a models=("rp000185" "rp000189") # extended RCP runs
# yr_span="0001_0009"

# declare -a models=("rp000135") # extended RCP runs
# yr_span="0001_0039"

# declare -a models=("rp000188" "rp000191" "rp000191b") # extended RCP runs
# yr_span="0040_0339"

declare -a models=("rp000188" "rp000191f") # extended RCP runs
# declare -a models=("rp000185" "rp000189") # extended RCP runs
yr_span="0040_0252"

# # declare -a models=("rp000190" "rp000190b" "rp000134" "rp000135") # extended RCP runs
# declare -a models=("rp000134" "rp000190f") # extended RCP runs
# # declare -a models=("rp000130" "rp000131" "rp000134" "rp000135") # extended RCP runs
# yr_span="0020_0039"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${model}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/ftoacs_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else

        full_traf0=${cwd}/${model}/traf0_${common}.nc
        full_sraf0=${cwd}/${model}/sraf0_${common}.nc
        full_ftoacs=${cwd}/${model}/ftoacs_${common}.nc

        python ${cwd}/make_ftoacs.py ${full_traf0} ${full_sraf0} ${full_ftoacs}
    fi

done
