#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

sim="echam"

# declare -a models=("rp000190") # climate name
# yr_span="0001_0009"

# declare -a models=("rp000134") # climate name
# yr_span="0001_0049"

# declare -a models=("rp000184") # extended RCP runs
# # declare -a models=("rp000185" "rp000189") # extended RCP runs
# yr_span="0040_0139"

declare -a models=("rp000190") # extended RCP runs
# declare -a models=("rp000130" "rp000131" "rp000134" "rp000135") # extended RCP runs
yr_span="0020_0039"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${model}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/r1_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else

        full_trad0=${cwd}/${model}/trad0_${common}.nc
        full_srad0=${cwd}/${model}/srad0_${common}.nc
        full_trads=${cwd}/${model}/trads_${common}.nc
        full_srads=${cwd}/${model}/srads_${common}.nc
        full_ra=${cwd}/${model}/ra_${common}.nc

        python ${cwd}/make_ra.py ${full_trad0} ${full_srad0} ${full_trads} ${full_srads} ${full_ra}
    fi

done
