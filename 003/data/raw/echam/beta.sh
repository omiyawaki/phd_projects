#!/usr/bin/env bash
# set -euo pipefail

# load python environment
# module unload python
# module load python

declare -a models=("rp000200") # extended RCP runs
yr_span="0040_0057"

# declare -a models=("rp000201") # extended RCP runs
# yr_span="0040_0083"

# declare -a models=("rp000188" "rp000191f") # extended RCP runs
# # declare -a models=("rp000185" "rp000189") # extended RCP runs
# yr_span="0040_0252"

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

    if ls $cwd/${model}/beta_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "beta was already created. Skipping..."
    else
        full_ps=${cwd}/${model}/aps_${common}.timmean.nc
        full_ta=${cwd}/${model}/t_${common}.nc
        full_beta=${cwd}/${model}/beta_${common}.timmean.nc

        python ${cwd}/make_beta.py ${full_ps} ${full_ta} ${full_beta}

    fi

done
