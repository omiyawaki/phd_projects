#!/usr/bin/env bash
# set -euo pipefail

# load python environment
# module unload python
# module load python

varnames=("aps")

# declare -a models=("rp000200") # extended RCP runs
# yr_span="0040_0057"

declare -a models=("rp000201") # extended RCP runs
yr_span="0040_0083"

# declare -a models=("rp000188" "rp000191f") # extended RCP runs
# # declare -a models=("rp000185" "rp000189") # extended RCP runs
# yr_span="0040_0252"

# # declare -a models=("rp000190" "rp000190b" "rp000134" "rp000135") # extended RCP runs
# declare -a models=("rp000134" "rp000190f") # extended RCP runs
# # declare -a models=("rp000130" "rp000131" "rp000134" "rp000135") # extended RCP runs
# yr_span="0020_0039"

cwd=$(pwd) # save current working directory

for varname in ${varnames[@]}; do
    for model in ${models[@]}; do
        
        model=${model%/}
        echo ${model}

        common=${model}_${yr_span}

        cd ${cwd}/${model}

        # if ls $cwd/${model}/${varname}_${common}.timmean.nc 1> /dev/null 2>&1; then # check if data is already there
        #     echo "${varname} timmean was already taken Skipping..."
        # else
            full_in=${cwd}/${model}/${varname}_${common}.nc
            full_tm=${cwd}/${model}/${varname}_${common}.timmean.nc

            cdo -O timmean ${full_in} ${full_tm}

        # fi
    done
done
