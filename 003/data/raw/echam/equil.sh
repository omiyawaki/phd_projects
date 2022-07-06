#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

sim="echam"

# declare -a models=("rp000190") # extended RCP runs
# # declare -a models=("rp000185" "rp000189") # extended RCP runs
# yr_span="0001_0009"

# declare -a models=("rp000134") # extended RCP runs
# yr_span="0001_0049"

# declare -a models=("rp000184") # extended RCP runs
# # declare -a models=("rp000185" "rp000189") # extended RCP runs
# yr_span="0040_0139"

declare -a models=("rp000135") # extended RCP runs
# declare -a models=("rp000131" "rp000134") # extended RCP runs
# declare -a models=("rp000130" "rp000131" "rp000134" "rp000135") # extended RCP runs
yr_span="0001_0039"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    plotdir=${cwd}/${model}/plot
    mkdir $plotdir

    common=${model}_${yr_span}

    cd ${cwd}/${model}

    full_ftoa=${cwd}/${model}/ftoa_${common}.gmean.nc
    plotname=$plotdir/check_equil_${common}.pdf

    python ${cwd}/check_equil.py ${full_ftoa} ${plotname}

done
