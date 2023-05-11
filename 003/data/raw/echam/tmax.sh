#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# # sims=("rp000130" "rp000131" "rp000134" "rp000135") # climate name
# sims=("rp000190") # climate name
# yr_span="0020_0039"

mean=""
sims=("rp000188" "rp000191" "rp000191b") # climate name
yr_span="0040_0339"

cwd=$(pwd)

for sim in ${sims[@]}; do
    
    indir=${cwd}/${sim}
    cd $indir

    common=${sim}_${yr_span}

    cd ${cwd}/${sim}

    if ls $cwd/${sim}/tmax_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "tmax was already created. Skipping..."
    else
        full_ps=${cwd}/${sim}/aps_${common}.nc
        full_ta=${cwd}/${sim}/t_${common}.nc
        full_tmax=${cwd}/${sim}/tmax_${common}.nc

        python ${cwd}/make_tmax.py ${full_ps} ${full_ta} ${full_tmax}
    fi

done
