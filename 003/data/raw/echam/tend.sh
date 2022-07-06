#!/usr/bin/env bash
# set -euo pipefail

# declare -a sims=("rp000190") # climate name
# yr_span="0001_0009"

# declare -a sims=("rp000134") # climate name
# yr_span="0001_0049"

# sims=("rp000184")
# # sims=("rp000185" "rp000188")
# yr_span="0040_0139"

sims=("rp000190")
# sims=("rp000130" "rp000131" "rp000134" "rp000135")
yr_span="0020_0039"

cwd=$(pwd) # save current working directory

for sim in ${sims[@]}; do

    echo ${sim}

    indir=${cwd}/${sim}
    cd $indir

    common=${sim}_${yr_span}


    if ls ${indir}/tend_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "MSE tendency was already created. Skipping..."
    else
        full_mse=${indir}/mse_${common}.nc
        full_ps=${indir}/aps_${common}.nc
        full_tend=${indir}/tend_${common}.nc

        # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_mse.py ${full_atm} ${full_mse}
        python ${cwd}/make_tend.py ${full_mse} ${full_ps} ${full_tend}
    fi

done
