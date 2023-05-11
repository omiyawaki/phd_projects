#!/usr/bin/env bash
# set -euo pipefail

# declare -a models=("rp000200") # extended RCP runs
# yr_span="0040_0057"

declare -a models=("rp000201") # extended RCP runs
yr_span="0040_0083"

# declare -a sims=("rp000188" "rp000191" "rp000191b") # climate name
# yr_span="0040_0339"

sims=("rp000188" "rp000191f")
yr_span="0040_0252"

# sims=("rp000195")
# # sims=("rp000130" "rp000131" "rp000134" "rp000135")
# yr_span="0020_0039"

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
