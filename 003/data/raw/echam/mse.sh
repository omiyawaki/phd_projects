#!/usr/bin/env bash
# set -euo pipefail

declare -a sims=("rp000200") # climate name
yr_span="0040_0057"

# declare -a sims=("rp000201") # climate name
# yr_span="0040_0083"

# declare -a sims=("rp000191b") # climate name
# # declare -a sims=("rp000188" "rp000191" "rp000191b") # climate name
# yr_span="0040_0339"

# sims=("rp000188" "rp000191f")
# # sims=("rp000188")
# yr_span="0040_0252"

# sims=("rp000195")
# # sims=("rp000130" "rp000131" "rp000134" "rp000135")
# yr_span="0020_0039"

cwd=$(pwd) # save current working directory

for sim in ${sims[@]}; do

    echo ${sim}

    indir=${cwd}/${sim}
    cd $indir

    common=${sim}_${yr_span}


    if ls ${indir}/mse_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "MSE was already created. Skipping..."
    else
        full_t=${indir}/t_${common}.nc
        full_q=${indir}/q_${common}.nc
        full_z=${indir}/geopoth_${common}.nc
        full_mse=${indir}/mse_${common}.nc

        # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_mse.py ${full_atm} ${full_mse}
        python ${cwd}/make_mse.py ${full_t} ${full_q} ${full_z} ${full_mse}
    fi

done
