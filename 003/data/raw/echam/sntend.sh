#!/usr/bin/env bash
# set -euo pipefail

# declare -a sims=("rp000190") # climate name
# yr_span="0001_0009"

# declare -a sims=("rp000188" "rp000191" "rp000191b") # climate name
# yr_span="0040_0339"

# sims=("rp000191b")
# # sims=("rp000185" "rp000188")
# yr_span="0040_0139"

sims=("rp000134" "rp000135" "rp000190" "rp000190a" "rp000190b")
# sims=("rp000134")
yr_span="0020_0039"

cwd=$(pwd) # save current working directory

for sim in ${sims[@]}; do

    echo ${sim}

    indir=${cwd}/${sim}
    cd $indir

    common=${sim}_${yr_span}


    if ls ${indir}/sntend_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "sni sntendency was already created. Skipping..."
    else
        full_sni=${indir}/sni_${common}.nc
        full_sntend=${indir}/sntend_${common}.nc

        # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_sni.py ${full_atm} ${full_sni}
        python ${cwd}/make_sntend.py ${full_sni} ${full_sntend}
    fi

done
