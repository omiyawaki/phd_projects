#!/usr/bin/env bash
# set -euo pipefail

# declare -a sims=("rp000190") # climate name
# yr_span="0001_0009"

# declare -a sims=("rp000188" "rp000191" "rp000191b") # climate name
# yr_span="0040_0339"

# sims=("rp000191b")
# # sims=("rp000185" "rp000188")
# yr_span="0040_0139"

# sims=("rp000134" "rp000135" "rp000190" "rp000190a" "rp000190b")
sims=("rp000134")
yr_span="0020_0039"

cwd=$(pwd) # save current working directory

for sim in ${sims[@]}; do

    echo ${sim}

    indir=${cwd}/${sim}
    cd $indir

    common=${sim}_${yr_span}


    if ls ${indir}/snmelt_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "sni snmeltency was already created. Skipping..."
    else
        full_tsurf=${indir}/tsurf_${common}.nc
        full_tsi=${indir}/tsi_${common}.nc
        full_sni=${indir}/sni_${common}.nc
        full_snmelt=${indir}/snmelt_${common}.nc

        # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_sni.py ${full_atm} ${full_sni}
        python ${cwd}/make_snmelt.py ${full_tsurf} ${full_tsi} ${full_sni} ${full_snmelt}
    fi

done
