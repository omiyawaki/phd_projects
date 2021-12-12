#!/usr/bin/env bash
# set -euo pipefail

sims=("echr0001/")
# varnames=("vmmmc" "vmse" "vmte")
varnames=("mse")
freq="6h"
lev="pl"
yr_begin="1030"
yr_end="1030"

cwd=$(pwd) # save current working directory

for sim in ${sims[@]}; do

    sim=${sim%/}
    echo ${sim}

    for varname in ${varnames[@]}; do
        echo ${varname}

        for yr in $(seq $yr_begin $yr_end); do
            echo ${yr}

            common=${varname}_${freq}_${lev}_${sim}_${yr}

            cd ${cwd}/${sim}

            if ls $cwd/${sim}/${varname}_${common}.nc 1> /dev/null 2>&1; then # check if data is already there
                echo "merged ${varname} was already created. Skipping..."
            else
                cdo -O mergetime ${common}{01..12}.nc ${common}.nc
            fi
        done

    done        

done
