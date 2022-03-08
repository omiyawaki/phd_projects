#!/usr/bin/env bash
# set -euo pipefail

# load python environment
# module unload python
# module load python

# declare -a varnames=("aht" "vmmmc" "vmse" "vmte")
declare -a varnames=("ps")

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
sim="historical"
freq="Amon"
ens="r1i1p1"
yr_span="186001-200512"

cwd=$(pwd) # save current working directory

for varname in ${varnames[@]}; do
    for model in ${models[@]}; do
        
        model=${model%/}
        echo ${model}

        common=${freq}_${model}_${sim}_${ens}_${yr_span}

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
