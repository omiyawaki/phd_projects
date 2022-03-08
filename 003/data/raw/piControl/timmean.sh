#!/usr/bin/env bash
# set -euo pipefail

# load python environment
# module unload python
# module load python

declare -a varnames=("aht" "vmmmc" "vmse" "vmte")
# declare -a varnames=("qaht")
# declare -a varnames=("ps")

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# models=("ACCESS1-0/")
# yr_begin="0750"
# yr_end="0799"
# models=("CCSM4/")
# yr_begin="1250"
# yr_end="1299"
models=("MPI-ESM-LR/")
yr_begin="2800"
yr_end="2849"
sim="piControl"
freq="Amon"
ens="r1i1p1"
yr_span="${yr_begin}01-${yr_end}12"

cwd=$(pwd) # save current working directory

for varname in ${varnames[@]}; do
    for model in ${models[@]}; do
        
        model=${model%/}
        echo ${model}

        common=${freq}_${model}_${sim}_${ens}_${yr_span}

        cd ${cwd}/${model}

        # if ls $cwd/${model}/${varname}_${common}.timmean.nc 1> /dev/null 2>&1; then # check if data is already there
        # if 0; then
        #     echo "${varname} timmean was already taken Skipping..."

            full_in=${cwd}/${model}/${varname}_${common}.nc
            full_tm=${cwd}/${model}/${varname}_${common}.timmean.nc

            cdo -O timmean ${full_in} ${full_tm}

        # fi
    done
done
