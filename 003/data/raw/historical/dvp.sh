#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# declare -a models=("HadGEM2-ES/" "bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
models=("HadGEM2-ES/")
# models=("IPSL-CM5A-LR/")
sim="historical"
freq="Amon"
ens="r1i1p1"
yr_span="186001-200512"
lev="92500"
mean=".ymonmean-30"
llo=70
clat=80
lhi=90

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/dvp${lev}_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "dvp was already created. Skipping..."
    else
        full_vmte=${cwd}/${model}/vmte_${common}.lat_${clat}${mean}.nc
        full_gmse=${cwd}/${model}/gmse${lev}_${common}.zonmean.amean_${llo}_${lhi}${mean}.nc
        full_diffv=${cwd}/${model}/dvp${lev}_${common}.lat_${clat}${mean}.nc

        python ${cwd}/make_dvp.py ${full_vmte} ${full_gmse} ${full_diffv} ${clat}
    fi

done
