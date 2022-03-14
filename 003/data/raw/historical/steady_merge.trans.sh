#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
declare -a models=("HadGEM2-ES") # extended RCP runs
sim="historical"
freq="Amon"
ens="r1i1p1"
yr_begin=1860
yr_end=2005
yr_span="${yr_begin}01-${yr_end}12"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}

    cd ${cwd}/${model}

    full_vmse=${cwd}/${model}/steady/vmse_*.nc
    full_merge=${cwd}/${model}/vmse_${common}.nc

    cdo -O mergetime $full_vmse $full_merge

    full_vmmmc=${cwd}/${model}/steady/vmmmc_*.nc
    full_merge=${cwd}/${model}/vmmmc_${common}.nc

    cdo -O mergetime $full_vmmmc $full_merge

done
