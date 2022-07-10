#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

declare -a ens="r1i1p1f1" # ensemble specification 
declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "CESM2-WACCM/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs

sim="ssp585"
freq="Amon"
yr_begin=2015
yr_end=2300
yr_span="${yr_begin}01-${yr_end}12"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=$(basename ${cwd}/${model}/rlut_${freq}_${model}_${sim}_${ens}_*${yr_span}.nc)
    common=${common#rlut_}
    common=${common%.nc}

    cd ${cwd}/${model}

    pr=${cwd}/${model}/pr_${common}.nc
    prc=${cwd}/${model}/prc_${common}.nc
    prl=${cwd}/${model}/prl_${common}.nc
    prfrac=${cwd}/${model}/prfrac_${common}.nc

    python ${cwd}/make_prfrac.py ${pr} ${prc} ${prl} ${prfrac}


done
