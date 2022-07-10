#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

declare -a ens="r1i1p1f1" # ensemble specification 
declare -a models=("MRI-ESM2-0/") # extended RCP runs

sim="ssp585"
freq="Amon"

yr_spans=("201501-206412" "206501-210012" "210101-215012" "215101-220012" "220101-225012" "225101-230012")

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    cd ${cwd}/${model}

    tend0=${cwd}/${model}/tend_Amon_MRI-ESM2-0_ssp585_r1i1p1f1_gn_201501-206412.nc
    tend1=${cwd}/${model}/tend_Amon_MRI-ESM2-0_ssp585_r1i1p1f1_gn_206501-210012.nc
    tend2=${cwd}/${model}/tend_Amon_MRI-ESM2-0_ssp585_r1i1p1f1_gn_210101-215012.nc
    tend3=${cwd}/${model}/tend_Amon_MRI-ESM2-0_ssp585_r1i1p1f1_gn_215101-220012.nc
    tend4=${cwd}/${model}/tend_Amon_MRI-ESM2-0_ssp585_r1i1p1f1_gn_220101-225012.nc
    tend5=${cwd}/${model}/tend_Amon_MRI-ESM2-0_ssp585_r1i1p1f1_gn_225101-230012.nc
    full_tend=${cwd}/${model}/tend_Amon_MRI-ESM2-0_ssp585_r1i1p1f1_gn_201501-230012.nc

    cdo mergetime $tend0 $tend1 $tend2 $tend3 $tend4 $tend5 $full_tend

done
