#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

ens="r1i1p1f1"
declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "CESM2-WACCM/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs

# declare -a ens="r1i1p1f2" # ensemble specification 
# declare -a models=("MIROC-ES2L/") # extended RCP runs
# declare -a ens="r3i1p1f2" # ensemble specification 
# declare -a models=("GISS-E2-1-G/" "GISS-E2-1-H/") # extended RCP runs
# declare -a ens="r4i1p1f2" # ensemble specification 
# declare -a models=("UKESM1-0-LL/") # extended RCP runs

sim="hist+ssp585"
clim=".djfmean"
# clim=".jjamean"
freq="Amon"
yr_begin=1860
yr_end=2299
yr_span="${yr_begin}01-${yr_end}12"

ref_sim="historical"
ref_clim=".ymonmean-30"
ref_yr_begin=1860
ref_yr_end=2014
ref_yr_span="${ref_yr_begin}01-${ref_yr_end}12"

cwd=$(pwd) # save current working directory
rfd=${cwd}/../${ref_sim} # reference data directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=$(basename ${cwd}/${model}/rlut_${freq}_${model}_${sim}_${ens}_*${yr_span}${clim}.nc)
    common=${common#rlut_}
    common=${common%.nc}

    ref_common=$(basename ${cwd}/${model}/rlut_${freq}_${model}_${ref_sim}_${ens}_*${ref_yr_span}${ref_clim}${clim}.nc)
    ref_common=${ref_common#rlut_}
    ref_common=${ref_common%.nc}


    cd ${cwd}/${model}

    if ls $cwd/${model}/vmmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_r1a=${cwd}/${model}/r1a_${common}.nc
        full_ra=${cwd}/${model}/ra_${common}.nc
        full_stf=${cwd}/${model}/stf_${common}.nc
        full_dyna=${cwd}/${model}/adv_${common}.nc
        ref_r1a=${rfd}/${model}/r1a_${ref_common}.nc
        ref_ra=${rfd}/${model}/ra_${ref_common}.nc
        ref_stf=${rfd}/${model}/stf_${ref_common}.nc
        ref_dyna=${rfd}/${model}/adv_${ref_common}.nc
        diff_r1a=${cwd}/${model}/dr1a_${common}.nc
        diff_ra=${cwd}/${model}/dcraa_${common}.nc
        diff_stf=${cwd}/${model}/dcstfa_${common}.nc
        diff_dyna=${cwd}/${model}/dcdyna_${common}.nc
        diff_res=${cwd}/${model}/dcresa_${common}.nc

        python ${cwd}/make_dr1_decomp.py ${full_r1a} ${full_ra} ${full_dyna} ${ref_r1a} ${ref_ra} ${ref_dyna} ${diff_r1a} ${diff_ra} ${diff_dyna} ${diff_res} ${full_stf} ${ref_stf} ${diff_stf}


    fi

done
