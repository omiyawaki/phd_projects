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

    grid=$(basename $cwd/../historical/${model}/rlut_${freq}_${model}_historical_${ens}_*_186001-201412.nc)
    echo $grid
    grid=${grid#rlut_${freq}_${model}_historical_${ens}_}
    grid=${grid%_186001-201412.nc}

    common=${freq}_${model}_${sim}_${ens}_${grid}_${yr_span}${clim}
    ref_common=${freq}_${model}_${ref_sim}_${ens}_${grid}_${ref_yr_span}${ref_clim}${clim}

    cd ${cwd}/${model}

    if ls $cwd/${model}/vmmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_r1=${cwd}/${model}/r1_${common}.nc
        full_ra=${cwd}/${model}/ra_${common}.nc
        full_dyn=${cwd}/${model}/stgadv_${common}.nc
        ref_r1=${rfd}/${model}/r1_${ref_common}.nc
        ref_ra=${rfd}/${model}/ra_${ref_common}.nc
        ref_dyn=${rfd}/${model}/stgadv_${ref_common}.nc
        diff_r1=${cwd}/${model}/dr1_${common}.nc
        diff_ra=${cwd}/${model}/dcra_${common}.nc
        diff_dyn=${cwd}/${model}/dcdyn_${common}.nc
        diff_res=${cwd}/${model}/dcres_${common}.nc

        python ${cwd}/make_decomp.py ${full_r1} ${full_ra} ${full_dyn} ${ref_r1} ${ref_ra} ${ref_dyn} ${diff_r1} ${diff_ra} ${diff_dyn} ${diff_res}


    fi

done
