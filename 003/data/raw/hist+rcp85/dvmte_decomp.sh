#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
declare -a models=("bcc-csm1-1") # extended RCP runs

lev=92500
llo=70
clat=80
lhi=90

sim="hist+rcp85"
clim=".zonmean.shsmooth.djfmean"
freq="Amon"
ens="r1i1p1"
yr_begin=1860
yr_end=2299
yr_span="${yr_begin}01-${yr_end}12"

ref_sim="historical"
ref_clim="zonmean.shsmooth.ymonmean-30.djfmean"
ref_yr_begin=1860
ref_yr_end=2005
ref_yr_span="${ref_yr_begin}01-${ref_yr_end}12"

cwd=$(pwd) # save current working directory
rfd=${cwd}/../${ref_sim} # reference data directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}
    ref_common=${freq}_${model}_${ref_sim}_${ens}_${ref_yr_span}

    cd ${cwd}/${model}

    if ls $cwd/${model}/dvmte_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "dvmte decomposition was already performed. Skipping..."
    else
        full_vmte=${cwd}/${model}/vmte_sm_${common}${clim}.lat_${clat}.nc
        full_dvp=${cwd}/${model}/dvp${lev}_sm_${common}${clim}.lat_${clat}.nc
        full_gmse=${cwd}/${model}/gmse${lev}_sm_${common}.zonmean${clim}.amean_${llo}_${lhi}.nc
        ref_vmte=${rfd}/${model}/vmte_${ref_common}.lat_${clat}${ref_clim}.nc
        ref_dvp=${rfd}/${model}/dvp${lev}_${ref_common}.lat_${clat}${ref_clim}.nc
        ref_gmse=${rfd}/${model}/gmse${lev}_${ref_common}.zonmean.amean_${llo}_${lhi}${ref_clim}.nc
        diff_vmte=${cwd}/${model}/dvmte_sm_${common}${clim}.lat_${clat}.nc
        diff_dvp=${cwd}/${model}/ddvp${lev}_sm_${common}${clim}.lat_${clat}.nc
        diff_gmse=${cwd}/${model}/dgmse${lev}_sm_${common}.zonmean${clim}.amean_${llo}_${lhi}.nc
        diff_res=${cwd}/${model}/ddvpres${lev}_sm_${common}${clim}.lat_${clat}.nc

        python ${cwd}/make_dvmte_decomp.py ${full_vmte} ${full_dvp} ${full_gmse} ${ref_vmte} ${ref_dvp} ${ref_gmse} ${diff_vmte} ${diff_dvp} ${diff_gmse} ${diff_res} ${clat}


    fi

done
