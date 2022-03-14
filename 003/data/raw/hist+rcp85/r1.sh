#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# declare -a models=("bcc-csm1-1") # extended RCP runs
sim="hist+rcp85"
clim=".djfmean"
freq="Amon"
ens="r1i1p1"
yr_begin=1860
yr_end=2299
yr_span="${yr_begin}01-${yr_end}12"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}${clim}

    cd ${cwd}/${model}

    if ls $cwd/${model}/vmmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_rsdt=${rfd}/${model}/rsdt_${ref_common}.nc
        full_rsut=${rfd}/${model}/rsut_${ref_common}.nc
        full_rlut=${rfd}/${model}/rlut_${ref_common}.nc
        full_rsds=${cwd}/${model}/rsds_${common}.nc
        full_rsus=${cwd}/${model}/rsus_${common}.nc
        full_rlus=${cwd}/${model}/rlus_${common}.nc
        full_rlds=${cwd}/${model}/rlds_${common}.nc
        full_hfls=${cwd}/${model}/hfls_${common}.nc
        full_hfss=${cwd}/${model}/hfss_${common}.nc
        full_r1=${cwd}/${model}/r1_${common}.nc
        full_ra=${cwd}/${model}/ra_${common}.nc
        full_stf=${cwd}/${model}/stf_${common}.nc

        python ${cwd}/make_decomp.py ${full_rsdt} ${ref_ra} ${ref_stf} ${diff_r1} ${diff_ra} ${diff_stf} ${diff_res} ${full_r1} ${full_ra} ${full_stf}


    fi

done
