#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

########## 30 m w/ ice
# declare -a models=("rp000184") # extended RCP runs
# declare -a ref_model=("rp000130") # extended RCP runs

########## 30 m w/o ice
# declare -a models=("rp000185") # extended RCP runs
# declare -a ref_model=("rp000131") # extended RCP runs

########## 40 m w/o ice QFLUX
# declare -a models=("rp000191") # extended RCP runs
# declare -a ref_model=("rp000190") # extended RCP runs

########## 40 m w/o ice QFLUX
# declare -a models=("rp000191b") # extended RCP runs
# declare -a ref_model=("rp000190b") # extended RCP runs

########## 40 m w/ ice
declare -a models=("rp000134") # extended RCP runs
ref_model=("rp000134")

# ########## 40 m w/o ice
# declare -a models=("rp000189") # extended RCP runs
# declare -a ref_model=("rp000135") # extended RCP runs

clim=".djfmean"
sim="echam"
yr_span="0020_0039"

ref_clim=".ymonmean-20"
ref_yr_span="0020_0039"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${model}_${yr_span}${clim}
    ref_common=${ref_model}_${ref_yr_span}${ref_clim}${clim}

    cd ${cwd}/${model}

    if ls $cwd/${model}/vmmmc_${common} 1> /dev/null 2>&1; then # check if data is already there
        echo "mmc and se was already created. Skipping..."
    else
        full_r1=${cwd}/${model}/r1_${common}.nc
        full_ra=${cwd}/${model}/ra_${common}.nc
        full_stf=${cwd}/${model}/stf_${common}.nc
        full_dyn=${cwd}/${model}/stgadv_${common}.nc
        ref_r1=${cwd}/${ref_model}/r1_${ref_common}.nc
        ref_ra=${cwd}/${ref_model}/ra_${ref_common}.nc
        ref_stf=${cwd}/${ref_model}/stf_${ref_common}.nc
        ref_dyn=${cwd}/${ref_model}/stgadv_${ref_common}.nc
        diff_r1=${cwd}/${model}/dr1_${common}.nc
        diff_ra=${cwd}/${model}/dcra_${common}.nc
        diff_stf=${cwd}/${model}/dcstf_${common}.nc
        diff_dyn=${cwd}/${model}/dcdyn_${common}.nc
        diff_res=${cwd}/${model}/dcres_${common}.nc

        python ${cwd}/make_decomp.py ${full_r1} ${full_ra} ${full_dyn} ${ref_r1} ${ref_ra} ${ref_dyn} ${diff_r1} ${diff_ra} ${diff_dyn} ${diff_res} ${full_stf} ${ref_stf} ${diff_stf}


    fi

done
