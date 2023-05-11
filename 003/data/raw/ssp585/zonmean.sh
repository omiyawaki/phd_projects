#!/bin/sh

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/historical/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "GISS-E2-H/" "GISS-E2-R/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# declare -a models=("bcc-csm1-1/" "HadGEM2-ES/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# declare -a models=("CNRM-CM5/" "HadGEM2-ES/") # extended RCP runs
# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "CESM2-WACCM/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# declare -a models=("MRI-ESM2-0/") # extended RCP runs
# varnames=("rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs")
# varnames=("rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "adv" "ftoa" "ftoacs" "fsfc" "tas" "ts" "ra" "racs")
# varnames=("rlut")
# varnames=("r1" "stgadv" "adv" "ra" "stf")
# varnames=("vvmmmc" "vvmse" "vvqmmc" "vvqse" "vvsmmc" "vvsse")
# varnames=("aht" "vmmmc" "vmse" "vmte")
# varnames=("vqse" "vqmmc" "vsse" "vsmmc")
# varnames=("aht" "qaht" "saht" "vmmmc" "vqmmc" "vsmmc" "vmse" "vqse" "vsse" "vmte" "vqte" "vste")
# varnames=("daht" "dqaht" "dsaht" "dvmmmc" "dvqmmc" "dvsmmc" "dvmse" "dvqse" "dvsse" "dvmte" "dvqte" "dvste")
sim="ssp585"

freq="Amon"
# declare -a varnames=("ftoacs" "racs" "ftoa" "fsfc" "r1" "stgadv" "adv" "ra" "stf" "huss" "hurs" "ps" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "prc" "evspsbl" "rlutcs" "rldscs" "rsuscs" "rsdscs" "rsutcs") # list of GCM variables that we want to process
declare -a varnames=("r1" "stgadv" "adv" "tend") # list of GCM variables that we want to process

# freq="SImon"
# declare -a varnames=("siconc") # list of GCM variables that we want to process

yr_span="201501-229912"
mean=""

# ens="r1i1p1f1"
# # declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# declare -a models=("CESM2-WACCM/") # extended RCP runs

# declare -a ens="r1i1p1f2" # ensemble specification 
# declare -a models=("MIROC-ES2L/") # extended RCP runs

declare -a ens="r3i1p1f2" # ensemble specification 
declare -a models=("GISS-E2-1-G/" "GISS-E2-1-H/") # extended RCP runs

# declare -a ens="r4i1p1f2" # ensemble specification 
# declare -a models=("UKESM1-0-LL/") # extended RCP runs

# save path to current directory
cwd=$(pwd)

for model in ${models[@]}; do

    model=${model%/}
    echo ${model}

    cd ${cwd}/${model}

    for varname in ${varnames[@]}; do

        echo ${varname}

        grid=$(basename ${cwd}/${model}/rlut_Amon_${model}_${sim}_${ens}_*_${yr_span}.nc)
        grid=${grid#rlut_Amon_${model}_${sim}_${ens}_}
        grid=${grid%_${yr_span}.nc}

        filename="${varname}_${freq}_${model}_${sim}_${ens}_${grid}_${yr_span}${mean}"

        if [[ "$varname" == "ta" || "$varname" == "hus" || "$varname" == "hur" ]]; then
            echo "Aborting: Don't use this script for 3D data!"
            exit
        else

        # create ymonmean mean file if it doesn't exist yet
        # if [ -f "${filename}.zonmean.nc" ]; then
        #     echo "zonmean already taken, skipping..."
        # else
            cdo zonmean ${filename}.nc ${filename}.zonmean.nc
        # fi

        fi

    done # varnames
done # models
