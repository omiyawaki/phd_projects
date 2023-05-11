#!/bin/sh

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/historical/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "GISS-E2-H/" "GISS-E2-R/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# declare -a models=("HadGEM2-ES/" "bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
declare -a models=("GISS-E2-H/" "GISS-E2-R/") # extended RCP runs

# declare -a varnames=("ftoacs" "racs" "ftoa" "fsfc" "r1" "stgadv" "adv" "ra" "stf" "huss" "hurs" "ps" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "prc" "evspsbl" "rlutcs" "rldscs" "rsuscs" "rsdscs" "rsutcs") # list of GCM variables that we want to process
varnames=("r1" "stgadv" "adv")
# varnames=("tas" "rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "tend" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "r1" "stgadv" "adv" "ra" "stf" "ftoa" "ftoacs" "fsfc" "ts" "racs")
# varnames=("tas")
# varnames=("r1" "stgadv" "adv" "ra" "stf")
# varnames=("vvmmmc" "vvmse" "vvqmmc" "vvqse" "vvsmmc" "vvsse")
# varnames=("aht" "vmmmc" "vmse" "vmte")
# varnames=("vqse" "vqmmc" "vsse" "vsmmc")
# varnames=("aht" "qaht" "saht" "vmmmc" "vqmmc" "vsmmc" "vmse" "vqse" "vsse" "vmte" "vqte" "vste")
# varnames=("daht" "dqaht" "dsaht" "dvmmmc" "dvqmmc" "dvsmmc" "dvmse" "dvqse" "dvsse" "dvmte" "dvqte" "dvste")
sim="rcp85"
freq="Amon"
ens="r1i1p1"
yr_span="200601-229912"
# mean=".ymonmean-30"
mean=""

# save path to current directory
cwd=$(pwd)

for model in ${models[@]}; do

    model=${model%/}
    echo ${model}

    cd ${cwd}/${model}

    for varname in ${varnames[@]}; do

        echo ${varname}

        filename="${varname}_${freq}_${model}_${sim}_${ens}_${yr_span}${mean}"

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
