#!/bin/sh

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/historical/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "GISS-E2-H/" "GISS-E2-R/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# declare -a models=("bcc-csm1-1/" "HadGEM2-ES/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# declare -a models=("CNRM-CM5/" "HadGEM2-ES/") # extended RCP runs
declare -a models=("ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "CESM2-WACCM/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# declare -a models=("ACCESS-CM2/") # extended RCP runs
varnames=("rlut")
# varnames=("rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "adv" "ftoa" "ftoacs" "fsfc" "tas" "ts" "ra" "racs")
# varnames=("racs" "ftoacs")
# varnames=("r1" "stgadv" "adv" "ra" "stf")
# varnames=("vvmmmc" "vvmse" "vvqmmc" "vvqse" "vvsmmc" "vvsse")
# varnames=("aht" "vmmmc" "vmse" "vmte")
# varnames=("vqse" "vqmmc" "vsse" "vsmmc")
# varnames=("aht" "qaht" "saht" "vmmmc" "vqmmc" "vsmmc" "vmse" "vqse" "vsse" "vmte" "vqte" "vste")
# varnames=("daht" "dqaht" "dsaht" "dvmmmc" "dvqmmc" "dvsmmc" "dvmse" "dvqse" "dvsse" "dvmte" "dvqte" "dvste")
sim="historical"
freq="Amon"
ens="r1i1p1f1"
yr_span="186001-201412"
mean=""

# save path to current directory
cwd=$(pwd)

for model in ${models[@]}; do

    model=${model%/}
    echo ${model}

    cd ${cwd}/${model}

    for varname in ${varnames[@]}; do

        echo ${varname}

        grid=$(basename ${cwd}/${model}/rlut_${freq}_${model}_${sim}_${ens}_*_${yr_span}.nc)
        grid=${grid#rlut_${freq}_${model}_${sim}_${ens}_}
        grid=${grid%_${yr_span}.nc}

        filename="${varname}_${freq}_${model}_${sim}_${ens}_${grid}_${yr_span}${mean}"

        if [[ "$varname" == "ta" || "$varname" == "hus" || "$varname" == "hur" ]]; then
            echo "Aborting: Don't use this script for 3D data!"
            exit
        else

        # create ymonmean mean file if it doesn't exist yet
        if [ -f "${filename}.zonmean.nc" ]; then
            echo "zonmean already taken, skipping..."
        else
            cdo zonmean ${filename}.nc ${filename}.zonmean.nc
        fi

        fi

    done # varnames
done # models
