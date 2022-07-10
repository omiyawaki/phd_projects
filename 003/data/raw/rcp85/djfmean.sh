#!/bin/sh

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# models=("IPSL-CM5A-LR/")
# models=("MPI-ESM-LR/")
# declare -a models=("CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# declare -a models=("bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "GISS-E2-H/" "GISS-E2-R/" "HadGEM2-ES/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# declare -a models=("bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "HadGEM2-ES/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
declare -a models=("bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "HadGEM2-ES/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# declare -a models=("CCSM4/") # extended RCP runs
# varnames=("rsdt" "rsut" "rlut" "rsds" "rsus" "rlds" "rlus" "hfls" "hfss")
# varnames=("divaht" "divte")
# varnames=("rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs")
# varnames=("divaht" "divmmc" "divse" "divte")
# varnames=("r1" "stgadv" "adv" "ra" "stf")
# varnames=("vmte_sm")
# varnames=("gmse92500")
varnames=("ftoa" "fsfc")
# varnames=("vvmmmc" "vvmse" "vvqmmc" "vvqse" "vvsmmc" "vvsse")
# varnames=("ta" "hur" "hus" "tendv")
# varnames=("tas" "rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "tend" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "r1" "stgadv" "adv" "ra" "stf" "ftoa" "ftoacs" "fsfc" "ts" "racs")
sim="rcp85"
freq="Amon"
ens="r1i1p1"
yr_span="200601-229912"
# yr_span="200601-230012"
# mean=".zonmean"
# mean=".lat_80"
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

        # create DJF mean file if it doesn't exist yet
        # if [ -f "${filename}.djfmean.nc" ]; then
        #     echo "DJF mean already taken, skipping..."
        # else
            cdo -seasmean -selseas,DJF ${filename}.nc ${filename}.djfmean.nc 
        # fi

    done # varnames
done # models
