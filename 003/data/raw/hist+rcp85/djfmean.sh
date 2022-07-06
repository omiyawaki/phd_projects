#!/bin/sh

models=("bcc-csm1-1")
# declare -a models=("HadGEM2-ES" "bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# varnames=("rsdt" "rsut" "rlut" "rsds" "rsus" "rlds" "rlus" "hfls" "hfss")
# varnames=("rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" )
# varnames=("aht" "vmmmc" "vmse" "vmte")
# varnames=("gmse92500")
declare -a varnames=("gmse92500" "mse92500" "vmte" "aht" "diffv92500" "tdiffv92500") # list of GCM variables that we want to process
# declare -a varnames=("qaht" "saht" "vqmmc" "vsmmc" "vqse" "vsse" "vqte" "vste") # list of GCM variables that we want to process
sim="hist+rcp85"
freq="Amon"
ens="r1i1p1"
yr_span="186001-229912"
mean=".zonmean.shsmooth"
# mean=""

# save path to current directory
cwd=$(pwd)

for model in ${models[@]}; do

    echo ${model}

    cd ${cwd}/${model}

    for varname in ${varnames[@]}; do

        echo ${varname}

        filename="${varname}_${freq}_${model}_${sim}_${ens}_${yr_span}${mean}"

        # create DJF mean file if it doesn't exist yet
        # if [ -f "${filename}.djfmean.nc" ]; then
            # echo "DJF mean already taken, skipping..."
        # else
            cdo -seasmean -selseas,DJF ${filename}.nc ${filename}.djfmean.nc 
        # fi

    done # varnames
done # models
