#!/bin/sh

freq="Amon"
# declare -a varnames=("stf" "evspsbl" "prfrac" "prc" "pr" "tend" "stgadv" "r1" "ps" "rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "adv" "ftoa" "ftoacs" "fsfc" "tas" "ts" "ra" "racs") # list of GCM variables that we want to process
# declare -a varnames=("stf" "evspsbl" "stgadv" "r1" "ps" "rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "adv" "ftoa" "ftoacs" "fsfc" "tas" "ts" "ra" "racs") # list of GCM variables that we want to process
# declare -a varnames=("clt" "clivi" "clwvi") # list of GCM variables that we want to process
declare -a varnames=("cl" "cli" "clw") # list of GCM variables that we want to process
# declare -a varnames=("zg") # list of GCM variables that we want to process

# varnames=("pr" "prc" "evspsbl" "ftoa" "fsfc" "rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "tend" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "tas" "ts" "r1" "stgadv" "adv" "ra" "stf")
# varnames=("ps" "zg" "ta" "hus" "hur")
# varnames=("prfrac")

# freq="SImon"
# varnames=("siconc")

# freq="OImon"
# varnames=("sic")

# mean=".zonmean.ymonmean-30"
# mean=".zonmean.amean_70_90.ymonmean-30"
# mean=".lat_80.ymonmean-30"

# models=("bcc-csm1-1")
# declare -a models=("HadGEM2-ES" "bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
declare -a models=("HadGEM2-ES" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# declare -a models=("GISS-E2-H/" "GISS-E2-R/") # extended RCP runs
ens="r1i1p1"
yr_span="186001-229912"

# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# # declare -a models=("CanESM5/") # extended RCP runs
# ens="r1i1p1"
# yr_span="186001-229912"

# mean=".zonmean.ymonmean-30"
# mean=".zonmean"
mean=""
sim="hist+rcp85"

# save path to current directory
cwd=$(pwd)

for model in ${models[@]}; do
    model=${model%/}
    echo ${model}

    cd ${cwd}/${model}

    for varname in ${varnames[@]}; do

        echo ${varname}

        filename=$(basename ${cwd}/${model}/${varname}_${freq}_${model}_${sim}_${ens}_${yr_span}.nc)
        filename=${filename%.nc}
        filename="${filename}${mean}"

        # create DJF mean file if it doesn't exist yet
        # if [ -f "${filename}.djfmean.nc" ]; then
        #     echo "DJF mean already taken, skipping..."
        if [[ $mean == *"ymonmean-30"* ]]; then
            cdo -selseas,DJF ${filename}.nc ${filename}.djfsel.nc 
            cdo -timmean ${filename}.djfsel.nc ${filename}.djfmean.nc 
            rm ${filename}.djfsel.nc 
        else
            cdo -seasmean -selseas,DJF ${filename}.nc ${filename}.djfmean.nc 
        fi

    done # varnames
done # models
