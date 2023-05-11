#!/bin/sh

freq="Amon"
# declare -a varnames=("stf" "evspsbl" "prfrac" "prc" "pr" "tend" "stgadv" "r1" "ps" "rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "adv" "ftoa" "ftoacs" "fsfc" "tas" "ts" "ra" "racs") # list of GCM variables that we want to process
# declare -a varnames=("stf" "evspsbl" "stgadv" "r1" "ps" "rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "adv" "ftoa" "ftoacs" "fsfc" "tas" "ts" "ra" "racs") # list of GCM variables that we want to process
# declare -a varnames=("r1" "stgadv" "adv" "tend" "tendv") # list of GCM variables that we want to process
# declare -a varnames=("cl" "clt" "cli" "clw" "clivi" "clwvi") # list of GCM variables that we want to process
declare -a varnames=("r1a") # list of GCM variables that we want to process
# declare -a varnames=("rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "ra" "racs") # list of GCM variables that we want to process
# declare -a varnames=("ta" "hur" "hus" "tendv" "rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "adv" "ftoa" "ftoacs" "fsfc" "tas" "ts" "ra" "racs") # list of GCM variables that we want to process
# declare -a varnames=("rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "adv" "ftoa" "ftoacs" "fsfc" "tas" "ts" "ra" "racs") # list of GCM variables that we want to process
# varnames=("evspsbl")

# freq="SImon"
# varnames=("siconc")

# freq="OImon"
# varnames=("sic")

# mean=".zonmean.ymonmean-30"
# mean=".zonmean.amean_70_90.ymonmean-30"
# mean=".lat_80.ymonmean-30"

# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "CESM2-WACCM/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# # declare -a models=("CanESM5/" "CESM2-WACCM/") # extended RCP runs
# ens="r1i1p1f1"

# declare -a ens="r1i1p1f2" # ensemble specification 
# declare -a models=("MIROC-ES2L/") # extended RCP runs
# declare -a ens="r3i1p1f2" # ensemble specification 
# declare -a models=("GISS-E2-1-G/" "GISS-E2-1-H/") # extended RCP runs
declare -a ens="r4i1p1f2" # ensemble specification 
declare -a models=("UKESM1-0-LL/") # extended RCP runs

yr_span="186001-229912"

# mean=".zonmean.ymonmean-30"
mean=""
# mean=".zonmean"
sim="hist+ssp585"

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

        filename=$(basename ${cwd}/${model}/${varname}_${freq}_${model}_${sim}_${ens}_${grid}_${yr_span}.nc)
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
