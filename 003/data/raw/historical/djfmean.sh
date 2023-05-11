#!/bin/sh

freq="Amon"
# declare -a varnames=("ftoacs" "racs" "ftoa" "fsfc" "r1" "stgadv" "adv" "ra" "stf" "huss" "hurs" "ps" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "prc" "evspsbl" "rlutcs" "rldscs" "rsuscs" "rsdscs" "rsutcs") # list of GCM variables that we want to process
# declare -a varnames=("stf" "evspsbl" "stgadv" "r1" "ps" "rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "adv" "ftoa" "ftoacs" "fsfc" "tas" "ts" "ra" "racs") # list of GCM variables that we want to process
# varnames=("tendv" "r1" "hus" "stgadv" "adv" "tend")
varnames=("r1a")

# freq="SImon"
# varnames=("siconc")

# freq="OImon"
# varnames=("sic")

# # declare -a models=("HadGEM2-ES" "bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# models=("GISS-E2-H/" "GISS-E2-R/")
# ens="r1i1p1"
# yr_span="186001-200512"

# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "CESM2-WACCM/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# # declare -a models=("MRI-ESM2-0/") # extended RCP runs
# ens="r1i1p1f1"

# declare -a models=("MIROC-ES2L/") # extended RCP runs
# ens="r1i1p1f2"

# declare -a models=("GISS-E2-1-G/" "GISS-E2-1-H/") # extended RCP runs
# ens="r3i1p1f2"

declare -a models=("UKESM1-0-LL/") # extended RCP runs
ens="r4i1p1f2"

yr_span="186001-201412"

# mean=".zonmean.ymonmean-30"
# mean=".ymonmean-30"
mean=".ymonmean-30"
# mean=".zonmean.ymonmean-30"
# mean=".zonmean.amean_70_90.ymonmean-30"
# mean=".lat_80.ymonmean-30"

sim="historical"

# save path to current directory
cwd=$(pwd)

for model in ${models[@]}; do
    model=${model%/}
    echo ${model}

    cd ${cwd}/${model}

    for varname in ${varnames[@]}; do

        echo ${varname}

        filename=$(basename ${cwd}/${model}/${varname}_${freq}_${model}_${sim}_${ens}_*${yr_span}.nc)
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
