#!/bin/sh

freq="Amon"
declare -a varnames=("ps" "ta" "hur" "hus" "tendv" "rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "adv" "ftoa" "ftoacs" "fsfc" "tas" "ts" "ra" "racs") # list of GCM variables that we want to process
# varnames=("pr" "prc" "evspsbl" "ftoa" "fsfc" "rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "tend" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "tas" "ts" "r1" "stgadv" "adv" "ra" "stf")
# varnames=("ps")

# freq="SImon"
# varnames=("siconc")

# freq="OImon"
# varnames=("sic")

# mean=".zonmean.ymonmean-30"
# mean=".zonmean.amean_70_90.ymonmean-30"
# mean=".lat_80.ymonmean-30"

# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# # declare -a models=("CanESM5/") # extended RCP runs
# ens="r1i1p1f1"

declare -a ens="r1i1p1f2" # ensemble specification 
declare -a models=("MIROC-ES2L/") # extended RCP runs
# declare -a ens="r3i1p1f2" # ensemble specification 
# declare -a models=("GISS-E2-1-G/" "GISS-E2-1-H/") # extended RCP runs
# declare -a ens="r4i1p1f2" # ensemble specification 
# declare -a models=("UKESM1-0-LL/") # extended RCP runs

yr_span="186001-229912"

# mean=".zonmean.ymonmean-30"
# mean=".ymonmean-30"
mean=""
sim="hist+ssp585"

n_yr="146" # average climatology from 1975-2005
n_myr_begin="30"

tstep_begin=$((1 + 12 * ($n_yr - $n_myr_begin)))
tstep_end=$((12 * $n_yr))

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
        cdo -seltimestep,$tstep_begin/$tstep_end ${filename}.nc ${filename}.last${n_myr_begin}.nc 
        cdo ymonmean ${filename}.last${n_myr_begin}.nc ${filename}.ymonmean-${n_myr_begin}.nc
        rm ${filename}.last${n_myr_begin}.nc

    done # varnames
done # models
