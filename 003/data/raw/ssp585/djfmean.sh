#!/bin/sh

freq="Amon"
# varnames=("ftoa" "fsfc" "rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "tend" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "tas" "ts" "r1" "stgadv" "adv" "ra" "stf")
# varnames=("rsdt" "rsut" "rlut" "rsds" "rsus" "rlds" "rlus" "hfls" "hfss" "pr" "prc" "evspsbl")
# varnames=("rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" )
# varnames=("ta" "zg" "ps")
varnames=("ftoacs" "racs" "tendv" "ta" "hur" "hus")

# freq="SImon"
# varnames=("siconc")

# mean=".zonmean.ymonmean-30"
# mean=".zonmean.amean_70_90.ymonmean-30"
# mean=".lat_80.ymonmean-30"

declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/") # extended RCP runs
# declare -a models=("MRI-ESM2-0/") # extended RCP runs
ens="r1i1p1f1"
yr_span="201501-230012"

mean=".zonmean"
# mean=".zonmean.ymonmean-30"
# mean=""
sim="ssp585"

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
