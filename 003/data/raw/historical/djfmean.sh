#!/bin/sh

models=("CCSM4")
# declare -a models=("HadGEM2-ES" "bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# varnames=("rsdt" "rsut" "rlut" "rsds" "rsus" "rlds" "rlus" "hfls" "hfss" "pr" "prc" "evspsbl")
# varnames=("rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" )
# varnames=("hur" "hus" "ta")
# varnames=("ts")
varnames=("hus" "ta")
# varnames=("rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "tend" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs")
# varnames=("gmse92500")
# varnames=("vmte" "dvp92500")
# varnames=("aht" "vmmmc" "vmse" "vmte" "qaht" "vqmmc" "vqse" "vqte" "saht" "vsmmc" "vsse" "vste")
# varnames=("r1" "stgadv" "adv" "ra" "stf")
# varnames=("vvmmmc" "vvmse" "vvqmmc" "vvqse" "vvsmmc" "vvsse")
# mean=".zonmean.ymonmean-30"
# mean=".zonmean.amean_70_90.ymonmean-30"
# mean=".lat_80.ymonmean-30"
mean=".zonmean"
# mean=".zonmean.ymonmean-30"
sim="historical"
freq="Amon"
ens="r1i1p1"
yr_span="186001-200512"

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
