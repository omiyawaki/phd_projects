#!/bin/sh

# freq="Amon"
# # varnames=("rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "tend" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "tas" "ts" "r1" "stgadv" "adv" "ra" "stf")
# # varnames=("hur" "hus" "ta")
# varnames=("pr" "prc" "evspsbl")
# # varnames=("prfrac")
# # varnames=("r1" "stgadv" "adv" "ra" "stf")

freq="SImon"
varnames=("siconc")

# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# # declare -a models=("HadGEM2-ES") # extended RCP runs
# ens="r1i1p1"
# yr_end=2005
# n_yr="146"

declare -a ens="r1i1p1f1" # ensemble specification 
declare -a models=("CanESM5/") # extended RCP runs
# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
yr_end=2014
# n_yr="155" # average climatology from 1984-2014
n_yr="146" # average climatology from 1975-2005

sim="historical"
yr_begin=1860
yr_span="${yr_begin}01-${yr_end}12"

# mean=".zonmean.amean_70_90"
# mean=".lat_80"
# mean=".zonmean"
mean=""

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
        echo $filename

        # create ymonmean mean file if it doesn't exist yet
        # if [ -f "${filename}.ymonmean-${n_myr_begin}.nc" ]; then
        #     echo "ymonmean-${n_myr_begin} already taken, skipping..."
        # else
            cdo -seltimestep,$tstep_begin/$tstep_end ${filename}.nc ${filename}.last${n_myr_begin}.nc 
            cdo ymonmean ${filename}.last${n_myr_begin}.nc ${filename}.ymonmean-${n_myr_begin}.nc
            rm ${filename}.last${n_myr_begin}.nc
        # fi

    done # varnames
done # models
