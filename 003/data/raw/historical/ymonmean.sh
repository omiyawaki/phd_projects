#!/bin/sh

freq="Amon"
# declare -a varnames=("prfrac" "ftoacs" "racs" "ftoa" "fsfc" "r1" "stgadv" "adv" "ra" "stf" "huss" "hurs" "ps" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "prc" "evspsbl" "rlutcs" "rldscs" "rsuscs" "rsdscs" "rsutcs") # list of GCM variables that we want to process
# # # varnames=("rlut" "rsut" "rsdt" "rsds" "rsus" "rlus" "rlds" "hfls" "hfss" "tend" "rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" "tas" "ts" "r1" "stgadv" "adv" "ra" "stf")
# varnames=("cl" "clt" "cli" "clw" "clivi" "clwvi")
# # varnames=("pr" "prc" "evspsbl")
varnames=("r1a")
# # # varnames=("r1" "stgadv" "adv" "ra" "stf")

# freq="OImon"
# varnames=("sic")

# freq="SImon"
# varnames=("siconc")

# # declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# declare -a models=("GISS-E2-H/" "GISS-E2-R/") # extended RCP runs
# ens="r1i1p1"
# yr_end=2005
# n_yr="146"

declare -a ens="r1i1p1f1" # ensemble specification 
# declare -a models=("CESM2-WACCM/") # extended RCP runs
declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/" "CESM2-WACCM/") # extended RCP runs

# declare -a models=("MIROC-ES2L/") # extended RCP runs
# ens="r1i1p1f2"

# declare -a models=("GISS-E2-1-G/" "GISS-E2-1-H/") # extended RCP runs
# ens="r3i1p1f2"

# declare -a models=("UKESM1-0-LL/") # extended RCP runs
# ens="r4i1p1f2"

yr_end=2014
n_yr="155" # average climatology from 1984-2014
# n_yr="146" # average climatology from 1975-2005

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
