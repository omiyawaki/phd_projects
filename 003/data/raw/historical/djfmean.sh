#!/bin/sh

models=("MPI-ESM-LR")
varnames=("rsdt" "rsut" "rlut" "rsds" "rsus" "rlds" "rlus" "hfls" "hfss" "pr" "prc" "evspsbl")
# varnames=("rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" )
# varnames=("pr" "prc")
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

        filename="${varname}_${freq}_${model}_${sim}_${ens}_${yr_span}"

        # create DJF mean file if it doesn't exist yet
        if [ -f "${filename}.djfmean.nc" ]; then
            echo "DJF mean already taken, skipping..."
        else
            cdo -seasmean -selseas,DJF ${filename}.nc ${filename}.djfmean.nc 
        fi

    done # varnames
done # models
