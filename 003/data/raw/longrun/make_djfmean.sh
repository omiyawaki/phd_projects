#!/bin/sh

models=("MPIESM12_abrupt32x")
# varnames=("rsdt" "rsut" "rlut" "rsds" "rsus" "rlds" "rlus" "hfls" "hfss")
varnames=("sic")
freq="mon"
n_yr="103"

# save path to current directory
cwd=$(pwd)

for model in ${models[@]}; do

    echo ${model}

    cd ${cwd}/${model}

    for varname in ${varnames[@]}; do

        echo ${varname}

        filename="${varname}_${freq}_${model}_${n_yr}"

        # create DJF mean file if it doesn't exist yet
        if [ -f "${filename}.djfmean.nc" ]; then
            echo "DJF mean already taken, skipping..."
        else
            cdo -seasmean -selseas,DJF ${filename}.nc ${filename}.djfmean.nc 
        fi

    done # varnames
done # models
