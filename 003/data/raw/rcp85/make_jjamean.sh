#!/bin/sh

models=("MPI-ESM-LR")
# varnames=("rsdt" "rsut" "rlut" "rsds" "rsus" "rlds" "rlus" "hfls" "hfss")
varnames=("tas" "ts" "ta" "zg")
sim="rcp85"
freq="Amon"
ens="r1i1p1"
yr_span="200601-230012"

# save path to current directory
cwd=$(pwd)

for model in ${models[@]}; do

    echo ${model}

    cd ${cwd}/${model}

    for varname in ${varnames[@]}; do

        echo ${varname}

        filename="${varname}_${freq}_${model}_${sim}_${ens}_${yr_span}"

        # create JJA mean file if it doesn't exist yet
        if [ -f "${filename}.jjamean.nc" ]; then
            echo "JJA mean already taken, skipping..."
        else
            cdo -seasmean -selseas,JJA ${filename}.nc ${filename}.jjamean.nc 
        fi

    done # varnames
done # models
