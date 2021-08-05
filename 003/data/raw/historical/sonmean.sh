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

        # create son mean file if it doesn't exist yet
        if [ -f "${filename}.sonmean.nc" ]; then
            echo "son mean already taken, skipping..."
        else
            cdo -seasmean -selseas,son ${filename}.nc ${filename}.sonmean.nc 
        fi

    done # varnames
done # models
