#!/bin/sh

varnames=("rsdt" "rsut" "rlut" "rsds" "rsus" "rlds" "rlus" "hfls" "hfss")
# varnames=("rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs" )
# varnames=("tas" "huss")
# varnames=("mses")
sim="piControl"
freq="Amon"
yr_span="280001-284912"
ens="r1i1p1"
models=("MPI-ESM-LR/")

# save path to current directory
cwd=$(pwd)

for model in ${models[@]}; do

    model=${model%/}
    echo ${model}

    cd ${cwd}/${model}

    for varname in ${varnames[@]}; do

        echo ${varname}

        filename=$(ls ${varname}_${freq}_${model}_${sim}_${ens}_${yr_span}.nc)
        filename=${filename%.nc}
        echo $filename

        # create DJF file if it doesn't exist yet
        if [ -f "${filename}.djf.nc" ]; then
            echo "DJF already taken, skipping..."
        else
            cdo -selseas,DJF ${filename}.nc ${filename}.djf.nc 
        fi

        # create MAM file if it doesn't exist yet
        if [ -f "${filename}.mam.nc" ]; then
            echo "MAM already taken, skipping..."
        else
            cdo -selseas,MAM ${filename}.nc ${filename}.mam.nc 
        fi

        # create JJA file if it doesn't exist yet
        if [ -f "${filename}.jja.nc" ]; then
            echo "JJA already taken, skipping..."
        else
            cdo -selseas,JJA ${filename}.nc ${filename}.jja.nc 
        fi

        # create SON file if it doesn't exist yet
        if [ -f "${filename}.son.nc" ]; then
            echo "SON already taken, skipping..."
        else
            cdo -selseas,SON ${filename}.nc ${filename}.son.nc 
        fi

    done # varnames
done # models
