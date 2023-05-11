#!/bin/sh

# models=("MPIESM12_abrupt32x")
models=("MPIESM12_abrupt2x")
# varnames=("rsdt" "rsut" "rlut" "rsds" "rsus" "rlds" "rlus" "hfls" "hfss" "tas" "rsutcs" "rlutcs" "rsuscs" "rsdscs" "rldscs" "rluscs")
varnames=("pr")
freq="mon"
n_yr="1000"

link_prefix="https://data.iac.ethz.ch/longrunmip/modeloutput/orig"

# save path to current directory
cwd=$(pwd)

for model in ${models[@]}; do

    echo ${model}

    # if the model directory doesn't exist, create it
    if [ ! -d "${cwd}/${model}" ]; then
        mkdir ${cwd}/${model}
    fi

    cd ${cwd}/${model}

    for varname in ${varnames[@]}; do

        echo ${varname}

        varname0="${varname}_${freq}"
        # if [[ "$varname" == "thetao" ]]; then
        #     varname0="${varname}_ann"
        # else
        #     varname0="${varname}_mon"
        # fi

        filename="${varname0}_${model}_${n_yr}"

        # download data if it has not been downloaded yet
        if [ -f "${filename}.nc" ]; then
            echo "File already downloaded, skipping..."
        else
            wget ${link_prefix}/${varname}/${filename}.nc --user=longrunmip --password=data4you
        fi

        # # create annual mean file if it doesn't exist yet
        # if [ -f "${filename}.yearmean.nc" ]; then
        #     echo "Annual mean already taken, skipping..."
        # else
        #     cdo yearmean ${filename}.nc ${filename}.yearmean.nc 
        # fi

    done # varnames
done # models
