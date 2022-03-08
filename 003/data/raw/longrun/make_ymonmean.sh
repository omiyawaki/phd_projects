#!/bin/sh

# models=("MPIESM12_piControl")
# n_yr="1237"
# n_myr_begin="30"

models=("MPIESM12_abrupt32x")
n_yr="103"
n_myr_begin="30"
# varnames=("rsdt" "rsut" "rlut" "rsds" "rsus" "rlds" "rlus" "hfls" "hfss")
varnames=("rlutcs" "rsutcs")
freq="mon"

tstep_begin=$((1 + 12 * ($n_yr - $n_myr_begin)))
tstep_end=$((12 * $n_yr))

# save path to current directory
cwd=$(pwd)

for model in ${models[@]}; do

    echo ${model}

    cd ${cwd}/${model}

    for varname in ${varnames[@]}; do

        echo ${varname}

        filename="${varname}_${freq}_${model}_${n_yr}"

        # create ymonmean mean file if it doesn't exist yet
        if [ -f "${filename}.ymonmean-${n_myr_begin}.nc" ]; then
            echo "ymonmean-${n_myr_begin} already taken, skipping..."
        else
            cdo -seltimestep,$tstep_begin/$tstep_end ${filename}.nc ${filename}.last${n_myr_begin}.nc 
            cdo ymonmean ${filename}.last${n_myr_begin}.nc ${filename}.ymonmean-${n_myr_begin}.nc
            rm ${filename}.last${n_myr_begin}.nc
        fi

    done # varnames
done # models
