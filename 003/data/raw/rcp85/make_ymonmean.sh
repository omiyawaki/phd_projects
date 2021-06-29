#!/bin/sh

models=("MPI-ESM-LR")
# varnames=("rsdt" "rsut" "rlut" "rsds" "rsus" "rlds" "rlus" "hfls" "hfss")
varnames=("pr")
sim="rcp85"
freq="Amon"
ens="r1i1p1"
yr_span="200601-230012"

n_yr="295"
n_myr_begin="30"

tstep_begin=$((1 + 12 * ($n_yr - $n_myr_begin)))
tstep_end=$((12 * $n_yr))

# save path to current directory
cwd=$(pwd)

for model in ${models[@]}; do

    echo ${model}

    cd ${cwd}/${model}

    for varname in ${varnames[@]}; do

        echo ${varname}

        filename="${varname}_${freq}_${model}_${sim}_${ens}_${yr_span}"

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
