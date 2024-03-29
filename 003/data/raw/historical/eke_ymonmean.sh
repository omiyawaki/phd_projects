#!/bin/sh

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/historical/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "GISS-E2-H/" "GISS-E2-R/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# declare -a models=("IPSL-CM5A-LR/") # extended RCP runs
declare -a models=("MPI-ESM-LR/")
varnames=("eke" "veke")
# varnames=("ta")
# varnames=("ps" "tas" "ta" "zg" "hus")
sim="historical"
freq="Amon"
ens="r1i1p1"
yr_span="195001-200512"

n_yr="56"
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
