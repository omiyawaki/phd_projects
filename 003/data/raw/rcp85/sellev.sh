#!/bin/sh

varnames=("mse")
sim="rcp85"
freq="Amon"
yr_span="200601-229912"
lev="92500" # level to extract in Pa
ens="r1i1p1"
mean=""
declare -a models=("HadGEM2-ES/") # extended RCP runs
# declare -a models=("HadGEM2-ES/" "bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs

# save path to current directory
cwd=$(pwd)

for model in ${models[@]}; do

    model=${model%/}
    echo ${model}

    cd ${cwd}/${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}${mean}

    for varname in ${varnames[@]}; do

        # if ls $cwd/${model}/${varname}${lev}_${common}.nc 1> /dev/null 2>&1; then # check if data is already there
        #     echo "sellev was already created. Skipping..."
        # else
            full_ps=${cwd}/${model}/ps_${common}.nc
            full_in=${cwd}/${model}/${varname}_${common}.nc
            full_out=${cwd}/${model}/${varname}${lev}_${common}.nc

            python ${cwd}/make_sellev.py ${full_in} ${full_out} ${full_ps} ${varname} ${lev}
        # fi


    done # varnames
done # models
