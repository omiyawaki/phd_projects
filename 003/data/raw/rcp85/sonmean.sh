#!/bin/sh

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# models=("IPSL-CM5A-LR/")
# models=("HadGEM2-ES/")
declare -a models=("bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "HadGEM2-ES/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# varnames=("rsdt" "rsut" "rlut" "rsds" "rsus" "rlds" "rlus" "hfls" "hfss" "pr" "prc" "evspsbl")
varnames=("tend" "r1" "stgadv" "adv" "ra" "stf" "rsdt" "rsut" "rlut" "rsds" "rsus" "rlds" "rlus" "hfls" "hfss" "pr" "prc" "evspsbl")
# varnames=("rsdt" "rsut" "rlut" "rsds" "rsus" "rlds" "rlus" "hfls" "hfss")
# varnames=("rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs")
# varnames=("divaht" "divmmc" "divse" "divte")
# varnames=("aht" "vmmmc" "vmse" "vmte")
# varnames=("ps" "tas" "ta" "zg")
# varnames=("pr" "prc")
# varnames=("eke" "veke")
# varnames=("sic")
sim="rcp85"
freq="Amon"
ens="r1i1p1"
yr_span="200601-229912"

# save path to current directory
cwd=$(pwd)

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    cd ${cwd}/${model}

    for varname in ${varnames[@]}; do

        echo ${varname}

        filename="${varname}_${freq}_${model}_${sim}_${ens}_${yr_span}"

        # create son mean file if it doesn't exist yet
        # if [ -f "${filename}.sonmean.nc" ]; then
        #     echo "son mean already taken, skipping..."
        # else
            cdo -seasmean -selseas,son ${filename}.nc ${filename}.sonmean.nc 
        # fi

    done # varnames
done # models
