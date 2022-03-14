#!/bin/sh

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/historical/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "GISS-E2-H/" "GISS-E2-R/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
declare -a models=("HadGEM2-ES/" "bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# declare -a models=("bcc-csm1-1/") # extended RCP runs
# varnames=("rsutcs" "rlutcs" "rsdscs" "rsuscs" "rldscs")
# varnames=("t850")
varnames=("r1" "stgadv" "adv" "ra" "stf")
# varnames=("aht" "vmmmc" "vmse" "vmte")
# varnames=("vqse" "vqmmc" "vsse" "vsmmc")
# varnames=("aht" "qaht" "saht" "vmmmc" "vqmmc" "vsmmc" "vmse" "vqse" "vsse" "vmte" "vqte" "vste")
# varnames=("daht" "dqaht" "dsaht" "dvmmmc" "dvqmmc" "dvsmmc" "dvmse" "dvqse" "dvsse" "dvmte" "dvqte" "dvste")
sim="historical"
freq="Amon"
ens="r1i1p1"
yr_span="186001-200512"

n_yr="146"
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
