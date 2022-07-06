#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */) # list of GCM models to process
# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
declare -a models=("bcc-csm1-1") # extended RCP runs
# declare -a varnames=("aht" "qaht" "saht" "vmmmc" "vqmmc" "vsmmc" "vmse" "vqse" "vsse" "vmte" "vqte" "vste" "daht" "dqaht" "dsaht" "dvmmmc" "dvqmmc" "dvsmmc" "dvmse" "dvqse" "dvsse" "dvmte" "dvqte" "dvste") # list of GCM variables that we want to process
# declare -a varnames=("gmse92500" "mse92500" "aht" "vmte")
# declare -a varnames=("mse92500")
declare -a varnames=("vmte")
sim="historical"
clim=".zonmean"
freq="Amon"
ens="r1i1p1"
yr_begin=1860
yr_end=2005
yr_span="${yr_begin}01-${yr_end}12"

cwd=$(pwd) # save current working directory

for model in ${models[@]}; do
    
    model=${model%/}
    echo ${model}

    common=${freq}_${model}_${sim}_${ens}_${yr_span}${clim}

    cd ${cwd}/${model}

    for varname in ${varnames[@]}; do
        noisy=${cwd}/${model}/${varname}_${common}.nc
        smooth=${cwd}/${model}/${varname}_${common}.shsmooth.nc

        ls $noisy

        rm $smooth
        ncl ifile=\"$noisy\" ofile=\"$smooth\" vname=\"$varname\" ${cwd}/make_shsmooth.ncl 
        cdo zonmean $smooth $smooth.new
        mv $smooth.new $smooth
    done

done
