#!/bin/sh

# models=("MPI-ESM-LR")
declare -a models=("bcc-csm1-1" "CCSM4" "CNRM-CM5" "CSIRO-Mk3-6-0" "GISS-E2-H" "GISS-E2-R" "HadGEM2-ES" "IPSL-CM5A-LR" "MPI-ESM-LR") # extended RCP runs
# varnames=("rsdt" "rsut" "rlut" "rsds" "rsus" "rlds" "rlus" "hfls" "hfss" "pr" "prc" "sic")
# varnames=("pr" "prc")
varnames=("hur")
sim="rcp85"
freq="Amon"
ens="r1i1p1"

# yr_span="200601-230012"
# n_yr="295"

# yr_span="200601-229912"
# n_yr="294"

# yr_span="200601-224912"
# n_yr="244"

# yr_span="200601-219912"
# n_yr="194"

# yr_span="200601-214912"
# n_yr="144"

# yr_span="200601-210012"
# n_yr="95"

# yr_span="200601-209912"
# n_yr="94"

# yr_span="200601-204912"
# n_yr="44"

# yr_span_list=("200601-229912" "200601-224912" "200601-219912" "200601-214912" "200601-209912" "200601-204912")
# n_yr_list=("294" "244" "194" "144" "94" "44")

yr_span_list=("200601-229912" "200601-224912")
n_yr_list=("294" "244")

n_myr_begin="30"

counter=1

cwd=$(pwd)

for yr_span in ${yr_span_list[@]}; do
    n_yr=${n_yr_list[$counter]}
    echo $counter
    echo $n_yr

    tstep_begin=$((1 + 12 * ($n_yr - $n_myr_begin)))
    tstep_end=$((12 * $n_yr))

    # save path to current directory

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

    counter=$((counter+1))
done
