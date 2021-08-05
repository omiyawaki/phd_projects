#!/usr/bin/env bash
set -euo pipefail

# declare -a vars_gcm=("zg" "ta" "hur" "ps" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "prc" "evspsbl") # list of GCM variables that we want to process
# declare -a vars_gcm=("rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss") # list of GCM variables that we want to process
# declare -a vars_gcm=("rlutcs" "rsutcs" "rldscs" "rsdscs" "rsuscs") # list of GCM variables that we want to process
# declare -a vars_gcm=("ps" "tas" "ta" "zg" "ts" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "evspsbl") # list of GCM variables that we want to process
# declare -a vars_gcm=("ps" "tas" "ta" "zg" "ts" "hus" "huss") # list of GCM variables that we want to process
# declare -a vars_gcm=("vas") # list of GCM variables that we want to process
declare -a vars_gcm=("pr" "prc" "evspsbl") # list of GCM variables that we want to process
declare -a realm=("atmos")
declare -a clim="historical" # climate name
declare -a freq="mon" # data output frequency (e.g. fx for fixed, mon for monthly, day for daily)
declare -a ens="r1i1p1" # ensemble specification 
# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/historical/ && ls -d */) # list of GCM models to process
declare -a models=("MPI-ESM-LR/") # list of GCM models to process
declare -a skip_models="FGOALS-s2/"
declare -a skip_files=("_eady.nc")
# declare -a skip_files=("185001-200512.nc _eady.nc")
# declare -a skip_files=("185001-201212.nc 185001-200512.nc _eady.nc")

out_yr_begin=1860
out_mn_begin=01
out_yr_end=2005
out_mn_end=12

cwd=$(pwd) # save current working directory

# for dirs in */; do # loop through all the models
for dirs in ${models[@]}; do # loop through models
    case $skip_models in *"$dirs"*)
        :
        ;; 
    *) 
        echo $dirs
        cd $dirs

        for vars in ${vars_gcm[@]}; do
            echo $vars
            # if ls $cwd/${dirs}${vars}_*historical*${out_yr_begin}${out_mn_begin}-${out_yr_end}${out_mn_end}*.nc 1> /dev/null 2>&1; then # check if data is already there
            #     echo "${vars} was already converted. Skipping..."
            # else

                pattern="${vars}_*${dirs%/}*200512.nc"
                files=( $pattern )

                if [[ ${files[@]} != *.nc ]]; then # check if this file exists
                    echo "File of type $pattern does not exist. Please download the file and place it in the corresponding directory."
                elif [ ${#files[@]} -eq 1 ]; then # check if there are multiple files of this variable
                    common="${files%_*-*}_" # common part of files name that doesn't contain the years
                    yr=${files#"${files%_*-*}_"} # extract just the year part of name containing last year
                    yr_start=${yr: 0:4}
                    yr_end=${yr: 7:4} # extract the last year and month
                    echo $yr_start
                    echo $yr_end
                    sel="${common}${out_yr_begin}${out_mn_begin}-${out_yr_end}${out_mn_end}" # write files name with last 5 years
                    merge=${files[@]%.nc}

                    if [ ! $yr_start -eq ${out_yr_begin} ]; then
                        selmerge="${common}${out_yr_begin}${out_mn_begin}-${out_yr_end}${out_mn_end}" # write file name with merged time
                        cdo -O seldate,${out_yr_begin}-${out_mn_begin}-01,${out_yr_end}-${out_mn_end}-31 $cwd/$dirs$merge.nc $cwd/$dirs$selmerge.nc
                        cdo -O yearmean $cwd/$dirs$selmerge.nc $cwd/$dirs$selmerge.yearmean.nc # combine multiple files into one
                    fi

                else
                    echo "There are multiple files ending with year 2005. Make sure there are no duplicates of this data."
                fi

        done
        
        ;;
    esac 

    cd $cwd # go back to base directory
    
done

cd $cwd
