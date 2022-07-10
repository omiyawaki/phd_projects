#!/usr/bin/env bash
set -euo pipefail

# declare -a vars_gcm=("zg" "ta" "ps" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "prc" "evspsbl" "rlutcs" "rldscs" "rsuscs" "rsdscs" "rsutcs") # list of GCM variables that we want to process
declare -a vars_gcm=("rlut") # list of GCM variables that we want to process
declare -a clim="ssp585" # climate name
declare -a freq="mon" # data output frequency (e.g. fx for fixed, mon for monthly, day for daily)
declare -a ens="r1i1p1f1" # ensemble specification 
# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
# declare -a models=("ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/") # extended RCP runs
declare -a models=("ACCESS-CM2/") # extended RCP runs
declare -a skip_models=("")
declare -a skip_files=("_eady.nc")

in_yr_out=2300

out_yr_begin=2015
out_mn_begin=01
# out_yr_end=2299
out_mn_end=12

# declare -a yrend_list=("2049" "2099" "2149" "2199" "2249")
declare -a yrend_list=("2299")

cwd=$(pwd) # save current working directory

for out_yr_end in ${yrend_list[@]}; do
    echo $out_yr_end

    # for dirs in */; do # loop through all the models
    for dirs in ${models[@]}; do # loop through models
        case $skip_models in *"$dirs"*)
            :
            ;; 
        *) 
            echo $dirs
            cd $dirs
            model=${dirs%/}

            for vars in ${vars_gcm[@]}; do
                echo $vars
                # if ls $cwd/${dirs}${vars}_*historical*${out_yr_begin}${out_mn_begin}-${out_yr_end}${out_mn_end}*.nc 1> /dev/null 2>&1; then # check if data is already there
                #     echo "${vars} was already converted. Skipping..."
                # else

                    # pattern="${vars}_*${dirs%/}*12.nc"
                    if [[ "$model" == "HadGEM2-ES" ]]; then
                        pattern="${vars}_*${dirs%/}*200601-229912.nc" # for HadGEM2-ES
                    elif [[ "$model" == "FGOALS-g2" ]]; then
                        pattern="${vars}_*${dirs%/}*210112.nc" # for FGOALS-g2
                    elif [[ "$model" == "GFDL-ESM2M" && ( ("$vars" == "ps") || ("$vars" == "tas" ) ) ]]; then
                        pattern="${vars}_*${dirs%/}*200601-220012.nc" # for GFDL-ESM2M
                    else
                        pattern="${vars}_*${dirs%/}*${in_yr_out}12.nc"
                    fi
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

                        if [ ! $yr_end -eq ${out_yr_end} ]; then
                            selmerge="${common}${out_yr_begin}${out_mn_begin}-${out_yr_end}${out_mn_end}" # write file name with merged time
                            cdo -O seldate,${out_yr_begin}-${out_mn_begin}-01,${out_yr_end}-${out_mn_end}-31 $cwd/$dirs$merge.nc $cwd/$dirs$selmerge.nc
                            cdo -O yearmean $cwd/$dirs$selmerge.nc $cwd/$dirs$selmerge.yearmean.nc # combine multiple files into one
                        fi

                    else
                        echo "There are multiple files beginning with year 2006. Make sure there are no duplicates of this data."
                    fi

            done
            
            ;;
        esac 

        cd $cwd # go back to base directory
        
    done

done

cd $cwd
