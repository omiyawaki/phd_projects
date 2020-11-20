#!/usr/bin/env bash
set -euo pipefail

# declare -a vars_gcm=("va" "zg" "wap" "ta" "hur" "ps" "hurs" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "rldscs" "rlutcs" "rsdscs" "rsuscs" "rsutcs" "hfls" "hfss" "pr" "prc" "evspsbl" "vas") # list of GCM variables that we want to process
declare -a vars_gcm=("ps" "tas" "ta" "zg" "ts" "hur") # list of GCM variables that we want to process
declare -a models=("MPI-ESM-LR") # list of GCMs to process

cwd=$(pwd) # save current working directory
cd ../historical_raw # switch to directory with raw data
rwd=$(pwd) # save raw data directory

# for dirs in */; do # loop through all the models
for dirs in ${models[@]}; do # loop through models
    mkdir -p $cwd/$dirs # make model directory in processed data folder if it doesn't exist yet
    cd $rwd/$dirs # go in the model directory
    for vars in ${vars_gcm[@]}; do
        if ls $cwd/${dirs}${vars}_*historical*.nc 1> /dev/null 2>&1; then # check if data is already there
            echo "${vars} was already converted. Skipping..."
        else
            pattern="${vars}_*"
            if [ $(find . -maxdepth 1 -name "$pattern" -printf '.' | wc -m) -eq 0 ]; then # check if this file exists
                echo "File of type $pattern does not exist. Skipping..."
            elif [ $(find . -maxdepth 1 -name "$pattern" -printf '.' | wc -m) -eq 1 ]; then # check if there are multiple files of this variable
                file=$(ls $pattern)
                common="${file%_*-*}_" # common part of file name that doesn't contain the years
                yr=${file#"${file%_*-*}_"} # extract just the year part of name containing last year
                yr_end=${yr: 7:4} # extract the last year and month
		yr_end_nlz=$(echo $yr_end | sed 's/^0*//') # remove leading zeros for arithmetic
		yr_endp1_nlz=$((yr_end_nlz+1)) # end year plus 1 (required because cdo seldate is non-inclusive in the end date)
                yr_m5_nlz=$((yr_end_nlz-4)) # 5 years before the end
		yr_m5=$(echo $yr_m5_nlz | sed -r ":r;s/\b[0-9]{1,$((4-1))}\b/0&/g;tr") # add leading zeros back in
                sel5="${common}${yr_m5}01-${yr_end}12" # write file name with last 5 years
                cdo seldate,${yr_m5_nlz}-01,$((yr_endp1_nlz))-01 $file $cwd/$dirs$sel5.nc # select last 5 years
                cdo ymonmean $cwd/$dirs$sel5.nc $cwd/$dirs$sel5.ymonmean.nc # take multi-year monthly climatology and save
                rm $cwd/$dirs/$sel5.nc # delete temporary 5 year file
            else
                files=( $pattern )
                first="${files[0]}" # first file of this variable
                last="${files[${#files[@]}-1]}" # last file of this variable
                common="${first%_*-*}_" # common part of file name that doesn't contain the years
                begin=${first#"${first%_*-*}_"} # extract just the year part of name containing first year
                end=${last#"${last%_*-*}_"} # extract just the year part of name containing last year
                yr_begin=${begin: 0:4} # extract the first year and month
                yr_end=${end: 7:4} # extract the last year and month
		yr_end_nlz=$(echo $yr_end | sed 's/^0*//') # remove leading zeros for arithmetic
		yr_endp1_nlz=$((yr_end_nlz+1)) # end year plus 1 (required because cdo seldate is non-inclusive in the end date)
                yr_m5_nlz=$((yr_end_nlz-4)) # 5 years before the end
		yr_m5=$(echo $yr_m5_nlz | sed -r ":r;s/\b[0-9]{1,$((4-1))}\b/0&/g;tr") # add leading zeros back in
                merge="${common}${yr_begin}01-${yr_end}12" # write file name with merged time
                sel5="${common}${yr_m5}01-${yr_end}12" # write file name with last 5 years
                cdo mergetime $(ls ${vars}_*) $cwd/$dirs/$merge.nc # combine multiple files into one
                cdo seldate,${yr_m5_nlz}-01,$((yr_endp1_nlz))-01 $cwd/$dirs/$merge.nc $cwd/$dirs/$sel5.nc # select last 5 years
                cdo ymonmean $cwd/$dirs/$sel5.nc $cwd/$dirs/$sel5.ymonmean.nc # take multi-year monthly climatology and save
                rm $cwd/$dirs/$merge.nc # delete temporary merged file
                rm $cwd/$dirs/$sel5.nc # delete temporary 5 year file
                if [ "$vars" = "wap" ]; then # we don't need full wap profile, just at 500 hPa
                    new5="${sel5//wap/w500}"
                    cdo sellevel,50000 $cwd/$dirs/$sel5.ymonmean.nc $cwd/$dirs/$new5.ymonmean.nc # extract 500 hPa data only
                    rm $cwd/$dirs/$sel5.ymonmean.nc # delete wap
                fi
            fi
        fi
    done
done

cd $cwd
