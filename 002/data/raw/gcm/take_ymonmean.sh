#!/usr/bin/env bash
set -euo pipefail

declare -a vars_gcm=("ta" "hur" "ps" "hurs" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss") # list of GCM variables that we want to process

cwd=$(pwd) # save current working directory
cd ../gcm_raw # switch to directory with raw data

# for dirs in */; do # loop through all the models
for dirs in MPI-ESM-LR/; do # for MPI-ESM-LR only (test run)
    mkdir -p $cwd/$dirs # make model directory in processed data folder if it doesn't exist yet
    cd $dirs # go in the model directory
    for vars in ${vars_gcm[@]}; do
        if ls $cwd/${dirs}${vars}_* 1> /dev/null 2>&1; then # check if data is already there
            echo "${vars} was already converted. Skipping..."
        else
            pattern="${vars}_*"
            if [ $(find . -maxdepth 1 -name "$pattern" -printf '.' | wc -m) -eq 0 ]; then # check if this file exists
                echo "File of type $pattern does not exist. Skipping..."
            elif [ $(find . -maxdepth 1 -name "$pattern" -printf '.' | wc -m) -eq 1 ]; then # check if there are multiple files of this variable
                file=$(ls $pattern)
                cdo ymonmean $file $cwd/$dirs$file.ymonmean.nc # take multi-year monthly climatology and save
            else
                files=( $pattern )
                first="${files[0]}" # first file of this variable
                last="${files[${#files[@]}-1]}" # last file of this variable
                common="${first%_*-*}_" # common part of file name that doesn't contain the years
                begin=${first#"${first%_*-*}_"} # extract just the year part of name containing first year
                end=${last#"${last%_*-*}_"} # extract just the year part of name containing last year
                yr_begin=${begin: 0:6} # extract the first year and month
                yr_end=${end: 7:6} # extract the last year and month
                merge="${common}${yr_begin}-${yr_end}.nc" # write file name with merged time
                cdo mergetime $(ls ${vars}_*) $cwd/$dirs/$merge # combine multiple files into one
                cdo ymonmean $cwd/$dirs/$merge $cwd/$dirs$merge.ymonmean.nc # take multi-year monthly climatology and save
                rm $cwd/$dirs/$merge # delete original merged file
            fi
        fi
    done
done

cd $cwd
