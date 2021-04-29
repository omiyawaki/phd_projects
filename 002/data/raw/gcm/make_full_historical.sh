#!/usr/bin/env bash
set -euo pipefail

declare -a vars_gcm=("ps" "ta" "zg" "hur") # list of GCM variables that we want to process
# declare -a vars_gcm=("zg") # list of GCM variables that we want to process
#declare -a models=$(cd /project2/tas1/ockham/data9/tas/CMIP5_RAW && ls -d */) # list of GCM models to process
declare -a models=("HadCM3/") # list of GCM models to process
declare -a clim="historical" # climate name
declare -a skip_models="EC-EARTH/ CanAM4/ CanCM4/ CESM1-CAM5-1-FV2/ CMCC-CMS/ FIO-ESM/ HadGEM2-A/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ NICAM-09/"
declare -a skip_files=("_eady.nc")

cwd=$(pwd) # save current working directory
cd ../historical_raw # switch to directory with raw data
rwd=$(pwd) # save raw data directory

#for dirs in */; do # loop through all the models
for dirs in ${models[@]}; do # loop through models
    case $skip_models in *"$dirs"*)
        :
        ;; 
    *) 
        echo $dirs
        mkdir -p $cwd/$dirs # make model directory in processed data folder if it doesn't exist yet
        cd $rwd/$dirs # go in the model directory
        for vars in ${vars_gcm[@]}; do
            echo $vars

            if ls $cwd/${dirs}${vars}_*historical*12.nc 1> /dev/null 2>&1; then # check if data is already there
                echo "${vars} was already merged. Skipping..."
            else

                cd ./${clim}/atmos/mon/${vars}/r1i1p1/
                pattern="${vars}_*${dirs%/}*"
                files=( $pattern )
                echo ${files[@]}
                # remove files that are not original (list of non-original data is contained in variable $delete )
                for (( i=0; i<${#files[@]}; i++ )); do 
                    case $skip_files in *"${files[i]#*_}"* | *"${files[i]:(-8)}")
                        files=( "${files[@]:0:$i}" "${files[@]:$((i + 1))}" )
                        i=$((i - 1))
                        ;;
                    *)
                        :
                        ;;
                    esac
                done
                echo ${files[@]}

                if [[ ${files[@]} != *.nc ]]; then # check if this file exists
                    echo "File of type $pattern does not exist. Please download the file and place it in the corresponding directory."
                elif [ ${#files[@]} -eq 1 ]; then # check if there are multiple files of this variable
                    cp $files $cwd/$dirs # no need to merge, just copy over
                else
                    first="${files[0]}" # first file of this variable
                    last="${files[${#files[@]}-1]}" # last file of this variable
                    common="${first%_*-*}_" # common part of file name that doesn't contain the years
                    begin=${first#"${first%_*-*}_"} # extract just the year part of name containing first year
                    end=${last#"${last%_*-*}_"} # extract just the year part of name containing last year
                    yr_begin=${begin: 0:4} # extract the first year and month
                    yr_end=${end: 7:4} # extract the last year and month
                    merge="${common}${yr_begin}01-${yr_end}12" # write file name with merged time
                    cdo -O mergetime ${files[@]} $cwd/$dirs$merge.nc # combine multiple files into one
                fi

            fi

            cd $rwd/$dirs # go in the model directory
        done
        
        ;;
    esac 
    
done

cd $cwd
