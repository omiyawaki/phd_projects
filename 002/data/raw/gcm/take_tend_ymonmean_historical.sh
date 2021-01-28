#!/usr/bin/env bash
set -euo pipefail

declare -a vars_gcm=("vimse") # list of GCM variables that we want to process
declare -a models=$(cd /project2/tas1/ockham/data9/tas/CMIP5_RAW && ls -d */) # list of GCM models to process
declare -a clim="historical" # climate name
declare -a skip_models="MPI-ESM-P/ FGOALS-g2/ CanAM4/ CanCM4/ CESM1-CAM5-1-FV2/ CMCC-CESM/ CMCC-CM/ CMCC-CMS/ FIO-ESM/ HadGEM2-A/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ NICAM-09/"
declare -a skip_files=("Amon_bcc-csm1-1_historical_r1i1p1_185001-200512.nc Amon_bcc-csm1-1-m_historical_r1i1p1_185001-200512.nc Amon_GFDL-CM3_historical_r1i1p1_185001-200512.nc Amon_GISS-E2-H-CC_historical_r1i1p1_185001-200512.nc Amon_GISS-E2-R-CC_historical_r1i1p1_185001-200512.nc Amon_MIROC5_historical_r1i1p1_185001-200512.nc")

cwd=$(pwd) # save current working directory

# for dirs in */; do # loop through all the models
for dirs in ${models[@]}; do # loop through models
    case $skip_models in *"$dirs"*)
        :
        ;; 
    *) 
        echo $dirs
        mkdir -p $cwd/$dirs # make model directory in processed data folder if it doesn't exist yet
        for vars in ${vars_gcm[@]}; do
            echo $vars
            if ls $cwd/${dirs}${vars}_*historical*.nc 1> /dev/null 2>&1; then # check if data is already there
                echo "${vars} was already converted. Skipping..."
            else
                cd ./${dirs}tend
                pattern="${vars}_*${dirs%/}*"
                files=( $pattern )

                if [ ${#files[@]} -eq 0 ]; then # check if this file exists
                    echo "File of type $pattern does not exist. Please download the file and place it in the corresponding directory."
                    exit
                elif [ ${#files[@]} -eq 1 ]; then # check if there are multiple files of this variable
                    common="${files%_*-*}_" # common part of files name that doesn't contain the years
                    yr=${files#"${files%_*-*}_"} # extract just the year part of name containing last year
                    yr_end=${yr: 7:4} # extract the last year and month
                    sel="${common}197901-200512" # write files name with last 5 years
                    cdo seldate,1979-01,2006-01 $files $cwd/$dirs$sel.nc # select last 5 years
                    cdo deltat $cwd/$dirs$sel.nc $cwd/$dirs$sel.deltat.nc
                    cdo ymonmean $cwd/$dirs$sel.deltat.nc $cwd/$dirs$sel.deltat.ymonmean.nc # take multi-year monthly climatology and save
                    rm $cwd/$dirs/$sel.nc # delete temporary 5 year files
                    rm $cwd/$dirs/$sel.deltat.nc # delete temporary 5 year files
                else
                    first="${files[0]}" # first file of this variable
                    last="${files[${#files[@]}-1]}" # last file of this variable
                    common="${first%_*-*}_" # common part of file name that doesn't contain the years
                    begin=${first#"${first%_*-*}_"} # extract just the year part of name containing first year
                    end=${last#"${last%_*-*}_"} # extract just the year part of name containing last year
                    yr_begin=${begin: 0:4} # extract the first year and month
                    yr_end=${end: 7:4} # extract the last year and month
                    merge="${common}${yr_begin}01-${yr_end}12" # write file name with merged time
                    sel="${common}197901-200512" # write file name with selected years
                    # cdo mergetime $(ls ${vars}_*) $cwd/$dirs/$merge.nc # combine multiple files into one
                    cdo mergetime ${files[@]} $cwd/$dirs/$merge.nc # combine multiple files into one
                    cdo seldate,1979-01,2006-01 $cwd/$dirs/$merge.nc $cwd/$dirs/$sel.nc # select years 1979-2005
                    cdo deltat $cwd/$dirs$sel.nc $cwd/$dirs$sel.deltat.nc
                    cdo ymonmean $cwd/$dirs/$sel.deltat.nc $cwd/$dirs/$sel.deltat.ymonmean.nc # take multi-year monthly climatology and save
                    rm $cwd/$dirs/$merge.nc # delete temporary merged file
                    rm $cwd/$dirs/$sel.nc # delete temporary selected year file
                    rm $cwd/$dirs/$sel.deltat.nc # delete temporary 5 year files
                fi
            fi
            cd $cwd/$dirs # go in the model directory
        done
        
        ;;
    esac 
    
done

cd $cwd
