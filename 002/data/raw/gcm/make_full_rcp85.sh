#!/usr/bin/env bash
set -euo pipefail

declare -a vars_gcm=("ps" "ta" "zg" "hus" "tas" "huss") # list of GCM variables that we want to process
# declare -a vars_gcm=("ps" "tas" "huss") # list of GCM variables that we want to process
# declare -a vars_gcm=("ta" "zg" "hus") # list of GCM variables that we want to process
# declare -a vars_gcm=("ps") # list of GCM variables that we want to process
#declare -a models=$(cd /project2/tas1/ockham/data9/tas/CMIP5_RAW && ls -d */) # list of GCM models to process
declare -a clim="rcp85" # climate name

###########################################################
# MODEL SUBSET 1
###########################################################
# declare -a models="bcc-csm1-1/"
# declare -a models="bcc-csm1-1-m/" # 2D VARS ONLY!!! (for 3D use subset 1.1)
# declare -a skip_files=("1-209912_eady.nc _eady.nc")

###########################################################
# MODEL SUBSET 1.1 
###########################################################
# declare -a models="bcc-csm1-1-m/" # 3D VARS ONLY!!! (for 2D use subset 1)
# declare -a skip_files=("200601-209912.nc 1-209912_eady.nc _eady.nc")

###########################################################
# MODEL SUBSET 2
###########################################################
declare -a models="MRI-ESM1/"
# declare -a models="ACCESS1-0/ ACCESS1-3/ BNU-ESM/ CanESM2/ CCSM4/ CESM1-BGC/ CESM1-CAM5/ CSIRO-Mk3-6-0/ GFDL-CM3/ FGOALS-g2/ inmcm4/ IPSL-CM5A-MR/ IPSL-CM5B-LR/ MIROC5/ MIROC-ESM/ MIROC-ESM-CHEM/ MPI-ESM-LR/ MPI-ESM-MR/ MRI-CGCM3/ MRI-ESM1/ NorESM1-M/ NorESM1-ME/"
declare -a skip_files=("200501-209912.nc 200501-210012.nc 200601-209912.nc _eady.nc")

###########################################################
# MODEL SUBSET 3
###########################################################
# declare -a models="CMCC-CESM/ CMCC-CM/ CMCC-CMS/ CNRM-CM5/ GFDL-ESM2G/ GFDL-ESM2M/ GISS-E2-H/ GISS-E2-H-CC/ GISS-E2-R/ GISS-E2-R-CC/ HadGEM2-CC/ HadGEM2-ES/ IPSL-CM5A-LR/"
# declare -a skip_files=("200601-209912.nc 200601-210012.nc _eady.nc")


###########################################################
# These models don't have RCP85 data
###########################################################
declare -a skip_models="CESM1-WACCM/ CanAM4/ CanCM4/ CESM1-CAM5-1-FV2/ CNRM-CM5-2/ EC-EARTH/ FGOALS-s2/ FIO-ESM/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ HadCM3/ HadGEM2-A/ HadGEM2-AO/ MPI-ESM-P/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ NICAM-09/"

cwd=$(pwd) # save current working directory
cd /project2/tas1/ockham/data9/tas/CMIP5_RAW/ # directory with data we will copy from
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

            # if ls $cwd/${dirs}${vars}_*${clim}*12.nc 1> /dev/null 2>&1; then # check if data is already there
            #     echo "${vars} was already merged. Skipping..."
            # else

                cd ./${clim}/atmos/mon/${vars}/r1i1p1/
                pattern="${vars}_*${dirs%/}*"
                files=( $pattern )
                echo ${files[@]}
                # remove files that are not original (list of non-original data is contained in variable $delete )
                for (( i=0; i<${#files[@]}; i++ )); do 
                    # case $skip_files in *"${files[i]#*_}"* | *"${files[i]:(-8)}"* | *"${files[i]:(-16)}"* )
                    case $skip_files in *"${files[i]#*_}"* | *"${files[i]:(-16)}"* )
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

            # fi

            cd $rwd/$dirs # go in the model directory
        done
        
        ;;
    esac 
    
done

cd $cwd
