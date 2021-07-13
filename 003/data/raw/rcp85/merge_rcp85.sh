#!/usr/bin/env bash
set -euo pipefail

# declare -a vars_gcm=("zg" "ta" "hur" "ps" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "prc" "evspsbl") # list of GCM variables that we want to process
declare -a vars_gcm=("rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss") # list of GCM variables that we want to process
# declare -a vars_gcm=("ps" "tas" "ta" "zg" "ts" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "evspsbl") # list of GCM variables that we want to process
# declare -a vars_gcm=("prc") # list of GCM variables that we want to process
declare -a realm=("atmos")
declare -a clim="rcp85" # climate name
declare -a freq="mon" # data output frequency (e.g. fx for fixed, mon for monthly, day for daily)
declare -a ens="r1i1p1" # ensemble specification 
#declare -a skip_models="MPI-ESM-LR/ MPI-ESM-MR/ MPI-ESM-P/ FGOALS-g2/ CanAM4/ CanCM4/ CESM1-CAM5-1-FV2/ CMCC-CESM/ CMCC-CM/ CMCC-CMS/ FIO-ESM/ HadGEM2-A/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ NICAM-09/"
# declare -a models=$(cd /project2/tas1/ockham/data9/tas/CMIP5_RAW && ls -d */) # list of GCM models to process
declare -a models=("CCSM4/") # list of GCM models to process
# declare -a models_noskip="bcc-csm1-1/ bcc-csm1-1-m/"
# declare -a skip_files=("_eady.nc")
# declare -a models_skip2099="ACCESS1-0/ ACCESS1-3/ BNU-ESM/ CanESM2/ CCSM4/ CESM1-BGC/ CESM1-CAM5/ CSIRO-Mk3-6-0/ GFDL-CM3/ FGOALS-g2/ inmcm4/ IPSL-CM5A-MR/ IPSL-CM5B-LR/ MIROC5/ MIROC-ESM/ MIROC-ESM-CHEM/ MPI-ESM-LR/ MPI-ESM-MR/ MRI-CGCM3/ MRI-ESM1/ NorESM1-M/ NorESM1-ME/"
declare -a skip_files=("200501-209912.nc 200601-209912.nc _eady.nc")
# declare -a models_skip2100="CMCC-CESM/ CMCC-CM/ CMCC-CMS/ CNRM-CM5/ GFDL-ESM2G/ GFDL-ESM2M/ GISS-E2-H/ GISS-E2-H-CC/ GISS-E2-R/ GISS-E2-R-CC/ HadGEM2-CC/ HadGEM2-ES/ IPSL-CM5A-LR/"
# declare -a skip_files=("200601-209912.nc 200601-210012.nc _eady.nc")
declare -a skip_models="CESM1-WACCM/ CanAM4/ CanCM4/ CESM1-CAM5-1-FV2/ CNRM-CM5-2/ EC-EARTH/ FGOALS-s2/ FIO-ESM/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ HadCM3/ HadGEM2-A/ HadGEM2-AO/ MPI-ESM-P/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ NICAM-09/"
# declare -a skip_models="ACCESS1-0/ ACCESS1-3/ bcc-csm1-1/ bcc-csm1-1-m/ BNU-ESM/ CanESM2/ CCSM4/ CESM1-BGC/ CESM1-CAM5/ CMCC-CESM/ CMCC-CM/ CMCC-CMS/ CNRM-CM5/ CSIRO-Mk3-6-0/ GFDL-CM3/ GFDL-ESM2G/ GFDL-ESM2M/ FGOALS-g2/ GISS-E2-H/ GISS-E2-H-CC/ GISS-E2-R/ GISS-E2-R-CC/ HadGEM2-CC/ HadGEM2-ES/ inmcm4/ IPSL-CM5A-LR/ IPSL-CM5A-MR/ IPSL-CM5B-LR/ MIROC5/ MIROC-ESM/ MIROC-ESM-CHEM/ MPI-ESM-LR/ MPI-ESM-MR/ MRI-CGCM3/ MRI-ESM1/ NorESM1-M/ CESM1-WACCM/ CanAM4/ CanCM4/ CESM1-CAM5-1-FV2/ CNRM-CM5-2/ EC-EARTH/ FGOALS-s2/ FIO-ESM/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ HadCM3/ HadGEM2-A/ HadGEM2-AO/ MPI-ESM-P/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ NICAM-09/"
# declare -a skip_models="ACCESS1-0/ ACCESS1-3/ bcc-csm1-1/ bcc-csm1-1-m/ BNU-ESM/ CanAM4/ CanCM4/ CanESM2/ CCSM4/ CESM1-BGC/ CESM1-CAM5/ CESM1-CAM5-1-FV2/ CESM1-WACCM/ CMCC-CESM/ CMCC-CM/ CMCC-CMS/ CNRM-CM5/ CNRM-CM5-2/ CSIRO-Mk3-6-0/ EC-EARTH/ FGOALS-g2/ FGOALS-s2/ FIO-ESM/ GFDL-CM3/ GFDL-ESM2G/ GFDL-ESM2M/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ GISS-E2-H/ GISS-E2-H-CC/ GISS-E2-R/ GISS-E2-R-CC/ HadCM3/ HadGEM2-A/ HadGEM2-AO/ HadGEM2-CC/ HadGEM2-ES/ inmcm4/ IPSL-CM5A-LR/ IPSL-CM5A-MR/ IPSL-CM5B-LR/ MIROC5/ MIROC-ESM/ MIROC-ESM-CHEM/ MPI-ESM-LR/ MPI-ESM-MR/ MPI-ESM-P/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ MRI-CGCM3/ MRI-ESM1/ NorESM1-M/ NICAM-09/"
# declare -a skip_models="bcc-csm1-1/ bcc-csm1-1-m/ CCSM4/ CESM1-BGC/ CESM1-WACCM/ FGOALS-g2/ NorESM1-M/ NorESM1-ME/ CanAM4/ CanCM4/ CMCC-CMS/ CESM1-CAM5-1-FV2/ EC-EARTH/ FIO-ESM/ HadGEM2-A/ HadGEM2-AO/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ NICAM-09/"

cwd=$(pwd) # save current working directory
cd /project2/tas1/ockham/data9/tas/CMIP5_RAW # switch to directory with raw data
rwd=$(pwd) # save raw data directory

# for dirs in */; do # loop through all the models
for dirs in ${models[@]}; do # loop through models
    case $skip_models in *"$dirs"*)
        :
        ;; 
    *) 
        echo $dirs
        mkdir -p $cwd/$dirs # make model directory in processed data folder if it doesn't exist yet
        # rm -f $cwd/$dirs*12.nc # remove nc files that are not ymonmean.nc
        # rm -f $cwd/$dirs*1231.nc
        # rm -f $cwd/$dirs*30.nc 
        cd $rwd/$dirs # go in the model directory
        echo $(pwd)
        for vars in ${vars_gcm[@]}; do
            echo $vars
            # if ls $cwd/${dirs}${vars}_*historical*${out_yr_begin}${out_mn_begin}-${out_yr_end}${out_mn_end}*.nc 1> /dev/null 2>&1; then # check if data is already there
            #     echo "${vars} was already converted. Skipping..."
            # else
                cd ./${clim}/${realm}/${freq}/${vars}/${ens}/
                pattern="${vars}_*${dirs%/}*$ens*"
                files=( $pattern )
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
                    # common="${files%_*-*}_" # common part of files name that doesn't contain the years
                    # yr=${files#"${files%_*-*}_"} # extract just the year part of name containing last year
                    # yr_end=${yr: 7:4} # extract the last year and month
                    # sel="${common}${out_yr_begin}${out_mn_begin}-${out_yr_end}${out_mn_end}" # write files name with last 5 years
                    cp ${files[@]} $cwd/$dirs
                    cdo -O yearmean $cwd/$dirs${files[@]} $cwd/$dirs${files[@]%.nc}.yearmean.nc # take annual mean
                else
                    first="${files[0]}" # first file of this variable
                    last="${files[${#files[@]}-1]}" # last file of this variable
                    common="${first%_*-*}_" # common part of file name that doesn't contain the years
                    begin=${first#"${first%_*-*}_"} # extract just the year part of name containing first year
                    end=${last#"${last%_*-*}_"} # extract just the year part of name containing last year
                    yr_begin=${begin: 0:4} # extract the first year and month
                    yr_end=${end: 7:4} # extract the last year and month
                    merge="${common}${yr_begin}01-${yr_end}12" # write file name with merged time
                    # cdo mergetime $(ls ${vars}_*) $cwd/$dirs/$merge.nc # combine multiple files into one
                    cdo -O mergetime ${files[@]} $cwd/$dirs$merge.nc # combine multiple files into one
                    cdo -O yearmean $cwd/$dirs$merge.nc $cwd/$dirs$merge.yearmean.nc # combine multiple files into one
                # fi
            fi
            cd $rwd/$dirs # go in the model directory
        done
        
        ;;
    esac 
    
done

cd $cwd
