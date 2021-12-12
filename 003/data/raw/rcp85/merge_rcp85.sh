#!/usr/bin/env bash
set -euo pipefail

# declare -a vars_gcm=("zg" "ta" "hur" "ps" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "prc" "evspsbl") # list of GCM variables that we want to process
# declare -a vars_gcm=("rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss") # list of GCM variables that we want to process
# declare -a vars_gcm=("rlutcs" "rsutcs" "rldscs" "rsdscs" "rsuscs") # list of GCM variables that we want to process
declare -a vars_gcm=("ta" "hus" "zg") # list of GCM variables that we want to process
# declare -a vars_gcm=("ps" "tas" "ts" "hurs") # list of GCM variables that we want to process
# declare -a vars_gcm=("clt" "clwvi") # list of GCM variables that we want to process
# declare -a vars_gcm=("hus") # list of GCM variables that we want to process
# declare -a vars_gcm=("tas" "ts") # list of GCM variables that we want to process
declare -a realm=("atmos")
declare -a clim="rcp85" # climate name
declare -a freq="mon" # data output frequency (e.g. fx for fixed, mon for monthly, day for daily)
declare -a ens="r1i1p1" # ensemble specification 
# declare -a ens="r0i0p0" # ensemble specification 
# declare -a models=("CSIRO-Mk3-6-0/") # extended RCP runs
# declare -a models=("HadGEM2-ES/" "bcc-csm1-1/" "CCSM4/" "CNRM-CM5/" "CSIRO-Mk3-6-0/" "IPSL-CM5A-LR/" "MPI-ESM-LR/") # extended RCP runs
# declare -a skip_files=("_eady.nc")

# declare -a models=$(cd /project2/tas1/ockham/data9/tas/CMIP5_RAW && ls -d */) # list of GCM models to process
# declare -a models=("MPI-ESM-LR/") # list of GCM models to process

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
# declare -a models="ACCESS1-0/ ACCESS1-3/ BNU-ESM/ CanESM2/ CCSM4/ CNRM-CM5/ CESM1-BGC/ CESM1-CAM5/ CSIRO-Mk3-6-0/ GFDL-CM3/ FGOALS-g2/ inmcm4/ IPSL-CM5A-MR/ IPSL-CM5B-LR/ MIROC5/ MIROC-ESM/ MIROC-ESM-CHEM/ MPI-ESM-LR/ MPI-ESM-MR/ MRI-CGCM3/ MRI-ESM1/ NorESM1-M/ NorESM1-ME/"
# declare -a models="CCSM4/ CNRM-CM5/ CESM1-BGC/ CESM1-CAM5/ CSIRO-Mk3-6-0/ GFDL-CM3/ FGOALS-g2/ inmcm4/ IPSL-CM5A-MR/ IPSL-CM5B-LR/ MIROC5/ MIROC-ESM/ MIROC-ESM-CHEM/ MPI-ESM-LR/ MPI-ESM-MR/ MRI-CGCM3/ MRI-ESM1/ NorESM1-M/ NorESM1-ME/"
# declare -a skip_files=("200501-209912.nc 200501-210012.nc 200601-209912.nc _eady.nc")

###########################################################
# MODEL SUBSET 3
###########################################################
# declare -a models="CMCC-CESM/ CMCC-CM/ CMCC-CMS/ GFDL-ESM2G/ GFDL-ESM2M/ GISS-E2-H/ GISS-E2-H-CC/ GISS-E2-R/ GISS-E2-R-CC/ HadGEM2-CC/ HadGEM2-ES/ IPSL-CM5A-LR/"
declare -a models="IPSL-CM5A-LR/"
declare -a skip_files=("200601-209912.nc 200601-210012.nc _eady.nc")


###########################################################
# These models don't have RCP85 data
###########################################################
declare -a skip_models="CESM1-WACCM/ CanAM4/ CanCM4/ CESM1-CAM5-1-FV2/ CNRM-CM5-2/ EC-EARTH/ FGOALS-s2/ FIO-ESM/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ HadCM3/ HadGEM2-A/ HadGEM2-AO/ MPI-ESM-P/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ NICAM-09/"

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
                    common="${files%_*-*}_" # common part of files name that doesn't contain the years
                    yr=${files#"${files%_*-*}_"} # extract just the year part of name containing last year
                    yr_begin=${yr: 0:4} # extract the first year and month
                    echo $yr_begin
                    yr_end=${yr: 7:4} # extract the last year and month
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

            #######################################################################
            # convert curvilinear to standard lat lon grid for sea ice data
            #######################################################################
            if [ ${vars} == "sic" ]; then
                ref_file=$(ls ${cwd}/${dirs}tas_*${yr_begin}01-${yr_end}12.nc)
                sic_file=$(ls ${cwd}/${dirs}sic_*${yr_begin}01-${yr_end}12.nc)

                # rename sic file in original grid
                mv ${sic_file} ${sic_file%.nc}.origgrid.nc

                # first create file containing standard lat-lon grid data (e.g., using tas file)
                cdo griddes ${ref_file} > ${cwd}/${dirs}grid_latlon
                sed -i "s/generic/lonlat/g" ${cwd}/${dirs}grid_latlon

                # convert to lat lon
                cdo -remapbil,${cwd}/${dirs}grid_latlon ${sic_file%.nc}.origgrid.nc ${sic_file}
                
            fi

            cd $rwd/$dirs # go in the model directory
        done
        
        ;;
    esac 
    
done

cd $cwd
