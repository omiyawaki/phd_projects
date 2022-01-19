#!/usr/bin/env bash
set -euo pipefail

# declare -a vars_gcm=("va" "zg" "wap" "ta" "hur" "ps" "hurs" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "rldscs" "rlutcs" "rsdscs" "rsuscs" "rsutcs" "hfls" "hfss" "pr" "prc" "evspsbl" "vas") # list of GCM variables that we want to process
# declare -a vars_gcm=("rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "prc" "evspsbl") # list of GCM variables that we want to process
# declare -a vars_gcm=("sfcWind") # list of GCM variables that we want to process
declare -a vars_gcm=("ps") # list of GCM variables that we want to process
# declare -a vars_gcm=("tas" "ts") # list of GCM variables that we want to process
# declare -a vars_gcm=("ps" "tas" "ta" "zg" "ts" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "evspsbl") # list of GCM variables that we want to process
# declare -a vars_gcm=("hurs") # list of GCM variables that we want to process
declare -a realm=("atmos")
#declare -a models=$(cd /project2/tas1/ockham/data9/tas/CMIP5_RAW && ls -d */) # list of GCM models to process
# declare -a models=$(cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85/ && ls -d */)
declare -a models=("CCSM4/") # list of GCM models to process
declare -a clim="rcp85" # climate name
#declare -a skip_models="MPI-ESM-LR/ MPI-ESM-MR/ MPI-ESM-P/ FGOALS-g2/ CanAM4/ CanCM4/ CESM1-CAM5-1-FV2/ CMCC-CESM/ CMCC-CM/ CMCC-CMS/ FIO-ESM/ HadGEM2-A/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ NICAM-09/"
declare -a skip_models="CESM1-WACCM/ CanAM4/ CanCM4/ CESM1-CAM5-1-FV2/ EC-EARTH/ FIO-ESM/ HadGEM2-A/ HadGEM2-AO/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ NICAM-09/ CNRM-CM5-2/ FGOALS-s2/ HadCM3/ MPI-ESM-P/"
# declare -a skip_models="ACCESS1-0/ ACCESS1-3/ bcc-csm1-1/ bcc-csm1-1-m/ BNU-ESM/ CanESM2/ CCSM4/ CESM1-BGC/ CESM1-CAM5/ CMCC-CESM/ CMCC-CM/ CMCC-CMS/ CNRM-CM5/ CSIRO-Mk3-6-0/ FGOALS-g2/ GFDL-CM3/ GFDL-ESM2G/ GFDL-ESM2M/ GISS-E2-H/ GISS-E2-H-CC/ GISS-E2-R/ GISS-E2-R-CC/ HadGEM2-CC/ HadGEM2-ES/ inmcm4/ IPSL-CM5A-LR/ IPSL-CM5A-MR/ IPSL-CM5B-LR/ MIROC5/ MIROC-ESM/ MIROC-ESM-CHEM/ MPI-ESM-LR/ MPI-ESM-MR/ MRI-CGCM3/ MRI-ESM1/ NorESM1-M/ NorESM1-ME/ CESM1-WACCM/ CanAM4/ CanCM4/ CESM1-CAM5-1-FV2/ EC-EARTH/ FIO-ESM/ HadGEM2-A/ HadGEM2-AO/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ NICAM-09/ CNRM-CM5-2/ FGOALS-s2/ HadCM3/ MPI-ESM-P/"
declare -a skip_files=("_eady.nc")

out_yr_begin="2070"
out_mn_begin="01"
out_dy_begin="01"
out_yr_end="2099"
out_mn_end="12"
out_dy_end="31"

cwd=$(pwd) # save current working directory
cd /project2/tas1/miyawaki/projects/003/data/raw/rcp85 # directory with data we will copy from
rwd=$(pwd) # save raw data directory

# for dirs in */; do # loop through all the models
for dirs in ${models[@]}; do # loop through models
    case $skip_models in *"$dirs"*)
        :
        ;; 
    *) 
        model=${dirs%/}
        echo $model
        mkdir -p $cwd/$dirs # make model directory in processed data folder if it doesn't exist yet
        # rm -f $cwd/$dirs*12.nc # remove nc files that are not ymonmean.nc
        # rm -f $cwd/$dirs*1231.nc
        # rm -f $cwd/$dirs*30.nc 
        cd $rwd/$dirs # go in the model directory
        echo $(pwd)
        for vars in ${vars_gcm[@]}; do
            echo $vars

            common=${vars}_Amon_${model}_${clim}_r1i1p1
            if [[ ("$model" == "HadGEM2-CC") && ( ("$vars" == "ta") || ("$vars" == "hur") || ("$vars" == "zg") ) ]]; then
                infile=${common}_200601-209912
            else
                infile=${common}_200601-210012
            fi
            outfile=${common}_${out_yr_begin}${out_mn_begin}-${out_yr_end}${out_mn_end}
            if [ ! -f "${infile}.nc" ]; then
                echo "${infile}.nc does not exist for ${model}. Skipping..."
            else
                cdo seldate,${out_yr_begin}-${out_mn_begin}-${out_dy_begin},${out_yr_end}-${out_mn_end}-${out_dy_end} ${infile}.nc ${cwd}/${model}/${outfile}.nc
                cdo ymonmean ${cwd}/${model}/${outfile}.nc ${cwd}/${model}/${outfile}.ymonmean.nc
                # rm ${cwd}/${model}/${outfile}.nc 
            fi

        done
        
        ;;
    esac 
    
done

cd $cwd
