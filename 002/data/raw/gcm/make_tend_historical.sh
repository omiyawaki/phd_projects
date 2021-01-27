#!/usr/bin/env bash
set -euo pipefail

#declare -a models=$(cd /project2/tas1/ockham/data9/tas/CMIP5_RAW && ls -d */) # list of GCM models to process
declare -a models=("MPI-ESM-LR/") # list of GCM models to process
declare -a clim="historical" # climate name
# declare -a skip_models="MPI-ESM-LR/ MPI-ESM-MR/ MPI-ESM-P/ FGOALS-g2/ CanAM4/ CanCM4/ CESM1-CAM5-1-FV2/ CMCC-CESM/ CMCC-CM/ CMCC-CMS/ FIO-ESM/ HadGEM2-A/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ NICAM-09/"
declare -a skip_models="MPI-ESM-P/ FGOALS-g2/ CanAM4/ CanCM4/ CESM1-CAM5-1-FV2/ CMCC-CESM/ CMCC-CM/ CMCC-CMS/ FIO-ESM/ HadGEM2-A/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ NICAM-09/"
declare -a skip_files=("Amon_bcc-csm1-1_historical_r1i1p1_185001-200512.nc Amon_bcc-csm1-1-m_historical_r1i1p1_185001-200512.nc Amon_GFDL-CM3_historical_r1i1p1_185001-200512.nc Amon_GISS-E2-H-CC_historical_r1i1p1_185001-200512.nc Amon_GISS-E2-R-CC_historical_r1i1p1_185001-200512.nc Amon_MIROC5_historical_r1i1p1_185001-200512.nc")

cwd=$(pwd) # save current working directory
cd ../${clim}_raw # switch to directory with raw data
rwd=$(pwd) # save raw data directory

# for dirs in */; do # loop through all the models
for dirs in ${models[@]}; do # loop through models
    case $skip_models in *"$dirs"*)
        :
        ;; 
    *) 
        echo $dirs
        mkdir -p $cwd/$dirs # make model directory in processed data folder if it doesn't exist yet
        cd $rwd/$dirs # go in the model directory

        if ls $cwd/${dirs}tend_*$clim*.nc 1> /dev/null 2>&1; then # check if data is already there
            echo "tend was already created. Skipping..."
        else
            cd ./${clim}/atmos/day/ta/r1i1p1/
            prepath_ta=$(pwd)
            pattern="ta_*${dirs%/}*"
            files=( $pattern )
            # remove files that are not original (list of non-original data is contained in variable $delete )
            for (( i=0; i<${#files[@]}; i++ )); do 
                case $skip_files in *"${files[i]#*_}"*)
                    files=( "${files[@]:0:$i}" "${files[@]:$((i + 1))}" )
                    i=$((i - 1))
                    ;;
                *)
                    :
                    ;;
                esac
            done

            if [ ${#files[@]} -eq 0 ]; then # check if this file exists
                echo "File of type $pattern does not exist. Please download the file and place it in the corresponding directory."
                exit
            fi
            
            for fn_ta in ${files[@]}; do
                
                prepath_hus="${prepath_ta/\/ta\///hus/}"
                prepath_zg="${prepath_ta/\/ta\///zg/}"
                
                fn_hus="${fn_ta/ta_/hus_}"
                fn_zg="${fn_ta/ta_/zg_}"
                fn_zg="${fn_ta/ta_/zg_}"
                
                full_ta=$prepath_ta/$fn_ta
                full_hus=$prepath_hus/$fn_hus
                full_zg=$prepath_zg/$fn_zg
                
                # create pressure data 
                fn_ps="${fn_ta/ta_/ps_}"
                common="${fn_ps%_*-*}_" # common part of file name that doesn't contain the years
                begin=${fn_ps#"${fn_ps%_*-*}_"} # extract just the year part of name containing first year
                yr_begin=${begin: 0:4} # extract the first year and month
                mn_begin=${begin: 4:2} # extract the first year and month
                yr_end=${begin: 9:4} # extract the last year and month
                mn_end=${begin: 13:2} # extract the last year and month
                
                ps_in=$(cd ${cwd}/${dirs} && ls ps_*${clim}*12.nc)
                begin_in=${ps_in#"${ps_in%_*-*}_"} # extract just the year part of name containing first year
                ps_out=${ps_in/${begin_in}/${begin}}
                full_ps=${cwd}/${dirs}$ps_out
                
                cdo seldate,${yr_begin}-${mn_begin},${yr_end}-${mn_end} ${cwd}/${dirs}${ps_in} ${full_ps} # extract out same year of ps data
                
                python ${cwd}/make_tend.py ${full_ps} ${full_ta} ${full_zg} ${full_hus}
                
                rm ${cwd}/${dirs}$ps_out
            done
        fi
        cd $rwd/$dirs # go in the model directory
        
        ;;
    esac 
    
done

cd $cwd
