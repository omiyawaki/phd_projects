#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python
source activate netcdf

#declare -a models=$(cd /project2/tas1/ockham/data9/tas/CMIP5_RAW && ls -d */) # list of GCM models to process
# declare -a models=("ACCESS1-0/") # list of GCM models to process
declare -a models=("MPI-ESM-LR/") # list of GCM models to process
declare -a clim="historical" # climate name
declare -a freq="day" # compute mse tendency using daily or monthly data?
# declare -a skip_models="FGOALS-s2/ FGOALS-g2/  CCSM4/ MPI-ESM-P/ FGOALS-g2/ CanAM4/ CanCM4/ CESM1-CAM5-1-FV2/ CMCC-CESM/ CMCC-CM/ CMCC-CMS/ FIO-ESM/ HadGEM2-A/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ NICAM-09/"
declare -a skip_models="ACCESS1-0/ ACCESS1-3/ bcc-csm1-1/ bcc-csm1-1-m/ BNU-ESM/ CanESM2/ CCSM4/ CESM1-BGC/ CESM1-CAM5/ CESM1-WACCM/ CNRM-CM5/ CNRM-CM5-2/ CSIRO-Mk3-6-0/ EC-EARTH/ MPI-ESM-P/ FGOALS-g2/ FGOALS-s2/ GFDL-CM3/ GFDL-ESM2G/ GFDL-ESM2M/ GISS-E2-H/ GISS-E2-H-CC/ GISS-E2-R/ GISS-E2-R-CC/ HadCM3/ HadGEM2-AO/ HadGEM2-CC/ HadGEM2-ES/ CanAM4/ CanCM4/ CESM1-CAM5-1-FV2/ CMCC-CESM/ CMCC-CM/ CMCC-CMS/ FIO-ESM/ HadGEM2-A/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ NICAM-09/"
declare -a skip_files=("Amon_bcc-csm1-1_historical_r1i1p1_185001-200512.nc Amon_bcc-csm1-1-m_historical_r1i1p1_185001-200512.nc Amon_GFDL-CM3_historical_r1i1p1_185001-200512.nc Amon_GISS-E2-H-CC_historical_r1i1p1_185001-200512.nc Amon_GISS-E2-R-CC_historical_r1i1p1_185001-200512.nc Amon_MIROC5_historical_r1i1p1_185001-200512.nc 197901-200512.nc")
declare -a no_huss=("MPI-ESM-LR/ MPI-ESM-MR/ MPI-ESM-P/") # models that don't have huss data

cwd=$(pwd) # save current working directory
cd ../${clim}_raw # switch to directory with raw data
rwd=$(pwd) # save raw data directory

for dirs in */; do # loop through all the models
# for dirs in ${models[@]}; do # loop through models
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
            cd ./${clim}/atmos/${freq}/ta/r1i1p1/
            prepath_ta=$(pwd)
            pattern="ta_*${dirs%/}*"
            files=( $pattern )
            # remove files that are not original (list of non-original data is contained in variable $delete )
            for (( i=0; i<${#files[@]}; i++ )); do 
                case $skip_files in *"${files[i]#*_}"* | *"${files[i]:(-16)}")
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
                echo "Creating tendency for dates ${fn_ta:(-20):(-3)}..."
                
                prepath_hus="${prepath_ta/\/ta\///hus/}"
                prepath_zg="${prepath_ta/\/ta\///zg/}"
                
                fn_hus="${fn_ta/ta_/hus_}"
                fn_zg="${fn_ta/ta_/zg_}"
                
                full_ta=$prepath_ta/$fn_ta
                full_hus=$prepath_hus/$fn_hus
                full_zg=$prepath_zg/$fn_zg
                full_orog=${rwd}/${dirs}${clim}/atmos/fx/orog/r0i0p0/orog_*.nc
                
                # import monthly pressure data 
                fn_ps="${fn_ta/ta_/ps_}"
                common="${fn_ps%_*-*}_" # common part of file name that doesn't contain the years
                begin=${fn_ps#"${fn_ps%_*-*}_"} # extract just the year part of name containing first year
                yr_begin=${begin: 0:4} # extract the first year and month
                mn_begin=${begin: 4:2} # extract the first year and month
                yr_end=${begin: 9:4} # extract the last year and month
                mn_end=${begin: 13:2} # extract the last year and month
                
                ps_in=$(cd ${cwd}/${dirs} && ls ps_*${clim}*12.nc | grep -v 197901-200512.nc)
                begin_in=${ps_in#"${ps_in%_*-*}_"} # extract just the year part of name containing first year
                ps_out=${ps_in/${begin_in}/${begin}}
                full_ps=${cwd}/${dirs}$ps_out
                cdo seldate,${yr_begin}-${mn_begin}-01,${yr_end}-${mn_end}-31 ${cwd}/${dirs}${ps_in} ${full_ps} # extract out same year of ps data

                # import monthly tas data
                fn_tas="${fn_ta/ta_/tas_}"
                common="${fn_tas%_*-*}_" # common part of file name that doesn't contain the years
                begin=${fn_tas#"${fn_tas%_*-*}_"} # extract just the year part of name containing first year
                yr_begin=${begin: 0:4} # extract the first year and month
                mn_begin=${begin: 4:2} # extract the first year and month
                yr_end=${begin: 9:4} # extract the last year and month
                mn_end=${begin: 13:2} # extract the last year and month
                
                tas_in=$(cd ${cwd}/${dirs} && ls tas_*${clim}*12.nc | grep -v 197901-200512.nc)
                begin_in=${tas_in#"${tas_in%_*-*}_"} # extract just the year part of name containing first year
                tas_out=${tas_in/${begin_in}/${begin}}
                full_tas=${cwd}/${dirs}$tas_out
                
                cdo seldate,${yr_begin}-${mn_begin}-01,${yr_end}-${mn_end}-31 ${cwd}/${dirs}${tas_in} ${full_tas} # extract out same year of tas data
                
                # import monthly huss data if available
                case ${no_huss} in *"${dirs}"*) # if not available then interpolate huss from 3d hus
                    full_huss=${full_hus}
                    :
                    ;;
                *)
                    fn_huss="${fn_ta/ta_/huss_}"
                    common="${fn_huss%_*-*}_" # common part of file name that doesn't contain the years
                    begin=${fn_huss#"${fn_huss%_*-*}_"} # extract just the year part of name containing first year
                    yr_begin=${begin: 0:4} # extract the first year and month
                    mn_begin=${begin: 4:2} # extract the first year and month
                    yr_end=${begin: 9:4} # extract the last year and month
                    mn_end=${begin: 13:2} # extract the last year and month
                    
                    huss_in=$(cd ${cwd}/${dirs} && ls huss_*${clim}*12.nc | grep -v 197901-200512.nc)
                    begin_in=${huss_in#"${huss_in%_*-*}_"} # extract just the year part of name containing first year
                    huss_out=${huss_in/${begin_in}/${begin}}
                    full_huss=${cwd}/${dirs}$huss_out

                    cdo seldate,${yr_begin}-${mn_begin}-01,${yr_end}-${mn_end}-31 ${cwd}/${dirs}${huss_in} ${full_huss} # extract out same year of huss data
                    ;;
                esac
                
                newdir=${cwd}/${dirs}tend
                if [ ! -d "$newdir" ]; then # create directory if it doesn't exist
                    mkdir -p $newdir
                fi
                fn_tend="${fn_ta/ta_/tend_}"
                full_tend=${newdir}/${fn_tend}
                
                # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_tend_mon.py ${full_ps} ${full_ta} ${full_zg} ${full_hus} ${full_tend} ${full_tas} ${full_orog} ${full_huss}
                python ${cwd}/make_tend_alt.py ${full_ps} ${full_ta} ${full_zg} ${full_hus} ${full_tend} ${full_tas} ${full_orog} ${full_huss}
                
            done
        fi
        cd $rwd/$dirs # go in the model directory
        
        ;;
    esac 
    
done

cd $cwd
