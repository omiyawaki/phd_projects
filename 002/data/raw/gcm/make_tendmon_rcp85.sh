#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python
source activate netcdf

#declare -a models=$(cd /project2/tas1/ockham/data9/tas/CMIP5_RAW && ls -d */) # list of GCM models to process
declare -a models=("ACCESS1-0/") # list of GCM models to process
# declare -a models=("bcc-csm1-1/ bcc-csm1-1-m/ BNU-ESM/ CanESM2/ CCSM4/ CESM1-BGC/ CESM1-CAM5/ CESM1-WACCM/ CMCC-CESM/ CMCC-CM/ CNRM-CM5/ CNRM-CM5-2/ CSIRO-Mk3-6-0/ FGOALS-g2/ FGOALS-s2/ GFDL-CM3/ GFDL-ESM2M/ GFDL-ESM2G/ GISS-E2-H/ GISS-E2-H-CC/ GISS-E2-R/ GISS-E2-R-CC/ HadGEM2-CC/ HadGEME2-ES/ inmcm4/ IPSL-CM5A-LR/ IPSL-CM5A-MR/ IPSL-CM5B-LR/ MIROC5/ MIROC-ESM/ MIROC-ESM-CHEM/ MPI-ESM-LR/ MPI-ESM-MR/ MPI-ESM-P/ MRI-CGCM3/ MRI-ESM1/ NorESM1-M/ NorESM1-ME/") # list of GCM models to process
declare -a clim="rcp85" # climate name
declare -a freq="mon" # compute mse tendency using daily or monthly data?
declare -a skip_models="CESM1-WACCM/ CanAM4/ CanCM4/ CESM1-CAM5-1-FV2/ EC-EARTH/ FIO-ESM/ HadGEM2-A/ HadGEM2-AO/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ NICAM-09/ CNRM-CM5-2/ FGOALS-s2/ HadCM3/ MPI-ESM-P/"
declare -a skip_files=("_eady.nc")
declare -a no_huss=("MPI-ESM-LR/ MPI-ESM-MR/ MPI-ESM-P/") # models that don't have huss data

yr_begin="2070" # extract the first year and month
mn_begin="01" # extract the first year and month
yr_end="2099" # extract the last year and month
mn_end="12" # extract the last year and month

cwd=$(pwd) # save current working directory
cd /project2/tas1/ockham/data9/tas/CMIP5_RAW/ # directory with data we will copy from
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

        if [[ $freq == "day" ]]; then
            checkstring=$cwd/${dirs}tend_*$clim*.nc
        elif [[ $freq == "mon" ]]; then
            checkstring=$cwd/${dirs}tendmon_*$clim*.nc
        fi

        if ls $checkstring 1> /dev/null 2>&1; then # check if data is already there
            echo "tend was already created. Skipping..."
        else
            cd ./${clim}/atmos/${freq}/ta/r1i1p1/
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
            
            fn_ta=ta_Amon_${dirs%/}_${clim}_r1i1p1_${yr_begin}${mn_begin}-${yr_end}${mn_end}.nc
            
            prepath_hus="${prepath_ta/\/ta\///hus/}"
            prepath_zg="${prepath_ta/\/ta\///zg/}"
            
            fn_hus="${fn_ta/ta_/hus_}"
            fn_zg="${fn_ta/ta_/zg_}"
            
            full_orog="null" #${rwd}/${dirs}${clim}/atmos/fx/orog/r0i0p0/orog_*.nc
            
            # import monthly pressure data 
            ps_in=$(cd ${cwd}/${dirs} && ls --ignore=$ignorestr ps_*${clim}*12.nc)
            begin_in=${ps_in#"${ps_in%_*-*}_"} # extract just the year part of name containing first year
            ps_out=ps_Amon_${dirs%/}_${clim}_r1i1p1_${yr_begin}${mn_begin}-${yr_end}${mn_end}.nc
            full_ps=${cwd}/${dirs}$ps_out

            if [[ ! -f $full_ps ]]; then
                cdo seldate,${yr_begin}-${mn_begin}-01,${yr_end}-${mn_end}-31 ${cwd}/${dirs}${ps_in} ${full_ps} # extract out same year of ps data
            fi 

            # import monthly tas data
            tas_in=$(cd ${cwd}/${dirs} && ls tas_*${clim}*12.nc)
            begin_in=${tas_in#"${tas_in%_*-*}_"} # extract just the year part of name containing first year
            tas_out=tas_Amon_${dirs%/}_${clim}_r1i1p1_${yr_begin}${mn_begin}-${yr_end}${mn_end}.nc
            full_tas=${cwd}/${dirs}$tas_out

            if [[ ! -f $full_tas ]]; then
                cdo seldate,${yr_begin}-${mn_begin}-01,${yr_end}-${mn_end}-31 ${cwd}/${dirs}${tas_in} ${full_tas} # extract out same year of tas data
            fi 

            # import monthly ta data
            ta_in=$(cd ${cwd}/${dirs} && ls --ignore=$ignorestr ta_*${clim}*12.nc)
            begin_in=${ta_in#"${ta_in%_*-*}_"} # extract just the year part of name containing first year
            ta_out=ta_Amon_${dirs%/}_${clim}_r1i1p1_${yr_begin}${mn_begin}-${yr_end}${mn_end}.nc
            full_ta=${cwd}/${dirs}$ta_out

            if [[ ! -f $full_ta ]]; then
                cdo seldate,${yr_begin}-${mn_begin}-01,${yr_end}-${mn_end}-31 ${cwd}/${dirs}${ta_in} ${full_ta} # extract out same year of ta data
            fi 

            # import monthly hus data
            hus_in=$(cd ${cwd}/${dirs} && ls --ignore=$ignorestr hus_*${clim}*12.nc)
            begin_in=${hus_in#"${hus_in%_*-*}_"} # extract just the year part of name containing first year
            hus_out=hus_Amon_${dirs%/}_${clim}_r1i1p1_${yr_begin}${mn_begin}-${yr_end}${mn_end}.nc
            full_hus=${cwd}/${dirs}$hus_out

            if [[ ! -f $full_hus ]]; then
                cdo seldate,${yr_begin}-${mn_begin}-01,${yr_end}-${mn_end}-31 ${cwd}/${dirs}${hus_in} ${full_hus} # extract out same year of hus data
            fi 

            # import monthly zg data
            zg_in=$(cd ${cwd}/${dirs} && ls --ignore=$ignorestr zg_*${clim}*12.nc)
            begin_in=${zg_in#"${zg_in%_*-*}_"} # extract just the year part of name containing first year
            zg_out=zg_Amon_${dirs%/}_${clim}_r1i1p1_${yr_begin}${mn_begin}-${yr_end}${mn_end}.nc
            full_zg=${cwd}/${dirs}$zg_out

            if [[ ! -f $full_zg ]]; then
                cdo seldate,${yr_begin}-${mn_begin}-01,${yr_end}-${mn_end}-31 ${cwd}/${dirs}${zg_in} ${full_zg} # extract out same year of zg data
            fi 
            
            newdir=${cwd}/${dirs}tend
            if [ ! -d "$newdir" ]; then # create directory if it doesn't exist
                mkdir -p $newdir
            fi

            # import monthly huss data if available
            case ${no_huss} in *"${dirs}"*) # if not available then interpolate huss from 3d hus
                full_huss=${full_hus}
                :
                ;;
            *)

                # import monthly huss data
                huss_in=$(cd ${cwd}/${dirs} && ls --ignore=$ignorestr huss_*${clim}*12.nc)
                begin_in=${huss_in#"${huss_in%_*-*}_"} # extract just the year part of name containing first year
                huss_out=huss_Amon_${dirs%/}_${clim}_r1i1p1_${yr_begin}${mn_begin}-${yr_end}${mn_end}.nc
                full_huss=${cwd}/${dirs}$huss_out

                if [[ ! -f $full_huss ]]; then
                    cdo seldate,${yr_begin}-${mn_begin}-01,${yr_end}-${mn_end}-31 ${cwd}/${dirs}${huss_in} ${full_huss} # extract out same year of huss data
                fi 

            esac

            fn_tend="${fn_ta/ta_/tendmon_}"
            full_tend=${newdir}/${fn_tend}
            # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_tend_mon.py ${full_ps} ${full_ta} ${full_zg} ${full_hus} ${full_tend} ${full_tas} ${full_orog} ${full_huss}
            python ${cwd}/make_tend_mon.py ${full_ps} ${full_ta} ${full_zg} ${full_hus} ${full_tend} ${full_tas} ${full_orog} ${full_huss}

        fi
        cd $rwd/$dirs # go in the model directory

        ;;
    esac 

    # take ymonmean
    cdo ymonmean ${full_tend} ${cwd}/${dirs}${fn_tend}.ymonmean.nc

done

cd $cwd
