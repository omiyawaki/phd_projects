#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python
source activate netcdf

declare -a models=("HadCM3/") # list of GCM models to process
declare -a clim="historical" # climate name
declare -a freq="mon" # compute mse tendency using daily or monthly data?

yr_begin="1979" # extract the first year and month
mn_begin="01" # extract the first year and month
yr_end="2005" # extract the last year and month
mn_end="12" # extract the last year and month

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

            if [ ${#files[@]} -eq 0 ]; then # check if this file exists
                echo "File of type $pattern does not exist. Please download the file and place it in the corresponding directory."
                exit
            fi
            
            fn_ta=ta_Amon_${dirs%/}_${clim}_r1i1p1_${yr_begin}${mn_begin}-${yr_end}${mn_end}.nc
            
            prepath_hur="${prepath_ta/\/ta\///hur/}"
            prepath_zg="${prepath_ta/\/ta\///zg/}"
            
            fn_hur="${fn_ta/ta_/hur_}"
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

            # import monthly hur data
            hur_in=$(cd ${cwd}/${dirs} && ls --ignore=$ignorestr hur_*${clim}*12.nc)
            begin_in=${hur_in#"${hur_in%_*-*}_"} # extract just the year part of name containing first year
            hur_out=hur_Amon_${dirs%/}_${clim}_r1i1p1_${yr_begin}${mn_begin}-${yr_end}${mn_end}.nc
            full_hur=${cwd}/${dirs}$hur_out

            if [[ ! -f $full_hur ]]; then
                cdo seldate,${yr_begin}-${mn_begin}-01,${yr_end}-${mn_end}-31 ${cwd}/${dirs}${hur_in} ${full_hur} # extract out same year of hur data
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

            # import monthly hurs data
            hurs_in=$(cd ${cwd}/${dirs} && ls --ignore=$ignorestr hurs_*${clim}*12.nc)
            begin_in=${hurs_in#"${hurs_in%_*-*}_"} # extract just the year part of name containing first year
            hurs_out=hurs_Amon_${dirs%/}_${clim}_r1i1p1_${yr_begin}${mn_begin}-${yr_end}${mn_end}.nc
            full_hurs=${cwd}/${dirs}$hurs_out

            if [[ ! -f $full_hurs ]]; then
                cdo seldate,${yr_begin}-${mn_begin}-01,${yr_end}-${mn_end}-31 ${cwd}/${dirs}${hurs_in} ${full_hurs} # extract out same year of hurs data
            fi 

            fn_tend="${fn_ta/ta_/tendmon_}"
            full_tend=${newdir}/${fn_tend}
            # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_tend_mon.py ${full_ps} ${full_ta} ${full_zg} ${full_hur} ${full_tend} ${full_tas} ${full_orog} ${full_hurs}
            python ${cwd}/make_tend_mon_hur.py ${full_ps} ${full_ta} ${full_zg} ${full_hur} ${full_tend} ${full_tas} ${full_orog} ${full_hurs}

        fi
        cd $rwd/$dirs # go in the model directory

        ;;
    esac 

done

cd $cwd
