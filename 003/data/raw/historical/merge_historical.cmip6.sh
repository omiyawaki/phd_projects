#!/usr/bin/env bash
set -euo pipefail

declare -a realm=("atmos")
declare -a vars_gcm=("ta") # list of GCM variables that we want to process
# declare -a vars_gcm=("ta" "zg" "hus") # list of GCM variables that we want to process
# declare -a vars_gcm=("huss" "hurs" "ps" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "prc" "evspsbl" "rlutcs" "rldscs" "rsuscs" "rsdscs" "rsutcs") # list of GCM variables that we want to process

# declare -a vars_gcm=("huss" "hurs" "hur" "hus" "zg" "ta" "ps" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "prc" "evspsbl" "rlutcs" "rldscs" "rsuscs" "rsdscs" "rsutcs") # list of GCM variables that we want to process
# declare -a vars_gcm=("huss" "hurs" "hur" "hus" "zg" "ta" "ps" "ts" "tas" "pr" "prc" "evspsbl") # list of GCM variables that we want to process
# declare -a vars_gcm=("rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "rlutcs" "rldscs" "rsuscs" "rsdscs" "rsutcs") # list of GCM variables that we want to process

# declare -a vars_gcm=("zg" "ta" "hur" "ps" "ts" "tas" "rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss" "pr" "prc" "evspsbl") # list of GCM variables that we want to process
# declare -a vars_gcm=("rlut" "rsut" "rsdt" "rlus" "rlds" "rsds" "rsus" "hfls" "hfss") # list of GCM variables that we want to process
# declare -a vars_gcm=("rlutcs" "rsutcs" "rldscs" "rsdscs" "rsuscs") # list of GCM variables that we want to process
# declare -a vars_gcm=("ps" "tas" "ta" "zg" "hus") # list of GCM variables that we want to process
# declare -a vars_gcm=("ta" "hus" "hur") # list of GCM variables that we want to process

# declare -a realm=("seaIce")
# declare -a vars_gcm=("siconc") # list of GCM variables that we want to process

declare -a clim="historical" # climate name
declare -a freq="mon" # data output frequency (e.g. fx for fixed, mon for monthly, day for daily)

declare -a ens="r1i1p1f1" # ensemble specification 
declare -a models=("MRI-ESM2-0/") # extended RCP runs
# declare -a models=("ACCESS-CM2/" "ACCESS-ESM1-5/" "CanESM5/" "IPSL-CM6A-LR/" "MRI-ESM2-0/" "CESM2-WACCM/") # extended RCP runs

# declare -a ens="r1i1p1f2" # ensemble specification 
# declare -a models=("MIROC-ES2L/") # extended RCP runs

# declare -a ens="r3i1p1f2" # ensemble specification 
# declare -a models=("GISS-E2-1-H/" "GISS-E2-1-G/") # extended RCP runs

# declare -a ens="r4i1p1f2" # ensemble specification 
# declare -a models=("UKESM1-0-LL/") # extended RCP runs

mean=""
# declare -a ens="r0i0p0" # ensemble specification 
# declare -a skip_files=("_eady.nc")
declare -a skip_files=("185001-201412.nc 198001-200012.nc 198001-200012_ave.nc")
# declare -a skip_files=("198001-200012.nc 198001-200012_ave.nc")

declare -a skip_models="inmcm4/ FGOALS-s2/ CanAM4/ CanCM4/ CESM1-CAM5-1-FV2/ EC-EARTH/ FIO-ESM/ HadGEM2-A/ HadGEM2-AO/ GFDL-HIRAM-C180/ GFDL-HIRAM-C360/ MRI-AGCM3-2H/ MRI-AGCM3-2S/ NICAM-09/ MPI-ESM-P/"

out_yr_begin=1850
out_mn_begin=01
out_yr_end=2014
out_mn_end=12

cwd=$(pwd) # save current working directory
cd /project2/tas1/ockham/data9/tas/CMIP6_RAW # switch to directory with raw data
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
                pattern="${vars}_*${dirs%/}*${mean}.nc"
                files=( $pattern )
                # remove files that are not original (list of non-original data is contained in variable $delete )
                for (( i=0; i<${#files[@]}; i++ )); do 
                    # case $skip_files in *"${files[i]#*_}"* | *"${files[i]:(-16)}"* )
                    # case $skip_files in *"${files[i]#*_}"* | *"${files[i]:(-8)}"* | *"${files[i]:(-16)}"* )
                    case $skip_files in *"${files[i]#*_}"* | *"${files[i]:(-16)}"* | *"${files[i]:(-20)}"* )
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
                    yr_end=${yr: 7:4} # extract the last year and month
                    sel="${common}${out_yr_begin}${out_mn_begin}-${out_yr_end}${out_mn_end}" # write files name with last 5 years
                    cp ${files[@]} $cwd/$dirs
                    # cdo -O yearmean $cwd/$dirs${files[@]} $cwd/$dirs${files[@]%.nc}.yearmean.nc # take annual mean
                    merge=${files[@]%.nc}
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
                    # cdo -O yearmean $cwd/$dirs$merge.nc $cwd/$dirs$merge.yearmean.nc # combine multiple files into one
                fi
                # fi

                # if [ ! $freq == "fx" ]; then
                #     if [ ! $yr_end -eq 2005 ]; then
                #         selmerge="${common}${out_yr_begin}${out_mn_begin}-${out_yr_end}${out_mn_end}" # write file name with merged time
                #         cdo -O seldate,${out_yr_begin}-${out_mn_begin}-01,${out_yr_end}-${out_mn_end}-31 $cwd/$dirs$merge.nc $cwd/$dirs$selmerge.nc
                #         cdo -O yearmean $cwd/$dirs$selmerge.nc $cwd/$dirs$selmerge.yearmean.nc # combine multiple files into one
                #     fi
                # fi

                #######################################################################
                # convert curvilinear to standard lat lon grid for sea ice data
                #######################################################################
                if [ ${vars} == "siconca" ]; then
                    sica_file=$(ls ${cwd}/${dirs}siconca_*${out_yr_begin}${out_mn_begin}-${out_yr_end}${out_mn_end}${mean}.nc)
                    sic_file="${sica_file//siconca/siconc}"
                    sic_file_rn=${sic_file%.nc}.rn.nc
                    mv ${sica_file} ${sic_file_rn}
                    cdo chname,siconca,siconc ${sic_file_rn} ${sic_file}
                    rm ${sic_file_rn}
                fi

                if [ ${vars} == "siconc" ]; then
                    ref_file=$(ls ${cwd}/${dirs}tas_*${out_yr_begin}${out_mn_begin}-${out_yr_end}${out_mn_end}${mean}.nc)
                    sic_file=$(ls ${cwd}/${dirs}siconc_*${out_yr_begin}${out_mn_begin}-${out_yr_end}${out_mn_end}${mean}.nc)

                    # extract siconc only for IPSL
                    if [ "${dirs%/}" == "IPSL-CM6A-LR" ]; then
                        echo "Extracting siconc only..."
                        cdo -selvar,siconc ${sic_file} ${sic_file%.nc}.origgrid.nc
                    # sic in canesm5 is not in percentage so correct by multiplying by factor of 100
                    elif [ "${dirs%/}" == "CanESM5" ]; then
                        cdo -mulc,100 ${sic_file} ${sic_file%.nc}.origgrid.nc
                    else
                        mv ${sic_file} ${sic_file%.nc}.origgrid.nc
                    fi


                    # first create file containing standard lat-lon grid data (e.g., using tas file)
                    cdo griddes ${ref_file} > ${cwd}/${dirs}grid_latlon
                    sed -i "s/generic/lonlat/g" ${cwd}/${dirs}grid_latlon

                    # convert to lat lon
                    cdo -remapbil,${cwd}/${dirs}grid_latlon ${sic_file%.nc}.origgrid.nc ${sic_file}
                    
                fi

                if [ ${vars} == "sithick" ]; then
                    ref_file=$(ls ${cwd}/${dirs}tas_*${out_yr_begin}${out_mn_begin}-${out_yr_end}${out_mn_end}${mean}.nc)
                    sic_file=$(ls ${cwd}/${dirs}sithick_*${out_yr_begin}${out_mn_begin}-${out_yr_end}${out_mn_end}${mean}.nc)

                    # extract sithick only for IPSL
                    if [ "${dirs%/}" == "IPSL-CM6A-LR" ]; then
                        echo "Extracting sithick only..."
                        cdo -selvar,sithick ${sic_file} ${sic_file%.nc}.origgrid.nc
                    # sic in canesm5 is not in percentage so correct by multiplying by factor of 100
                    elif [ "${dirs%/}" == "CanESM5" ]; then
                        cdo -mulc,100 ${sic_file} ${sic_file%.nc}.origgrid.nc
                    else
                        mv ${sic_file} ${sic_file%.nc}.origgrid.nc
                    fi


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
