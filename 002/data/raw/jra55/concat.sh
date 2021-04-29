#!/bin/sh

declare -a vtypes=("rad" "hydro" "stf" "temp" "temp_trop" "srfc" "zg" "zg_trop" "hus" "hus_trop" "lfrac" "rh") # variable type
#declare -a vtypes=("temp_trop" "zg_trop")
yr_start="1980"
mon_start="1"
day_start="1"
yr_end="2005"
mon_end="12"
day_end="31"

cwd=$(pwd)

for vtype in ${vtypes[@]}; do
    cd $vtype

    if [ "$vtype" == "rad" ]; then
        declare -a vnames=("dswrf" "dlwrf" "uswrf" "ulwrf") # variable names of this type
    elif [ "$vtype" == "radcs" ]; then
        declare -a vnames=("csusf" "csdsf" "csulf" "csdlf") # variable names of this type
    elif [ "$vtype" == "hydro" ]; then
        declare -a vnames=("evp" "tprat" "cprat")
    elif [ "$vtype" == "stf" ]; then
        declare -a vnames=("lhtfl" "shtfl")
    elif [ "$vtype" == "temp" ]; then
        declare -a vnames=("tmp")
    elif [ "$vtype" == "temp_trop" ]; then
        declare -a vnames=("tmp")
    elif [ "$vtype" == "temp_ml" ]; then
        declare -a vnames=("tmp")
    elif [ "$vtype" == "srfc" ]; then
        declare -a vnames=("pres" "tmp" "rh")
    elif [ "$vtype" == "zg" ]; then
        declare -a vnames=("hgt")
    elif [ "$vtype" == "zg_trop" ]; then
        declare -a vnames=("hgt")
    elif [ "$vtype" == "lfrac" ]; then
        declare -a vnames=("land")
    elif [ "$vtype" == "sice" ]; then
        declare -a vnames=("icec")
    elif [ "$vtype" == "rh" ]; then
        declare -a vnames=("rh")
    elif [ "$vtype" == "hus" ]; then
        declare -a vnames=("spfh")
    elif [ "$vtype" == "hus_trop" ]; then
        declare -a vnames=("spfh")
    fi

    for vname in ${vnames[@]}; do
        filename=jra55_${vname}_${yr_start}_${yr_end}

#        if test -f "$filename.ymonmean.nc"; then
#            echo "$filename.ymonmean.nc exists. Skipping..."
#        else
            cdo -O -b F64 mergetime *${vname}*1979_2005.nc $filename.merge.nc
            cdo -O seldate,$yr_start-$mon_start-$day_start,$yr_end-$mon_end-$day_end $filename.merge.nc $filename.nc
            rm $filename.merge.nc

            cdo -O ymonmean $filename.nc $filename.ymonmean.nc
#        fi
    done

    cd $cwd

done
