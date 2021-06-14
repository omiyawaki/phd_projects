#!/bin/sh

module unload python
module load python

# declare -a vtypes=("rad" "hydro" "stf" "temp" "temp_ml" "srfc" "zg" "w500" "tend" "div" "vas" "albedo") # variable type
declare -a vtypes=("rad" "hydro" "stf" "temp" "temp_ml" "srfc" "zg") # variable type
#declare -a vtypes=("hus") # variable type
yr_start=1979
mon_start=1
day_start=1
yr_end=2005
mon_end=12
day_end=31

cwd=$(pwd)

for vtype in ${vtypes[@]}; do
    cd $vtype

    for yr in $(seq $yr_start $yr_end); do
        file=interim_${vtype}_${yr}
        if test -f "$file.nc"; then
            echo "${file}.nc already exists. Skipping..."
        else
            python get_interim_${vtype}.py $yr
        fi
    done

    $cwd/concat.sh $vtype $yr_start $yr_end $mon_start $mon_end $day_start $day_end

    cd $cwd

done

