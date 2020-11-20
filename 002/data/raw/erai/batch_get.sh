#!/bin/sh

module unload python
module load python

declare -a vtypes=("rad" "hydro" "stf" "temp" "temp_ml" "srfc" "zg" "w500" "tend" "div" "vas" "albedo") # variable type
yr_start=2000
mon_start=3
day_start=1
yr_end=2018
mon_end=2
day_end=1

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

