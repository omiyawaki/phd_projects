#!/bin/sh

declare -a vtypes=("rad" "pe" "stf" "temp") # variable type
yr_start=2000
yr_end=2012

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

    $cwd/concat.sh $vtype $yr_start $yr_end

    cd $cwd

done

