#!/bin/sh

# declare -a forcings=("0")
# declare -a forcings=("0" "10")
declare -a forcings=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10")
# declare -a forcings=("-1" "-2" "-3" "-4" "-5" "-6" "-7" "-8" "-9" "-10")
# declare -a forcings=("-10")

cwd=$(pwd)

for f in ${forcings[@]}; do

    echo -e "\nSolving $f W m**-2 Forcing...\n"
    folder="$cwd/frad_$f"
    mkdir -p $folder

    ################################################################################
    # SIMULATION
    ################################################################################
    # python we15.py $f $folder

    ################################################################################
    # DATA PROCESSING
    ################################################################################
    # python yearmean.py $f $folder
    python djfmean.py $f $folder

    ################################################################################
    # PLOTS
    ################################################################################
    # python f2_we15.py $f $folder

    # python equil.py $f $folder

    echo -e "\nDone.\n"
done
