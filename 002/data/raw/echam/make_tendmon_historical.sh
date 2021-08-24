#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python
source activate p003

# declare -a models=("rp000172" "rp000086" "rp000034" "rp000141" "rp000133" "rp000145" "rp000131" "rp000147" "rp000135" "rp000149" "rp000046" "rp000134") # list of GCM models to process
declare -a models=("rp000126" "rp000148" "rp000146" "rp000130" "rp000144" "rp000132" "rp000140" "rp000124") # list of GCM models to process
declare -a clim="historical" # climate name
declare -a freq="mon" # compute mse tendency using daily or monthly data?

cwd=$(pwd) # save current working directory
cd /project2/tas1/ockham/data11/tas/echam-aiv_rcc_6.1.00p1/ # switch to directory with raw data
rwd=$(pwd) # save raw data directory

# for dirs in */; do # loop through all the models
for dirs in ${models[@]}; do # loop through models

        echo $dirs
        mkdir -p $cwd/$dirs # make model directory in processed data folder if it doesn't exist yet
        cd $rwd/$dirs # go in the model directory

        if [[ $freq == "day" ]]; then
            checkstring=$cwd/${dirs}tend_*$clim*39.nc
        elif [[ $freq == "mon" ]]; then
            checkstring=$cwd/${dirs}tendmon_*$clim*39.nc
        fi

        if ls $checkstring 1> /dev/null 2>&1; then # check if data is already there
            echo "tend was already created. Skipping..."
        else

            if [[ $dirs == "rp000172" ]]; then
                for yr in {0009..0014}; do
                    echo ${yr}

                    full_ps=$(pwd)/BOT_${dirs}_${yr}.nc
                    full_tas=$(pwd)/BOT_${dirs}_${yr}.nc
                    full_orog=$(pwd)/BOT_${dirs}_${yr}.nc
                    full_huss=$(pwd)/BOT_${dirs}_${yr}.nc

                    full_ta=$(pwd)/ATM_${dirs}_${yr}.nc
                    full_zg=$(pwd)/ATM_${dirs}_${yr}.nc
                    full_hus=$(pwd)/ATM_${dirs}_${yr}.nc

                    fn_tend=tendmon_${dirs}_${yr}.nc
                    full_tend=${cwd}/${dirs}/${fn_tend}

                    # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_tend_mon.py ${full_ps} ${full_ta} ${full_zg} ${full_hus} ${full_tend} ${full_tas} ${full_orog} ${full_huss}
                    python ${cwd}/make_tend_mon.py ${full_ps} ${full_ta} ${full_zg} ${full_hus} ${full_tend} ${full_tas} ${full_orog} ${full_huss}
                done

                cdo -O mergetime ${cwd}/${dirs}/tendmon*{0009..0014}.nc ${cwd}/${dirs}/tendmon_${dirs}_0009_14.nc
                cdo -O shifttime,-15day ${cwd}/${dirs}/tendmon_${dirs}_0009_14.nc ${cwd}/${dirs}/tendmon_${dirs}_0009_14.shifttime.nc # shift back timestamps by 15 days because February data is reported as March otherwise
                cdo -O ymonmean ${cwd}/${dirs}/tendmon_${dirs}_0009_14.shifttime.nc ${cwd}/${dirs}/tendmon_${dirs}_0009_14.ymonmean.nc
            else
                for yr in {0020..0039}; do
                    echo ${yr}

                    full_ps=$(pwd)/BOT_${dirs}_${yr}.nc
                    full_tas=$(pwd)/BOT_${dirs}_${yr}.nc
                    full_orog=$(pwd)/BOT_${dirs}_${yr}.nc
                    full_huss=$(pwd)/BOT_${dirs}_${yr}.nc

                    full_ta=$(pwd)/ATM_${dirs}_${yr}.nc
                    full_zg=$(pwd)/ATM_${dirs}_${yr}.nc
                    full_hus=$(pwd)/ATM_${dirs}_${yr}.nc

                    fn_tend=tendmon_${dirs}_${yr}.nc
                    full_tend=${cwd}/${dirs}/${fn_tend}

                    # srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_tend_mon.py ${full_ps} ${full_ta} ${full_zg} ${full_hus} ${full_tend} ${full_tas} ${full_orog} ${full_huss}
                    python ${cwd}/make_tend_mon.py ${full_ps} ${full_ta} ${full_zg} ${full_hus} ${full_tend} ${full_tas} ${full_orog} ${full_huss}
                done

                cdo -O mergetime ${cwd}/${dirs}/tendmon*{0020..0039}.nc ${cwd}/${dirs}/tendmon_${dirs}_0020_39.nc
                cdo -O shifttime,-15day ${cwd}/${dirs}/tendmon_${dirs}_0020_39.nc ${cwd}/${dirs}/tendmon_${dirs}_0020_39.shifttime.nc # shift back timestamps by 15 days because February data is reported as March otherwise
                cdo -O ymonmean ${cwd}/${dirs}/tendmon_${dirs}_0020_39.shifttime.nc ${cwd}/${dirs}/tendmon_${dirs}_0020_39.ymonmean.nc
            fi

        fi
        cd $rwd/$dirs # go in the model directory

done

cd $cwd
