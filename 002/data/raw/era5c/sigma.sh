#!/usr/bin/env bash
# set -euo pipefail

sim="era5c"
yr_begin="1980"
yr_end="2005"

cwd=$(pwd) # save current working directory

full_srfc=${cwd}/srfc/${sim}_srfc_${yr_begin}_${yr_end}.nc
full_temp=${cwd}/temp/${sim}_temp_${yr_begin}_${yr_end}.nc
full_tempsi=${cwd}/temp/${sim}_tempsi_${yr_begin}_${yr_end}.nc

# srun --partition=tas1 --time=6:00:00 --exclusive --pty python ${cwd}/make_sigma.py ${full_srfc} ${full_temp} ${full_tempsi}
python ${cwd}/make_sigma.py ${full_srfc} ${full_temp} ${full_tempsi}

