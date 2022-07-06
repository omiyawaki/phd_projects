#!/usr/bin/env bash
# set -euo pipefail

cwd=$(pwd)
fdir="${d003}/data/raw/echam"

infile=${cwd}/qflux_zero.nc
ina32=${cwd}/q_flux_final_32_mon_MPIESM12_abrupt32x_103.timmean.nc
oqflux=${cwd}/qflux_abrupt32x.nc

python make_qflux_a32.py $infile $ina32 $oqflux 
