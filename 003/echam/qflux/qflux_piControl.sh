#!/usr/bin/env bash
# set -euo pipefail

cwd=$(pwd)
fdir="${d003}/data/raw/echam"

infile=${cwd}/qflux_zero.nc
inpic=${cwd}/q_flux_final_c_Amon_MPI-ESM-LR_piControl_r1i1p1_280001-284912.timmean.nc
oqflux=${cwd}/qflux_piControl.nc

python make_qflux_pic.py $infile $inpic $oqflux 
