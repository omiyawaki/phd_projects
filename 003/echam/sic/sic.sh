#!/usr/bin/env bash
# set -euo pipefail

cwd=$(pwd)
fdir="${d003}/data/raw/echam"

infile=${cwd}/T63_sic_aqua.nc
inice=${fdir}/rp000134/seaice_rp000134_0020_0039.ymonmean-20.nc
ofile=${cwd}/T63_sic_ice40m.nc

python make_sic.py $infile $inice $ofile
