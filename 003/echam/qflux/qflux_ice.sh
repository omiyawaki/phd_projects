#!/usr/bin/env bash
# set -euo pipefail

cwd=$(pwd)
fdir="${d003}/data/raw/echam"

infile=${cwd}/qflux_zero.nc
sicfile=${fdir}/rp000134/seaice_rp000134_0020_0039.ymonmean-20.nc
sitfile=${fdir}/rp000134/siced_rp000134_0020_0039.ymonmean-20.nc
inice=${fdir}/rp000134/fsfc_rp000134_0020_0039.ymonmean-20.nc
innoice=${fdir}/rp000135/fsfc_rp000135_0020_0039.ymonmean-20.nc
oqflux=${cwd}/qflux_ice40m.old.nc

rm $oqflux
python make_qflux_ice.py $infile $sicfile $sitfile $inice $innoice $oqflux 
