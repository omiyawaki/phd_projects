#!/usr/bin/env bash
# set -euo pipefail

cwd=$(pwd)
fdir="${d003}/data/raw/echam"

sicfile=${fdir}/rp000134/seaice_rp000134_0020_0039.ymonmean-20.nc
yicefile=${fdir}/rp000134/fsfc_rp000134_0020_0039.ymonmean-20.nc
nicefile=${fdir}/rp000135/fsfc_rp000135_0020_0039.ymonmean-20.nc
qflxfile=${cwd}/qflux_ice40m.nc

python make_amp_simple.py $sicfile $yicefile $nicefile $qflxfile
