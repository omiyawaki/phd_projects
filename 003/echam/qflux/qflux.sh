#!/usr/bin/env bash
# set -euo pipefail

cwd=$(pwd)
fdir="${d003}/data/raw/echam"

temp=${fdir}/rp000134/tsurf_rp000134_0020_0039.ymonmean-20.nc
fsfc=${fdir}/rp000134/fsfc_rp000134_0020_0039.ymonmean-20.nc
swni=${fdir}/rp000134/srads_rp000134_0020_0039.ymonmean-20.nc
swnq=${fdir}/rp000190a/srads_rp000190a_0020_0039.ymonmean-20.nc

oqflux=${cwd}/qflux_ice40m.nc

python make_qflux.py $temp $fsfc $swni $swnq $oqflux
