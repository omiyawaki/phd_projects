#!/usr/bin/env bash
# set -euo pipefail

cwd=$(pwd)
fdir="${d003}/data/raw/echam"

temp=${fdir}/rp000134/tsurf_rp000134_0020_0039.ymonmean-20.nc
fsfc=${fdir}/rp000134/fsfc_rp000134_0020_0039.ymonmean-20.nc
qflx=${cwd}/qflux_ice40m_a.nc

python make_qflux_a.py $temp $fsfc $qflx
