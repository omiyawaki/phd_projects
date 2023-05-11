#!/usr/bin/env bash
# set -euo pipefail

cwd=$(pwd)
fdir="${d003}/data/raw/echam"

temp=${fdir}/rp000134/tsurf_rp000134_0020_0039.ymonmean-20.nc
fsfc=${fdir}/rp000134/fsfc_rp000134_0020_0039.ymonmean-20.nc
swref=${fdir}/rp000134/swref_rp000134_0020_0039.ymonmean-20.nc
albedo=${fdir}/rp000134/albedo_rp000134_0020_0039.ymonmean-20.nc
sradsd=${fdir}/rp000134/sradsd_rp000134_0020_0039.ymonmean-20.nc
qflx=${cwd}/qflux_ice40m_e.nc
qsw=${cwd}/qsw_ice40m_e.nc

python make_qflux_e.py $temp $fsfc $swref $albedo $sradsd $qflx $qsw
