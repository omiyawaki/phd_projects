#!/usr/bin/env bash
# set -euo pipefail

cwd=$(pwd)
fdir="${d003}/data/raw/echam"

temp=${fdir}/rp000134/tsurf_rp000134_0020_0039.ymonmean-20.nc
fsfc=${fdir}/rp000134/fsfc_rp000134_0020_0039.ymonmean-20.nc
swrefu=${fdir}/rp000135/swrefu_rp000135_0020_0039.ymonmean-20.nc
albedo=${fdir}/rp000134/albedo_rp000134_0020_0039.ymonmean-20.nc
sradsd=${fdir}/rp000134/sradsd_rp000134_0020_0039.ymonmean-20.nc
qflx=${cwd}/qflux_ice40m_h.nc
qsw=${cwd}/qsw_ice40m_h.nc

python make_qflux_f.py $temp $fsfc $swrefu $albedo $sradsd $qflx $qsw
