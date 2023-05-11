#!/usr/bin/env bash
# set -euo pipefail

cwd=$(pwd)
fdir="${d003}/data/raw/echam"

temp=${fdir}/rp000134/tsurf_rp000134_0020_0039.ymonmean-20.nc
fsfc=${fdir}/rp000134/fsfc_rp000134_0020_0039.ymonmean-20.nc
srads=${fdir}/rp000134/srads_rp000134_0020_0039.ymonmean-20.nc
sradsu=${fdir}/rp000134/sradsu_rp000134_0020_0039.ymonmean-20.nc
albedo=${fdir}/rp000134/albedo_rp000134_0020_0039.ymonmean-20.nc
qflx=${cwd}/qflux_ice40m_d.nc
qsw=${cwd}/qsw_ice40m_d.nc

python make_qflux_d.py $temp $fsfc $srads $sradsu $albedo $qflx $qsw
