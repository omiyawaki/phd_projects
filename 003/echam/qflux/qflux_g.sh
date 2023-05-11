#!/usr/bin/env bash
# set -euo pipefail

cwd=$(pwd)
fdir="${d003}/data/raw/echam"

temp=${fdir}/rp000134/tsurf_rp000134_0020_0039.ymonmean-20.nc
fsfc=${fdir}/rp000134/fsfc_rp000134_0020_0039.ymonmean-20.nc
swtau=${fdir}/rp000134/swtau_rp000134_0020_0039.ymonmean-20.nc
swtauni=${fdir}/rp000135/swtau_rp000135_0020_0039.ymonmean-20.nc
albedo=${fdir}/rp000134/albedo_rp000134_0020_0039.ymonmean-20.nc
srad0d=${fdir}/rp000134/srad0d_rp000134_0020_0039.ymonmean-20.nc
qflx=${cwd}/qflux_ice40m_g.nc
qsw=${cwd}/qsw_ice40m_g.nc

python make_qflux_g.py $temp $fsfc $swtau $swtauni $albedo $srad0d $qflx $qsw
