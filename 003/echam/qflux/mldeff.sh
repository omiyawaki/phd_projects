#!/usr/bin/env bash
# set -euo pipefail

cwd=$(pwd)
fdir="${d003}/data/raw/echam"

sicfile=${fdir}/rp000134/seaice_rp000134_0020_0039.ymonmean-20.nc
sitfile=${fdir}/rp000134/siced_rp000134_0020_0039.ymonmean-20.nc
meltfile=${fdir}/rp000134/ahfres_rp000134_0020_0039.ymonmean-20.nc
yicefile=${fdir}/rp000134/fsfc_rp000134_0020_0039.ymonmean-20.nc
nicefile=${fdir}/rp000135/fsfc_rp000135_0020_0039.ymonmean-20.nc
qflxfile=${cwd}/qflux_ice40m.nc
defffile=${cwd}/deff_ice40m.nc
ampffile=${cwd}/ampf_ice40m.nc

python make_mldeff.py $sicfile $sitfile $meltfile $yicefile $nicefile $qflxfile $defffile $ampffile
