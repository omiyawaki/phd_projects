#!/usr/bin/env bash
# set -euo pipefail

cwd=$(pwd)
fdir="${d003}/data/raw/echam"

infile=${cwd}/qflux_zero.nc

dtt=${fdir}/rp000134/tsurf_rp000134_0020_0039.ymonmean-20.nc
swn=${fdir}/rp000134/srads_rp000134_0020_0039.ymonmean-20.nc
swu=${fdir}/rp000134/sradsu_rp000134_0020_0039.ymonmean-20.nc
lw=${fdir}/rp000134/trads_rp000134_0020_0039.ymonmean-20.nc
lh=${fdir}/rp000134/ahfl_rp000134_0020_0039.ymonmean-20.nc
sh=${fdir}/rp000134/ahfs_rp000134_0020_0039.ymonmean-20.nc

oqflux=${cwd}/qflux_ice40m.nc

rm $oqflux
python make_qflux_malte.py $infile $dtt $swn $swu $lw $lh $sh $oqflux 
