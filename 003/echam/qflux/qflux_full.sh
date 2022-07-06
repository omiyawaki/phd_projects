#!/usr/bin/env bash
# set -euo pipefail

cwd=$(pwd)
fdir="${d003}/data/raw/echam"

infile=${cwd}/qflux_zero.nc

alice=${fdir}/rp000134/albedo_rp000134_0020_0039.ymonmean-20.nc
alnoice=${fdir}/rp000135/albedo_rp000135_0020_0039.ymonmean-20.nc

sic=${fdir}/rp000134/seaice_rp000134_0020_0039.ymonmean-20.nc
sit=${fdir}/rp000134/siced_rp000134_0020_0039.ymonmean-20.nc
swn=${fdir}/rp000134/srads_rp000134_0020_0039.ymonmean-20.nc
swu=${fdir}/rp000134/sradsu_rp000134_0020_0039.ymonmean-20.nc
lw=${fdir}/rp000134/trads_rp000134_0020_0039.ymonmean-20.nc
lh=${fdir}/rp000134/ahfl_rp000134_0020_0039.ymonmean-20.nc
sh=${fdir}/rp000134/ahfs_rp000134_0020_0039.ymonmean-20.nc

oqflux=${cwd}/qflux_ice40m.nc

rm $oqflux
python make_qflux_full.py $infile $alice $alnoice $sic $sit $swn $swu $lw $lh $sh $oqflux 
