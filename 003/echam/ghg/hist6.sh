#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

cwd=$(pwd)
ircp45=${cwd}/greenhouse_rcp45.nc
ircp85=${cwd}/mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-2-0_gr1-GMNHSH_0000-2014.nc
orcp85=${cwd}/greenhouse_hist6.nc

rm $orcp85
ncl i45=\"$ircp45\" i85=\"$ircp85\" o85=\"$orcp85\" ${cwd}/make_hist6.ncl 
