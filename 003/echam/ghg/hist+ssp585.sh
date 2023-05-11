#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

cwd=$(pwd)
ircp45=${cwd}/greenhouse_rcp45.nc
ircp85=${cwd}/ghg_hist+ssp585_1765-2500.nc
orcp85=${cwd}/greenhouse_hist+ssp585.nc

rm $orcp85
ncl i45=\"$ircp45\" i85=\"$ircp85\" o85=\"$orcp85\" ${cwd}/make_hist+ssp585.ncl 
