#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

cwd=$(pwd)
ircp45=${cwd}/greenhouse_rcp45.nc
ircp85=${cwd}/ghg_rcp85_1765-2500_c100203.nc
orcp85=${cwd}/greenhouse_rcp85.nc

rm $orcp85
ncl i45=\"$ircp45\" i85=\"$ircp85\" o85=\"$orcp85\" ${cwd}/make_rcp85.ncl 
