#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

cwd=$(pwd)
infile=${cwd}/greenhouse_aqua_timeind.nc
outfile=${cwd}/greenhouse_abrupt16xCO2.nc

python ${cwd}/make_abrupt16x.py $infile $outfile
