#!/usr/bin/env bash
# set -euo pipefail

# load python environment
module unload python
module load python

cwd=$(pwd)
infile=${cwd}/greenhouse_aqua_timeind.nc
outfile=${cwd}/greenhouse_piControl.nc

python ${cwd}/make_piControl.py $infile $outfile
