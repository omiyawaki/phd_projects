#!/bin/sh

cwd=$(pwd)

folder="$cwd"
mkdir -p $folder

# AGGREGATED PLOTS
# python f3_we15.py $folder
python fluxes.py $folder
