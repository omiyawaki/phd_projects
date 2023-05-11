#!/bin/sh

in="ghg_hist+ssp585_0000-2500.nc"
out="ghg_hist+ssp585_1765-2500.nc"

cdo seldate,1765-01-01,2500-12-31 $in $out
