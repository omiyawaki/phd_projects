#!/bin/sh

in1="ghg_hist6_0000-2014.nc"
in2="ghg_ssp585_2015-2500.nc"
out="ghg_hist+ssp585_0000-2500.nc"

cdo mergetime $in1 $in2 $out
