#!/bin/bash

# sequence of scripts to create diffusivity-related variables 

# select 925 hPa mse
bash sellev.sh # make sure variable is set to mse
bash ymonmean.sh # take ymonmean 
bash msez.sh # take zonal mean
