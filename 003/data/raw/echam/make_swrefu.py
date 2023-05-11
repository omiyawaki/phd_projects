import sys
import numpy as np
from scipy import interpolate, integrate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps and mse files
path_sradsd = sys.argv[1]
path_srafsd = sys.argv[2]
path_swrefu = sys.argv[3]

# open files
file_sradsd = Dataset(path_sradsd, 'r')
file_srafsd = Dataset(path_srafsd, 'r')

# read data
sradsd = file_sradsd.variables['sradsd'][:] # (mon x lat x lon)
srafsd = file_srafsd.variables['srafsd'][:] # (mon x lat x lon)

# compute atmos reflectivity to upwelling shortwave flux
swref_up = 0.05+0.85*(1-sradsd/srafsd)

# save file as netCDF
file_swrefu = Dataset(path_swrefu, "w", format='NETCDF4_CLASSIC')

# copy attributes from sradsd file
file_swrefu.setncatts(file_sradsd.__dict__)

# copy dimensions from sradsd file
for name, dimension in file_sradsd.dimensions.items():
    file_swrefu.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from sradsd file
for name, variable in file_sradsd.variables.items():
    if any(name in s for s in ['sradsd']):
        continue
    
    x = file_swrefu.createVariable(name, variable.datatype, variable.dimensions)
    file_swrefu[name].setncatts(file_sradsd[name].__dict__)
    file_swrefu[name][:] = file_sradsd[name][:]
    
swrefu = file_swrefu.createVariable('swrefu', 'f4', ("time","lat","lon"))
swrefu.units = "unitless"
swrefu.long_name = "atmospheric shortwave reflectivity (upwelling)"
swrefu[:] = swref_up

file_swrefu.close()
