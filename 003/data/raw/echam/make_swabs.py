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
path_srad0d = sys.argv[2]
path_alb = sys.argv[3]
path_swref = sys.argv[4]
path_swabs = sys.argv[5]

# open files
file_sradsd = Dataset(path_sradsd, 'r')
file_srad0d = Dataset(path_srad0d, 'r')
file_alb = Dataset(path_alb, 'r')
file_swref = Dataset(path_swref, 'r')

# read data
sradsd = file_sradsd.variables['sradsd'][:] # (mon x lat x lon)
srad0d = file_srad0d.variables['srad0d'][:] # (mon x lat x lon)
alb = file_alb.variables['alb'][:] # (mon x lat x lon)
swref = file_swref.variables['swref'][:] # (mon x lat x lon)

# compute atmospheric sw absorptivity following Donohoe and Battisti (2010)
ref = 1-swref-(sradsd*(1-alb*swref))/srad0d

# save file as netCDF
file_swabs = Dataset(path_swabs, "w", format='NETCDF4_CLASSIC')

# copy attributes from sradsd file
file_swabs.setncatts(file_sradsd.__dict__)

# copy dimensions from sradsd file
for name, dimension in file_sradsd.dimensions.items():
    file_swabs.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from sradsd file
for name, variable in file_sradsd.variables.items():
    if any(name in s for s in ['sradsd']):
        continue
    
    x = file_swabs.createVariable(name, variable.datatype, variable.dimensions)
    file_swabs[name].setncatts(file_sradsd[name].__dict__)
    file_swabs[name][:] = file_sradsd[name][:]
    
swabs = file_swabs.createVariable('swabs', 'f4', ("time","lat","lon"))
swabs.units = "unitless"
swabs.long_name = "atmospheric sw absorptivity"
swabs[:] = ref

file_swabs.close()
