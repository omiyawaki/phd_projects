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
path_albedo = sys.argv[3]
path_swref0 = sys.argv[4]
path_swabs0 = sys.argv[5]

# open files
file_sradsd = Dataset(path_sradsd, 'r')
file_srad0d = Dataset(path_srad0d, 'r')
file_albedo = Dataset(path_albedo, 'r')
file_swref0 = Dataset(path_swref0, 'r')

# read data
sradsd = file_sradsd.variables['sradsd'][:] # (mon x lat x lon)
srad0d = file_srad0d.variables['srad0d'][:] # (mon x lat x lon)
albedo = file_albedo.variables['albedo'][:] # (mon x lat x lon)
swref0 = file_swref0.variables['swref0'][:] # (mon x lat x lon)

# compute atmospheric sw absorptivity following Donohoe and Battisti (2010)
ref = 1-swref0-(sradsd*(1-albedo*swref0))/srad0d

# save file as netCDF
file_swabs0 = Dataset(path_swabs0, "w", format='NETCDF4_CLASSIC')

# copy attributes from sradsd file
file_swabs0.setncatts(file_sradsd.__dict__)

# copy dimensions from sradsd file
for name, dimension in file_sradsd.dimensions.items():
    file_swabs0.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from sradsd file
for name, variable in file_sradsd.variables.items():
    if any(name in s for s in ['sradsd']):
        continue
    
    x = file_swabs0.createVariable(name, variable.datatype, variable.dimensions)
    file_swabs0[name].setncatts(file_sradsd[name].__dict__)
    file_swabs0[name][:] = file_sradsd[name][:]
    
swabs0 = file_swabs0.createVariable('swabs0', 'f4', ("time","lat","lon"))
swabs0.units = "unitless"
swabs0.long_name = "atmospheric sw absorptivity"
swabs0[:] = ref

file_swabs0.close()
