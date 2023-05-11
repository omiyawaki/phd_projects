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
path_sradsu = sys.argv[2]
path_srad0d = sys.argv[3]
path_srad0u = sys.argv[4]
path_alb = sys.argv[5]
path_swref = sys.argv[6]

# open files
file_sradsd = Dataset(path_sradsd, 'r')
file_sradsu = Dataset(path_sradsu, 'r')
file_srad0d = Dataset(path_srad0d, 'r')
file_srad0u = Dataset(path_srad0u, 'r')
file_alb = Dataset(path_alb, 'r')

# read data
sradsd = file_sradsd.variables['sradsd'][:] # (mon x lat x lon)
sradsu = file_sradsu.variables['sradsu'][:] # (mon x lat x lon)
srad0d = file_srad0d.variables['srad0d'][:] # (mon x lat x lon)
srad0u = file_srad0u.variables['srad0u'][:] # (mon x lat x lon)
alb = file_alb.variables['alb'][:] # (mon x lat x lon)

# compute atmospheric sw reflectivity following Donohoe and Battisti (2010)
ref = (-srad0d*srad0u-alb*sradsd**2)/(srad0d**2-alb**2*sradsd**2)
# ref = (-srad0d*srad0u+sradsd*sradsu)/(srad0d**2-sradsu**2)

# save file as netCDF
file_swref = Dataset(path_swref, "w", format='NETCDF4_CLASSIC')

# copy attributes from sradsd file
file_swref.setncatts(file_sradsd.__dict__)

# copy dimensions from sradsd file
for name, dimension in file_sradsd.dimensions.items():
    file_swref.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from sradsd file
for name, variable in file_sradsd.variables.items():
    if any(name in s for s in ['sradsd']):
        continue
    
    x = file_swref.createVariable(name, variable.datatype, variable.dimensions)
    file_swref[name].setncatts(file_sradsd[name].__dict__)
    file_swref[name][:] = file_sradsd[name][:]
    
swref = file_swref.createVariable('swref', 'f4', ("time","lat","lon"))
swref.units = "unitless"
swref.long_name = "atmospheric sw reflectivity"
swref[:] = ref

file_swref.close()
