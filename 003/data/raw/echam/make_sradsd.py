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
path_srads = sys.argv[1]
path_sradsu = sys.argv[2]
path_sradsd = sys.argv[3]

# open files
file_srads = Dataset(path_srads, 'r')
file_sradsu = Dataset(path_sradsu, 'r')

# read data
srads = file_srads.variables['srads'][:] # (mon x lat x lon)
sradsu = file_sradsu.variables['sradsu'][:] # (mon x lat x lon)

# compute downwelling sw surface flux
sw_sfc_down = srads - sradsu

# save file as netCDF
file_sradsd = Dataset(path_sradsd, "w", format='NETCDF4_CLASSIC')

# copy attributes from srads file
file_sradsd.setncatts(file_srads.__dict__)

# copy dimensions from srads file
for name, dimension in file_srads.dimensions.items():
    file_sradsd.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from srads file
for name, variable in file_srads.variables.items():
    if any(name in s for s in ['srads']):
        continue
    
    x = file_sradsd.createVariable(name, variable.datatype, variable.dimensions)
    file_sradsd[name].setncatts(file_srads[name].__dict__)
    file_sradsd[name][:] = file_srads[name][:]
    
sradsd = file_sradsd.createVariable('sradsd', 'f4', ("time","lat","lon"))
sradsd.units = "W m**-2"
sradsd.long_name = "downwelling surface shortwave flux"
sradsd[:] = sw_sfc_down

file_sradsd.close()
