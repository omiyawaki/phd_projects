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
path_ahfl = sys.argv[1]
path_orog = sys.argv[2]

# open files
file_ahfl = Dataset(path_ahfl, 'r')

# read data
ahfl = file_ahfl.variables['ahfl'][:] # (mon x lat x lon)

# compute surface turbulent fluxes
zsurf = np.zeros([ahfl.shape[1], ahfl.shape[2]])

# save file as netCDF
file_orog = Dataset(path_orog, "w", format='NETCDF4_CLASSIC')

# copy attributes from ahfl file
file_orog.setncatts(file_ahfl.__dict__)

# copy dimensions from ahfl file
for name, dimension in file_ahfl.dimensions.items():
    file_orog.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from ahfl file
for name, variable in file_ahfl.variables.items():
    if any(name in s for s in ['ahfl']):
        continue
    
    x = file_orog.createVariable(name, variable.datatype, variable.dimensions)
    file_orog[name].setncatts(file_ahfl[name].__dict__)
    file_orog[name][:] = file_ahfl[name][:]
    
orog = file_orog.createVariable('orog', 'f4', ("lat","lon"))
orog.units = "m"
orog.long_name = "surface orography"
orog[:] = zsurf

file_orog.close()
