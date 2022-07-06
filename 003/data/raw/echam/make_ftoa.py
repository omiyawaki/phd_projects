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
path_trad0 = sys.argv[1]
path_srad0 = sys.argv[2]
path_ftoa = sys.argv[3]

# open files
file_srad0 = Dataset(path_srad0, 'r')
file_trad0 = Dataset(path_trad0, 'r')

# read data
srad0 = file_srad0.variables['srad0'][:] # (mon x lat x lon)
trad0 = file_trad0.variables['trad0'][:] # (mon x lat x lon)

# compute net TOA radiative flux
toa = trad0 + srad0

# save file as netCDF
file_ftoa = Dataset(path_ftoa, "w", format='NETCDF4_CLASSIC')

# copy attributes from trad0 file
file_ftoa.setncatts(file_trad0.__dict__)

# copy dimensions from trad0 file
for name, dimension in file_trad0.dimensions.items():
    file_ftoa.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from trad0 file
for name, variable in file_trad0.variables.items():
    if any(name in s for s in ['trad0']):
        continue
    
    x = file_ftoa.createVariable(name, variable.datatype, variable.dimensions)
    file_ftoa[name].setncatts(file_trad0[name].__dict__)
    file_ftoa[name][:] = file_trad0[name][:]
    
ftoa = file_ftoa.createVariable('ftoa', 'f4', ("time","lat","lon"))
ftoa.units = "W m**-2"
ftoa.long_name = "net TOA radiative flux"
ftoa[:] = toa

file_ftoa.close()
