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
path_rlutcs = sys.argv[1]
path_rsdt = sys.argv[2]
path_rsutcs = sys.argv[3]
path_ftoa = sys.argv[4]

# open files
print('Loading rlutcs file...')
file_rlutcs = Dataset(path_rlutcs, 'r')
print('Done.\n')

print('Loading rsdt file...')
file_rsdt = Dataset(path_rsdt, 'r')
print('Done.\n')

print('Loading rsutcs file...')
file_rsutcs = Dataset(path_rsutcs, 'r')
print('Done.\n')

# read data
rlutcs = file_rlutcs.variables['rlutcs'][:] # (mon x lat x lon)
rsdt = file_rsdt.variables['rsdt'][:] # (mon x lat x lon)
rsutcs = file_rsutcs.variables['rsutcs'][:] # (mon x lat x lon)

# compute net TOA flux
toa = rsdt - rsutcs - rlutcs # positive heats planet

# save file as netCDF
file_ftoa = Dataset(path_ftoa, "w", format='NETCDF4_CLASSIC')

# copy attributes from rlutcs file
file_ftoa.setncatts(file_rlutcs.__dict__)

# copy dimensions from rlutcs file
for name, dimension in file_rlutcs.dimensions.items():
    file_ftoa.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from rlutcs file
for name, variable in file_rlutcs.variables.items():
    if any(name in s for s in ['rlutcs']):
        continue
    
    x = file_ftoa.createVariable(name, variable.datatype, variable.dimensions)
    file_ftoa[name].setncatts(file_rlutcs[name].__dict__)
    file_ftoa[name][:] = file_rlutcs[name][:]
    
ftoa = file_ftoa.createVariable('ftoacs', 'f4', ("time","lat","lon"))
ftoa.units = "W m**-2"
ftoa.long_name = "net clear sky TOA flux"
ftoa[:] = toa

file_ftoa.close()
