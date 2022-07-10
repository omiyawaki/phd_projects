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
path_rlut = sys.argv[1]
path_rsdt = sys.argv[2]
path_rsut = sys.argv[3]
path_ftoa = sys.argv[4]

# open files
print('Loading rlut file...')
file_rlut = Dataset(path_rlut, 'r')
print('Done.\n')

print('Loading rsdt file...')
file_rsdt = Dataset(path_rsdt, 'r')
print('Done.\n')

print('Loading rsut file...')
file_rsut = Dataset(path_rsut, 'r')
print('Done.\n')

# read data
rlut = file_rlut.variables['rlut'][:] # (mon x lat x lon)
rsdt = file_rsdt.variables['rsdt'][:] # (mon x lat x lon)
rsut = file_rsut.variables['rsut'][:] # (mon x lat x lon)

# compute net TOA flux
toa = rsdt - rsut - rlut # positive heats planet

# save file as netCDF
file_ftoa = Dataset(path_ftoa, "w", format='NETCDF4_CLASSIC')

# copy attributes from rlut file
file_ftoa.setncatts(file_rlut.__dict__)

# copy dimensions from rlut file
for name, dimension in file_rlut.dimensions.items():
    file_ftoa.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from rlut file
for name, variable in file_rlut.variables.items():
    if any(name in s for s in ['rlut']):
        continue
    
    x = file_ftoa.createVariable(name, variable.datatype, variable.dimensions)
    file_ftoa[name].setncatts(file_rlut[name].__dict__)
    file_ftoa[name][:] = file_rlut[name][:]
    
ftoa = file_ftoa.createVariable('ftoa', 'f4', ("time","lat","lon"))
ftoa.units = "W m**-2"
ftoa.long_name = "net TOA flux"
ftoa[:] = toa

file_ftoa.close()
