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
path_ahfs = sys.argv[2]
path_stf = sys.argv[3]

# open files
print('Loading ahfl file...')
file_ahfl = Dataset(path_ahfl, 'r')
print('Done.\n')

print('Loading ahfs file...')
file_ahfs = Dataset(path_ahfs, 'r')
print('Done.\n')

# read data
ahfl = file_ahfl.variables['ahfl'][:] # (mon x lat x lon)
ahfs = file_ahfs.variables['ahfs'][:] # (mon x lat x lon)

# compute surface turbulent fluxes
lhsh = - (ahfl + ahfs)

# save file as netCDF
file_stf = Dataset(path_stf, "w", format='NETCDF4_CLASSIC')

# copy attributes from ahfl file
file_stf.setncatts(file_ahfl.__dict__)

# copy dimensions from ahfl file
for name, dimension in file_ahfl.dimensions.items():
    file_stf.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from ahfl file
for name, variable in file_ahfl.variables.items():
    if any(name in s for s in ['ahfl']):
        continue
    
    x = file_stf.createVariable(name, variable.datatype, variable.dimensions)
    file_stf[name].setncatts(file_ahfl[name].__dict__)
    file_stf[name][:] = file_ahfl[name][:]
    
stf = file_stf.createVariable('stf', 'f4', ("time","lat","lon"))
stf.units = "W m**-2"
stf.long_name = "surface latent plus sensible heat flux"
stf[:] = lhsh

file_stf.close()
