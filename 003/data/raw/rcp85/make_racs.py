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
path_rsutcs = sys.argv[2]
path_rsdt = sys.argv[3]
path_rlus = sys.argv[4]
path_rldscs = sys.argv[5]
path_rsdscs = sys.argv[6]
path_rsuscs = sys.argv[7]
path_racs = sys.argv[8]

# open files
print('Loading rlutcs file...')
file_rlutcs = Dataset(path_rlutcs, 'r')
print('Done.\n')

print('Loading rsutcs file...')
file_rsutcs = Dataset(path_rsutcs, 'r')
print('Done.\n')

print('Loading rsdt file...')
file_rsdt = Dataset(path_rsdt, 'r')
print('Done.\n')

print('Loading rlus file...')
file_rlus = Dataset(path_rlus, 'r')
print('Done.\n')

print('Loading rldscs file...')
file_rldscs = Dataset(path_rldscs, 'r')
print('Done.\n')

print('Loading rsdscs file...')
file_rsdscs = Dataset(path_rsdscs, 'r')
print('Done.\n')

print('Loading rsuscs file...')
file_rsuscs = Dataset(path_rsuscs, 'r')
print('Done.\n')

# read data
rlutcs = file_rlutcs.variables['rlutcs'][:] # (mon x lat x lon)
rsutcs = file_rsutcs.variables['rsutcs'][:] # (mon x lat x lon)
rsdt = file_rsdt.variables['rsdt'][:] # (mon x lat x lon)
rlus = file_rlus.variables['rlus'][:] # (mon x lat x lon)
rldscs = file_rldscs.variables['rldscs'][:] # (mon x lat x lon)
rsdscs = file_rsdscs.variables['rsdscs'][:] # (mon x lat x lon)
rsuscs = file_rsuscs.variables['rsuscs'][:] # (mon x lat x lon)

# compute radiative cooling
racool = rsdt - rsutcs - rlutcs + rsuscs - rsdscs + rlus - rldscs

# save file as netCDF
file_racs = Dataset(path_racs, "w", format='NETCDF4_CLASSIC')

# copy attributes from rlutcs file
file_racs.setncatts(file_rlutcs.__dict__)

# copy dimensions from rlutcs file
for name, dimension in file_rlutcs.dimensions.items():
    file_racs.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from rlutcs file
for name, variable in file_rlutcs.variables.items():
    if any(name in s for s in ['rlutcs']):
        continue
    
    x = file_racs.createVariable(name, variable.datatype, variable.dimensions)
    file_racs[name].setncatts(file_rlutcs[name].__dict__)
    file_racs[name][:] = file_rlutcs[name][:]
    
racs = file_racs.createVariable('racs', 'f4', ("time","lat","lon"))
racs.units = "W m**-2"
racs.long_name = "clear sky radiative cooling"
racs[:] = racool

file_racs.close()
