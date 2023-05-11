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
path_traf0 = sys.argv[1]
path_sraf0 = sys.argv[2]
path_trafs = sys.argv[3]
path_srafs = sys.argv[4]
path_ra = sys.argv[5]

# open files
print('Loading traf0 file...')
file_traf0 = Dataset(path_traf0, 'r')
print('Done.\n')

print('Loading sraf0 file...')
file_sraf0 = Dataset(path_sraf0, 'r')
print('Done.\n')

print('Loading trafs file...')
file_trafs = Dataset(path_trafs, 'r')
print('Done.\n')

print('Loading srafs file...')
file_srafs = Dataset(path_srafs, 'r')
print('Done.\n')

# read data
traf0 = file_traf0.variables['traf0'][:] # (mon x lat x lon)
sraf0 = file_sraf0.variables['sraf0'][:] # (mon x lat x lon)
trafs = file_trafs.variables['trafs'][:] # (mon x lat x lon)
srafs = file_srafs.variables['srafs'][:] # (mon x lat x lon)

# compute radiative cooling
racool = traf0 - trafs + sraf0 - srafs

# save file as netCDF
file_ra = Dataset(path_ra, "w", format='NETCDF4_CLASSIC')

# copy attributes from traf0 file
file_ra.setncatts(file_traf0.__dict__)

# copy dimensions from traf0 file
for name, dimension in file_traf0.dimensions.items():
    file_ra.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from traf0 file
for name, variable in file_traf0.variables.items():
    if any(name in s for s in ['traf0']):
        continue
    
    x = file_ra.createVariable(name, variable.datatype, variable.dimensions)
    file_ra[name].setncatts(file_traf0[name].__dict__)
    file_ra[name][:] = file_traf0[name][:]
    
ra = file_ra.createVariable('racs', 'f4', ("time","lat","lon"))
ra.units = "W m**-2"
ra.long_name = "radiative cooling (clear sky)"
ra[:] = racool

file_ra.close()
