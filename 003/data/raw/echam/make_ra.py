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
path_trads = sys.argv[3]
path_srads = sys.argv[4]
path_ra = sys.argv[5]

# open files
print('Loading trad0 file...')
file_trad0 = Dataset(path_trad0, 'r')
print('Done.\n')

print('Loading srad0 file...')
file_srad0 = Dataset(path_srad0, 'r')
print('Done.\n')

print('Loading trads file...')
file_trads = Dataset(path_trads, 'r')
print('Done.\n')

print('Loading srads file...')
file_srads = Dataset(path_srads, 'r')
print('Done.\n')

# read data
trad0 = file_trad0.variables['trad0'][:] # (mon x lat x lon)
srad0 = file_srad0.variables['srad0'][:] # (mon x lat x lon)
trads = file_trads.variables['trads'][:] # (mon x lat x lon)
srads = file_srads.variables['srads'][:] # (mon x lat x lon)

# compute radiative cooling
racool = trad0 - trads + srad0 - srads

# save file as netCDF
file_ra = Dataset(path_ra, "w", format='NETCDF4_CLASSIC')

# copy attributes from trad0 file
file_ra.setncatts(file_trad0.__dict__)

# copy dimensions from trad0 file
for name, dimension in file_trad0.dimensions.items():
    file_ra.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from trad0 file
for name, variable in file_trad0.variables.items():
    if any(name in s for s in ['trad0']):
        continue
    
    x = file_ra.createVariable(name, variable.datatype, variable.dimensions)
    file_ra[name].setncatts(file_trad0[name].__dict__)
    file_ra[name][:] = file_trad0[name][:]
    
ra = file_ra.createVariable('ra', 'f4', ("time","lat","lon"))
ra.units = "W m**-2"
ra.long_name = "radiative cooling"
ra[:] = racool

file_ra.close()
