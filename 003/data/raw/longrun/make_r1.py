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
path_ra = sys.argv[1]
path_stf = sys.argv[2]
path_r1 = sys.argv[3]

# open files
print('Loading ra file...')
file_ra = Dataset(path_ra, 'r')
print('Done.\n')

print('Loading stf file...')
file_stf = Dataset(path_stf, 'r')
print('Done.\n')

# read data
ra = file_ra.variables['ra'][:] # (mon x lat x lon)
stf = file_stf.variables['stf'][:] # (mon x lat x lon)

lat2d = file_ra.variables['lat']
# # fill mask with nans                                                           
# lat3d = lat3d.filled(np.nan)
# lat2d = lat2d.filled(np.nan)                                                    
# take zonal means
ra_z = np.nanmean(ra, axis=2)
stf_z = np.nanmean(stf, axis=2)

# compute zonal mean r1
r1_z = 1+stf_z/ra_z

# save files as netCDF
file_r1 = Dataset(path_r1, "w", format='NETCDF4_CLASSIC')

# copy attributes from ra file
file_r1.setncatts(file_ra.__dict__)

# copy dimensions from ra file
for name, dimension in file_ra.dimensions.items():
    if any(name in s for s in ['lon']):
        continue
    file_r1.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

for name, variable in file_ra.variables.items():
    if any(name in s for s in ['ra', 'lon', 'lon_bnds']):
        continue

    x = file_r1.createVariable(name, variable.datatype, variable.dimensions)
    file_r1[name].setncatts(file_ra[name].__dict__)
    file_r1[name][:] = file_ra[name][:]

r1 = file_r1.createVariable('r1', 'f4', ("time","lat"))
r1.units = "unitless"
r1.long_name = "relative importance of advective heating to radiative cooling"
r1[:,:] = r1_z
file_r1.close()
