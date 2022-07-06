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
path_tend = sys.argv[3]
path_stgadv = sys.argv[4]
path_adv = sys.argv[5]
path_r1 = sys.argv[6]

# open files
print('Loading ra file...')
file_ra = Dataset(path_ra, 'r')
print('Done.\n')

print('Loading stf file...')
file_stf = Dataset(path_stf, 'r')
print('Done.\n')

print('Loading tend file...')
file_tend = Dataset(path_tend, 'r')
print('Done.\n')

# read data
ra = file_ra.variables['ra'][:] # (mon x lat x lon)
stf = file_stf.variables['stf'][:] # (mon x lat x lon)
tend = file_tend.variables['tend'][:] # (mon x lat x lon)

# interpolate tend (on 3d grid) to 2d grid if not the same                                
print(tend.shape)                                                       
lat3d = file_tend.variables['lat']
lat2d = file_ra.variables['lat']
# # fill mask with nans                                                           
# lat3d = lat3d.filled(np.nan)
# lat2d = lat2d.filled(np.nan)                                                    
if not np.array_equal(lat2d,lat3d):                                             
    tend = tend.filled(np.nan)                                      
    f_tend = interpolate.interp1d(lat3d, tend, axis=1, fill_value='extrapolate')
    tend = f_tend(lat2d)                                                
    print(tend.shape)

# compute storage plus advection
tend_p_divfm = ra + stf

# compute advection
divfm = tend_p_divfm - tend

# take zonal means
ra_z = np.nanmean(ra, axis=2)
stf_z = np.nanmean(stf, axis=2)
tend_p_divfm_z = np.nanmean(tend_p_divfm, axis=2)

# compute zonal mean r1
r1_z = tend_p_divfm_z / ra_z

# save files as netCDF
file_stgadv = Dataset(path_stgadv, "w", format='NETCDF4_CLASSIC')
file_adv = Dataset(path_adv, "w", format='NETCDF4_CLASSIC')
file_r1 = Dataset(path_r1, "w", format='NETCDF4_CLASSIC')

# copy attributes from ra file
file_stgadv.setncatts(file_ra.__dict__)
file_adv.setncatts(file_ra.__dict__)
file_r1.setncatts(file_ra.__dict__)

# copy dimensions from ra file
for name, dimension in file_ra.dimensions.items():
    file_stgadv.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_adv.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

for name, dimension in file_ra.dimensions.items():
    if any(name in s for s in ['lon']):
        continue
    file_r1.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

# copy all variables except time from ra file
for name, variable in file_ra.variables.items():
    if any(name in s for s in ['ra']):
        continue

    x = file_stgadv.createVariable(name, variable.datatype, variable.dimensions)
    file_stgadv[name].setncatts(file_ra[name].__dict__)
    file_stgadv[name][:] = file_ra[name][:]

    x = file_adv.createVariable(name, variable.datatype, variable.dimensions)
    file_adv[name].setncatts(file_ra[name].__dict__)
    file_adv[name][:] = file_ra[name][:]

stgadv = file_stgadv.createVariable('stgadv', 'f4', ("time","lat","lon"))
stgadv.units = "W m**-2"
stgadv.long_name = "MSE storage plus MSE flux divergence"
stgadv[:] = tend_p_divfm
file_stgadv.close()

adv = file_adv.createVariable('adv', 'f4', ("time","lat","lon"))
adv.units = "W m**-2"
adv.long_name = "MSE flux divergence"
adv[:] = divfm
file_adv.close()

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
