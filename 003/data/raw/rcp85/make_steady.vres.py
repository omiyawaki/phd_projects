import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps and mse files
path_ps = sys.argv[1]
path_va = sys.argv[2]
path_vas = sys.argv[3]
path_mse = sys.argv[4]
path_mses = sys.argv[5]
path_vmmmc = sys.argv[6]
path_vmse = sys.argv[7]

# open files
file_ps = Dataset(path_ps, 'r')
file_va = Dataset(path_va, 'r')
file_vas = Dataset(path_vas, 'r')
file_mse = Dataset(path_mse, 'r')
file_mses = Dataset(path_mses, 'r')

# read data
ps = file_ps.variables['ps'][:] # (mon x lat x lon)
va = file_va.variables['va'][:] # (mon x lev x lat x lon)
vas = file_vas.variables['vas'][:] # (mon x lev x lat x lon)
mse = file_mse.variables['mse'][:] # (mon x lev x lat x lon)
mses = file_mses.variables['mses'][:] # (mon x lev x lat x lon)
plev = file_mse.variables['plev'][:]

# for datasets that fill data below surface as missing data, fill with nans
va = va.filled(fill_value=np.nan)
vas = vas.filled(fill_value=np.nan)
mse = mse.filled(fill_value=np.nan)
mses = mses.filled(fill_value=np.nan)

# check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
lat2d = file_ps.variables['lat'][:] # (mon x lat x lon)
lat3d = file_mse.variables['lat'][:] # (mon x lev x lat x lon)
if not np.array_equal(lat2d,lat3d):
    filledps = ps.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filledps, axis=1)
    ps = f(filledlat3d)
    f = interpolate.interp1d(filledlat2d, vas, axis=1)
    vas = f(filledlat3d)
    filledps = None; f = None;

# replace subsurface data in general with nans (important for computing dhdt near the surface)
ps3d = np.tile(ps, [plev.size, 1, 1, 1])
ps3d = np.transpose( ps3d, [1, 0, 2, 3] )
pa3d = np.tile(plev, [ps.shape[0], ps.shape[1], ps.shape[2], 1])
pa3d = np.transpose( pa3d, [0, 3, 1, 2] )
idx_below = pa3d > ps3d
ps3d = None; pa3d = None;

va[idx_below]=np.nan
mse[idx_below]=np.nan

# take zonal mean
vas_z = np.mean(vas, axis=2)
mses_z = np.mean(mses, axis=2)
va_z = np.nanmean(va, axis=3)
mse_z = np.nanmean(mse, axis=3)

# compute mean meridional circulation transport
vm_mmc = va_z * mse_z

# compute stationary eddy transport
vms_se = np.squeeze(np.mean((vas - vas_z[...,np.newaxis]) * (mses - mses_z[...,np.newaxis]), axis=2))
vm_se = np.squeeze(np.mean((va - va_z[...,np.newaxis]) * (mse - mse_z[...,np.newaxis]), axis=3))

# compute vertical deviation of v and mse following Marshall et al. (2014)
ps_z = np.mean(ps, axis=2)
ps_za = np.mean(ps, axis=(0,2))

# # compute MSE flux divergences
rlat = np.radians(lat3d)
clat = np.cos(rlat)

# save file as netCDF
file_vmmmc = Dataset(path_vmmmc, "w", format='NETCDF4_CLASSIC')

# # copy attributes from mse file
# file_vmmmc.setncatts(file_mse.__dict__)

# copy dimensions from mse file
for name, dimension in file_mse.dimensions.items():
    file_vmmmc.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables from mse file
for name, variable in file_mse.variables.items():
    if any(name in s for s in ['mse']):
        continue
    
    x = file_vmmmc.createVariable(name, variable.datatype, variable.dimensions)
    file_vmmmc[name].setncatts(file_mse[name].__dict__)
    file_vmmmc[name][:] = file_mse[name][:]
    
vmmmc = file_vmmmc.createVariable('vmmmc', 'f4', ("time","plev","lat"))
vmmmc.units = "C m s**-1"
vmmmc.long_name = "northward transport of MSE/cp"
vmmmc[:,:,:] = vm_mmc

file_vmmmc.close()
