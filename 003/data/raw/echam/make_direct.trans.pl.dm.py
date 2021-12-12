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
path_atm = sys.argv[1]
path_mse = sys.argv[2]
path_te = sys.argv[3]
path_ps = sys.argv[4]

# open files
print('Loading ATM file...')
file_atm = Dataset(path_atm, 'r')
print('Done.\n')

print('Loading mse file...')
file_mse = Dataset(path_mse, 'r')
print('Done.\n')

print('Loading ps file...')
file_ps = Dataset(path_ps, 'r')
print('Done.\n')

# read data
ps = file_ps.variables['aps'][:] # (mon x lat x lon)
va = file_atm.variables['v'][:] # (mon x lev x lat x lon)
mse = file_mse.variables['mse'][:] # (mon x lev x lat x lon)
lat = file_atm.variables['lat'][:] # (lat)
plev = file_atm.variables['lev'][:] # (lev)

# # for datasets that fill data below surface as missing data, fill with nans
plev = plev.filled(fill_value=np.nan)
ps = ps.filled(fill_value=np.nan)
va = va.filled(fill_value=np.nan)
mse = mse.filled(fill_value=np.nan)

# Zonal and monthly means of surface pressure
ps_z = np.nanmean(ps, axis=2)
ps_t = np.nanmean(ps, axis=(0))
ps_za = np.nanmean(ps, axis=(0,2))

ps_t_tile = np.tile(ps_t, [ps.shape[0], 1, 1])

# mask out data below surface
pa = np.transpose(np.tile(plev, [ps.shape[0], ps.shape[1], ps.shape[2], 1]), [0,3,1,2])

# subsurf = pa > np.transpose(np.tile(ps_a_tile, [len(plev),1,1,1]),[1,0,2,3])
# subsurf = pa > np.transpose(np.tile(ps_t_tile, [len(plev),1,1,1]),[1,0,2,3])
subsurf = pa > np.transpose(np.tile(ps, [len(plev),1,1,1]),[1,0,2,3])

va[subsurf] = np.nan
mse[subsurf] = np.nan
# pa[subsurf] = np.nan

# take monthly mean
va_t = np.nanmean(va, axis=0)
mse_t = np.nanmean(mse, axis=0)
pa_t = np.nanmean(pa, axis=0)

# take zonal mean
pa_zt = np.nanmean(pa_t, axis=2)
va_zt = np.nanmean(va_t, axis=2)
mse_zt = np.nanmean(mse_t, axis=2)

# deviation from monthly mean zonal mean
va_dt = va - va_t
mse_dt = mse - mse_t
va_dzt = va_dt - np.nanmean(va_dt, axis=3)[...,np.newaxis]
mse_dzt = mse_dt - np.nanmean(mse_dt, axis=3)[...,np.newaxis]

# compute transient eddy transport
vm_te = np.nanmean(va_dzt * mse_dzt, axis=(0,3))

rlat = np.radians(lat)
clat = np.cos(rlat)

# replace nans with zeros in preparation for taking vertical integrals
vm_te[np.isnan(vm_te)] = 0

mask_zt = np.ones_like(va_zt)
mask_zt[np.isnan(va_zt)] = 0

va_zt[np.isnan(va_zt)] = 0
mse_zt[np.isnan(mse_zt)] = 0

# compute vertical integral of heat transports (simplified)
if plev[1]-plev[0]>0: # if pressure increases with index
    vm_te_vint = 2*np.pi*a*clat/g*np.trapz(vm_te, plev, axis=0)
else:
    vm_te_vint = -2*np.pi*a*clat/g*np.trapz(vm_te, plev, axis=0)

# save file as netCDF
file_te = Dataset(path_te, "w", format='NETCDF4_CLASSIC')

# copy attributes from mse file
file_te.setncatts(file_mse.__dict__)

# copy dimensions from mse file
for name, dimension in file_mse.dimensions.items():
    if any(name in s for s in ['lev', 'hyai', 'hybi', 'hyam', 'hybm']):
        continue
    file_te.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_mse.variables.items():
    if any(name in s for s in ['mse' 'lev', 'hyai', 'hybi', 'hyam', 'hybm']):
        continue
    
    x = file_te.createVariable(name, variable.datatype, variable.dimensions)
    file_te[name].setncatts(file_mse[name].__dict__)
    file_te[name][:] = file_mse[name][:]

te = file_te.createVariable('vmte', 'f4', ("time","lat"))
te.units = "W"
te.long_name = "vertically integrated moist static energy flux transport due to transient eddies"
te[:,:] = vm_te_vint

file_te.close()
