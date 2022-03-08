import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;
p0 = 1100.e2 # bottom integral bound [hPa]
pt = 0. # top integral bound [hPa]

# paths to ps and va files
path_ps = sys.argv[1]
path_va = sys.argv[2]
path_beta = sys.argv[3]
path_psi = sys.argv[4]

# open files
file_ps = Dataset(path_ps, 'r')
file_va = Dataset(path_va, 'r')
file_beta = Dataset(path_beta, 'r')

# read data
ps = np.squeeze(file_ps.variables['ps'][:]) # (lat x lon)
va = file_va.variables['va'][:] # (mon x lev x lat x lon)
beta = np.squeeze(file_beta.variables['beta'][:]) # (lev x lat x lon)

plev = file_va.variables['plev'][:]
plev_half = 1/2 * (plev[1:] + plev[:-1])
plev_half = np.sort(np.append(plev_half, [pt, p0]))
plev_full = np.sort(np.concatenate((plev, plev_half)))
if plev[1]-plev[0]<0:
    plev_half = plev_half[::-1]
    plev_full = plev_full[::-1]
dplev = plev_half[1:] - plev_half[:-1]

# for datasets that fill data below surface as missing data, fill with nans
va = va.filled(fill_value=np.nan)

# check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
lat2d = file_ps.variables['lat'][:] 
lat3d = file_va.variables['lat'][:] 
if not np.array_equal(lat2d,lat3d):
    filledps = ps.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filledps, axis=0)
    ps = f(filledlat3d)
    filledps = None; f = None;

# replace subsurface data in general with nans (important for computing dhdt near the surface)
ps_tile = np.tile(ps, [va.shape[0], va.shape[1], 1, 1])
pa = np.transpose(np.tile(plev, [va.shape[0], va.shape[2], va.shape[3], 1]), [0,3,1,2])
subsurf = pa > ps_tile
va[subsurf] = np.nan

# take zonal mean
beta = beta[None,...] # create time dim in beta
deriv_z = np.nansum(beta, axis=3)
va_z = np.divide( np.nansum(beta*va, axis=3), deriv_z, out=np.zeros([va.shape[0],va.shape[1],va.shape[2]]), where=deriv_z!=0)
beta_z = np.mean(beta, axis=3)

# compute streamfunction
rlat = np.radians(lat3d)
clat = np.cos(rlat)
if plev[1]-plev[0]>0: # if pressure increases with index
    psi = 2*np.pi*a*clat[None,None,:]/g * np.cumsum( beta_z * va_z * dplev[None,...,None], axis=1)
else:
    psi = -2*np.pi*a*clat[None,None,:]/g * np.cumsum( beta_z[:,::-1,:] * va_z[:,::-1,:] * dplev[None,::-1,None], axis=1)
    psi = psi[:,::-1,:]

# save file as netCDF
file_psi = Dataset(path_psi, "w", format='NETCDF4_CLASSIC')

# copy attributes from va file
file_psi.setncatts(file_va.__dict__)

# copy dimensions from mse file
for name, dimension in file_va.dimensions.items():
    if any(name in s for s in ['lon']):
        continue
    file_psi.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_va.variables.items():
    if any(name in s for s in ['va', 'lon', 'lon_bnds']):
        continue
    
    x = file_psi.createVariable(name, variable.datatype, variable.dimensions)
    file_psi[name].setncatts(file_va[name].__dict__)
    file_psi[name][:] = file_va[name][:]
    
sfunc = file_psi.createVariable('psi', 'f4', ("time","plev","lat"))
sfunc.units = "kg/s"
sfunc.long_name = "mean meridional mass streamfunction"
sfunc[:,:,:] = psi

file_psi.close()
