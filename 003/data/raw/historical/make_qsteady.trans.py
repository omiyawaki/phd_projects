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

# Paths to ps and hus files
path_ps = sys.argv[1]
path_va = sys.argv[2]
path_hus = sys.argv[3]
path_beta = sys.argv[4]
path_mmc = sys.argv[5]
path_se = sys.argv[6]

# Open files
file_ps = Dataset(path_ps, 'r')
file_va = Dataset(path_va, 'r')
file_hus = Dataset(path_hus, 'r')
file_beta = Dataset(path_beta, 'r')

# Read data
ps = np.squeeze(file_ps.variables['ps'][:]) # (lat x lon)
va = file_va.variables['va'][:] # (mon x lev x lat x lon)
hus = L * file_hus.variables['hus'][:] # (mon x lev x lat x lon)
beta = file_beta.variables['beta'][:] # (1 x lev x lat x lon)

# Uniform hus v offset test
# hus = hus + 1e7
# va = va+ 1e2

# For datasets that fill data below surface as missing data, set as nans
va = va.filled(fill_value=np.nan)
hus = hus.filled(fill_value=np.nan)

# Grid data
lat2d = file_ps.variables['lat'][:] # (mon x lat x lon)
lat3d = file_va.variables['lat'][:] # (mon x lev x lat x lon)
rlat = np.radians(lat3d)
clat = np.cos(rlat)

lon = file_va.variables['lon'][:]
dlon = np.pi / 180 * (lon[1] - lon[0])

plev = file_va.variables['plev'][:]
plev_half = 1/2 * (plev[1:] + plev[:-1])
plev_half = np.sort(np.append(plev_half, [pt, p0]))
plev_full = np.sort(np.concatenate((plev, plev_half)))
if plev[1]-plev[0]<0:
    plev_half = plev_half[::-1]
    plev_full = plev_full[::-1]
dplev = plev_half[1:] - plev_half[:-1]

# Check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
if not np.array_equal(lat2d,lat3d):
    print('\nInterpolating ps to 3d lat grid...\n')
    filledps = ps.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filledps, axis=0)
    ps = f(filledlat3d)
    filledps = None; f = None;
    print('\nDone.\n')

# Zonal surface pressure
ps_z = np.nanmean(ps, axis=1)
ps_tile = np.tile(ps, [va.shape[0], va.shape[1], 1, 1])

# Mask out data below surface using long-term annual mean surface pressure
pa = np.transpose(np.tile(plev, [va.shape[0], va.shape[2], va.shape[3], 1]), [0,3,1,2])
subsurf = pa > ps_tile
va[subsurf] = np.nan
hus[subsurf] = np.nan
beta[subsurf[[0],...]] = np.nan
pa = None

# take zonal mean
va_z = np.nanmean(va, axis=3, keepdims=True)
hus_z = np.nanmean(hus, axis=3, keepdims=True)
beta_z = np.nanmean(beta, axis=3, keepdims=True)

# compute stationary eddy (SE) hus transport
vm_se = np.nanmean((va - va_z) * (hus - hus_z), axis=3)

# vertically and zonally integrate SE hus transport
clat_ext = clat[None,:] # extend clat dimensions (mon x lat)
dplev_ext = dplev[None,:,None,None] # extend dplev dimensions (mon x lev x lat x lon)
vm_se_z_vint = a*clat_ext*dlon/g * np.nansum(beta*vm_se[...,None]*dplev_ext, axis=(1,3))

# compute deviation of vertical mean zonal mean v and hus (used to compute MMC transport)
denom_zv = np.nansum(beta*dplev_ext, axis=1, keepdims=True)
va_zv = np.nansum(beta*va_z*dplev_ext, axis=1, keepdims=True) / denom_zv
hus_zv = np.nansum(beta*hus_z*dplev_ext, axis=1, keepdims=True) / denom_zv

# compute mean meridional circulation (MMC) transport following Marshall et al. (2014)
vm_mmc = (va_z - va_zv) * (hus_z - hus_zv)

# vertically and zonally integrate MMC hus transport
vm_mmc_z_vint = a*clat_ext*dlon/g * np.nansum(beta*vm_mmc*dplev_ext, axis=(1,3))

# reverse the sign of the integration if data is stored in plev descending order
if plev[1]-plev[0]<0: # if pressure increases with index
    vm_se_z_vint = -vm_se_z_vint
    vm_mmc_z_vint = -vm_mmc_z_vint

# save file as netCDF
file_mmc = Dataset(path_mmc, "w", format='NETCDF4_CLASSIC')
file_se = Dataset(path_se, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_mmc.setncatts(file_hus.__dict__)
file_se.setncatts(file_hus.__dict__)

# copy dimensions from va file
for name, dimension in file_va.dimensions.items():
    if any(name in s for s in ['plev']):
        continue

    file_mmc.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_se.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from va file
for name, variable in file_va.variables.items():
    if any(name in s for s in ['va', 'plev', 'plev_bnds']):
        continue
    
    x = file_mmc.createVariable(name, variable.datatype, variable.dimensions)
    file_mmc[name].setncatts(file_va[name].__dict__)
    file_mmc[name][:] = file_va[name][:]

    x = file_se.createVariable(name, variable.datatype, variable.dimensions)
    file_se[name].setncatts(file_va[name].__dict__)
    file_se[name][:] = file_va[name][:]
    
mmc = file_mmc.createVariable('vqmmc', 'f4', ("time","lat"))
mmc.units = "W"
mmc.long_name = "vertically integrated latent energy flux energy flux transport due to mean meridional circulation"
mmc[:,:] = vm_mmc_z_vint

se = file_se.createVariable('vqse', 'f4', ("time","lat"))
se.units = "W"
se.long_name = "vertically integrated latent energy flux transport due to stationary eddies"
se[:,:] = vm_se_z_vint

file_mmc.close()
file_se.close()
