import sys
import numpy as np
import xarray as xr
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

################################################################################
# Constants and parameters
################################################################################
cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;
p0 = 1100.e2 # bottom integral bound [hPa]
pt = 0. # top integral bound [hPa]

# Uniform mse v offset test (acts as a sanity check because transport should be insensitive to constant V or MSE offsets if analysis is implemented correctly)
do_uniform_offset = False
mse_offset = 1e7
va_offset=1e2

# Choose method for treating subsurface data
# Method 1 discards all subsurface data (any existing subsurface mask plus p > ps)
# Method 2 linearly extrapolates for surface data
method_subsurface = 1

# Choose method for computing transport across latitude band
# Method 1 vertically and zonally integrates the transport at each longtitude and level (accounts for zonal variations in surface and its covariation with stationary eddies)
# Method 2 vertically integrates the zonal mean transport (using the zonal mean beta mask)
method_integral = 2

################################################################################
# Import data
################################################################################

# Paths to ps and mse files
path_va = sys.argv[1]
path_mse = sys.argv[2]
path_cbeta = sys.argv[3] 
path_mmc = sys.argv[4]
path_se = sys.argv[5]

# Open files
file_va = Dataset(path_va, 'r')
file_mse = Dataset(path_mse, 'r')
file_cbeta = Dataset(path_cbeta, 'r')

# Read data
va = file_va.variables['va'][:] # (mon x lev x lat x lon)
mse = file_mse.variables['mse'][:] # (mon x lev x lat x lon)
beta = file_cbeta.variables['beta'][:] # (1 x lev x lat x lon) climatological beta

################################################################################
# Uniform mse v offset test 
################################################################################
if do_uniform_offset:
    print('Adding uniform %g J/kg to MSE field and %g m/s to v field...' % (mse_offset, va_offset))
    mse = mse + mse_offset
    va = va + va_offset

################################################################################
# Grid data
################################################################################
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

################################################################################
# What to do if subsurface data is already masked
################################################################################
if method_subsurface == 1:
    #
    # METHOD 1: discard subsurface data
    # make sure that the beta mask is consistent with the existing subsurface mask
    mask = ~np.isnan(mse.filled(fill_value=np.nan))
    beta = beta * mask

elif method_subsurface == 2:
    #
    # METHOD 2: linearly extrapolate to determine surface data
    # If the dataset already has a subsurface mask, extrapolate (subsurface values will be treated by interpolating the surface value based on the climatological annual mean surface pressure field
    #
    va = va.filled(fill_value=np.nan)
    if np.any(np.isnan(va)):
        # reopen data using xarray
        va = xr.open_dataset(path_va).va # (mon x plev x lat x lon)
        # fill missing subsurface data with extrapolation
        va = va.reindex(plev=va['plev'][::-1])
        va = va.interpolate_na(dim="plev", method='linear', fill_value="extrapolate")
        va = va.reindex(plev=va['plev'][::-1])
        va = va.values

    mse = mse.filled(fill_value=np.nan)
    if np.any(np.isnan(mse)):
        mse = xr.open_dataset(path_mse).mse # (mon x plev x lat x lon)
        mse = mse.reindex(plev=mse['plev'][::-1])
        mse = mse.interpolate_na(dim="plev", method='linear', fill_value="extrapolate")
        mse = mse.reindex(plev=mse['plev'][::-1])
        mse = mse.values

################################################################################
# Take zonal mean
################################################################################
deriv_z = np.nansum(beta, axis=3, keepdims=True)
va_z = np.divide( np.nansum(beta*va, axis=3, keepdims=True), deriv_z, out=np.zeros_like(va), where=deriv_z!=0)
mse_z = np.divide( np.nansum(beta*mse, axis=3, keepdims=True), deriv_z, out=np.zeros_like(mse), where=deriv_z!=0)
beta_z = np.mean(beta, axis=3, keepdims=True)

################################################################################
# Compute stationary eddy (SE) MSE transport
################################################################################
vm_se = (va - va_z) * (mse - mse_z)

clat_ext = clat[None,:] # extend clat dimensions (mon x lat)
dplev_ext = dplev[None,:,None,None] # extend dplev dimensions (mon x lev x lat x lon)

if method_integral == 1:
    # METHOD 1:
    # zonally and vertically integrate the zonally asymmetric VM_SE transport
    # (takes into account the covariation in zonal asymmetric VM_SE and beta)

    vm_se_z_vint = a*clat_ext*dlon/g * np.nansum(beta*vm_se*dplev_ext, axis=(1,3))

elif method_integral == 2:
    # METHOD 2:
    # vertically integrate the zonal mean VM_SE transport

    vm_se_z = np.divide( np.nansum(beta*vm_se, axis=3, keepdims=True), deriv_z, out=np.zeros([vm_se.shape[0], vm_se.shape[1], vm_se.shape[2], 1]), where=deriv_z!=0)
    vm_se_z_vint = 2*np.pi*a*clat_ext/g * np.nansum(np.squeeze( beta_z*vm_se_z*dplev_ext) , axis=1)

################################################################################
# Compute mean meridional circulation (MMC) MSE transport
################################################################################

# compute deviation of vertical mean zonal mean v and mse
denom_zv = np.nansum(beta*dplev_ext, axis=(1,3), keepdims=True)
va_zv = np.divide( np.nansum(beta*va*dplev_ext, axis=(1,3), keepdims=True), denom_zv, out=np.zeros_like(va), where=denom_zv!=0)
mse_zv = np.divide( np.nansum(beta*mse*dplev_ext, axis=(1,3), keepdims=True), denom_zv, out=np.zeros_like(mse), where=denom_zv!=0)

# compute mean meridional circulation (MMC) transport following Marshall et al. (2014)
vm_mmc = (va_z - va_zv) * (mse_z - mse_zv)

if method_integral == 1:
    #
    # METHOD 1:
    # zonally and vertically integrate VM_MMC transport accounting for zonal asymmetries in beta
    #
    vm_mmc_z_vint = a*clat_ext*dlon/g * np.nansum(beta*vm_mmc*dplev_ext, axis=(1,3))

elif method_integral == 2:
    #
    # METHOD 2:
    # vertically integrate the VM_MMC transport using the zonal mean beta
    #
    vm_mmc_z = np.divide( np.nansum(beta*vm_mmc, axis=3, keepdims=True), deriv_z, out=np.zeros([vm_mmc.shape[0], vm_mmc.shape[1], vm_mmc.shape[2], 1]), where=deriv_z!=0)
    vm_mmc_z_vint = 2*np.pi*a*clat_ext/g * np.nansum(np.squeeze(beta_z*vm_mmc_z*dplev_ext), axis=1)

################################################################################
# Correct sign of the transport if necessary
################################################################################
# reverse the sign of the integration if data is stored in plev descending order
if plev[1]-plev[0]<0: # if pressure increases with index
    vm_se_z_vint = -vm_se_z_vint
    vm_mmc_z_vint = -vm_mmc_z_vint

################################################################################
# SAVE DATA
################################################################################
# save file as netCDF
file_mmc = Dataset(path_mmc, "w", format='NETCDF4_CLASSIC')
file_se = Dataset(path_se, "w", format='NETCDF4_CLASSIC')

# copy attributes from mse file
file_mmc.setncatts(file_mse.__dict__)
file_se.setncatts(file_mse.__dict__)

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
    
mmc = file_mmc.createVariable('vmmmc', 'f4', ("time","lat"))
mmc.units = "W"
mmc.long_name = "vertically integrated moist static energy flux transport due to mean meridional circulation"
mmc[:,:] = vm_mmc_z_vint

se = file_se.createVariable('vmse', 'f4', ("time","lat"))
se.units = "W"
se.long_name = "vertically integrated moist static energy flux transport due to stationary eddies"
se[:,:] = vm_se_z_vint

file_mmc.close()
file_se.close()
