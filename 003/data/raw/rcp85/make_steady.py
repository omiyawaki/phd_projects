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
path_mmc = sys.argv[6]
path_se = sys.argv[7]

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

# # for datasets that fill data below surface as missing data, fill with nans
va = va.filled(fill_value=np.nan)
mse = mse.filled(fill_value=np.nan)

# check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
lat2d = file_ps.variables['lat'][:] # (mon x lat x lon)
lat3d = file_mse.variables['lat'][:] # (mon x lev x lat x lon)
if not np.array_equal(lat2d,lat3d):
    filledps = ps.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filledps, axis=1)
    ps = f(filledlat3d)
    filledps = None; f = None;

# Zonal and annual means of surface pressure
ps_z = np.mean(ps, axis=2)
ps_a = np.mean(ps, axis=(0))
ps_za = np.mean(ps, axis=(0,2))

ps_a_tile = np.tile(ps_a, [ps.shape[0], 1, 1])

# remove subsurface data
pa = np.transpose(np.tile(plev, [ps.shape[0], ps.shape[1], ps.shape[2], 1]), [0,3,1,2])
subsurf = pa > np.transpose(np.tile(ps, [len(plev),1,1,1]),[1,0,2,3])

va[subsurf] = np.nan
mse[subsurf] = np.nan

# take zonal mean
va_z = np.nanmean(va, axis=3)
vas_z = np.nanmean(vas, axis=2)
mse_z = np.nanmean(mse, axis=3)
mses_z = np.nanmean(mses, axis=2)

# compute stationary eddy transport
vms_mmc = vas_z * mses_z
vm_mmc = va_z * mse_z
vms_se = np.squeeze(np.nanmean((vas - vas_z[...,np.newaxis]) * (mses - mses_z[...,np.newaxis]), axis=2))
vm_se = np.squeeze(np.nanmean((va - va_z[...,np.newaxis]) * (mse - mse_z[...,np.newaxis]), axis=3))

# take meridional derivative
rlat = np.radians(lat3d)
clat = np.cos(rlat)

div_vms_mmc = np.empty_like(vms_mmc)
div_vm_mmc = np.empty_like(vm_mmc)
div_vms_se = np.empty_like(vms_se)
div_vm_se = np.empty_like(vm_se)

div_vms_mmc[...,1:-1] = 1/(a * 1/2*(clat[2:]+clat[:-2]) ) * (vms_mmc[...,2:] - vms_mmc[...,:-2]) / (rlat[2:] - rlat[:-2])
div_vm_mmc[...,1:-1] = 1/(a * 1/2*(clat[2:]+clat[:-2]) ) * (vm_mmc[...,2:] - vm_mmc[...,:-2]) / (rlat[2:] - rlat[:-2])
div_vms_se[...,1:-1] = 1/(a * 1/2*(clat[2:]+clat[:-2]) ) * (vms_se[...,2:] - vms_se[...,:-2]) / (rlat[2:] - rlat[:-2])
div_vm_se[...,1:-1] = 1/(a * 1/2*(clat[2:]+clat[:-2]) ) * (vm_se[...,2:] - vm_se[...,:-2]) / (rlat[2:] - rlat[:-2])
div_vms_mmc[...,0] = np.nan; div_vms_mmc[...,-1] = np.nan
div_vm_mmc[...,0] = np.nan; div_vm_mmc[...,-1] = np.nan
div_vms_se[...,0] = np.nan; div_vms_se[...,-1] = np.nan
div_vm_se[...,0] = np.nan; div_vm_se[...,-1] = np.nan

# compute vertical integral
div_vm_mmc_vint = np.empty_like(ps_z)
div_vm_se_vint = np.empty_like(ps_z)
for itime in tqdm(range(ps_z.shape[0])):
    for ilat in range(ps_z.shape[1]):
        div_vms_mmc_local = div_vms_mmc[itime, ilat]
        div_vms_se_local = div_vms_se[itime, ilat]

        # remove subsurface data
        abovesurf = (plev < ps_za[ilat])
        plev_local = plev[abovesurf]
        div_vm_mmc_local = div_vm_mmc[itime, abovesurf, ilat]
        div_vm_se_local = div_vm_se[itime, abovesurf, ilat]

        if plev[1]-plev[0]>0: # if pressure increases with index
            plev_local = np.append(plev_local, ps_za[ilat]) 
            div_vm_mmc_local = np.append(div_vm_mmc_local, div_vms_mmc_local)
            div_vm_se_local = np.append(div_vm_se_local, div_vms_se_local)

            div_mmc_vint[itime, ilat] = 1/g*np.trapz(div_vm_mmc_local, plev_local)
            div_se_vint[itime, ilat] = 1/g*np.trapz(div_vm_se_local, plev_local)
        else:
            plev_local = np.insert(plev_local, 0, ps_za[ilat]) 
            div_vm_mmc_local = np.insert(div_vm_mmc_local, 0, div_vms_mmc_local)
            div_vm_se_local = np.insert(div_vm_se_local, 0, div_vms_se_local)

            div_vm_mmc_vint[itime, ilat] = -1/g*np.trapz(div_vm_mmc_local, plev_local)
            div_vm_se_vint[itime, ilat] = -1/g*np.trapz(div_vm_se_local, plev_local)

# save file as netCDF
file_mmc = Dataset(path_mmc, "w", format='NETCDF4_CLASSIC')
file_se = Dataset(path_se, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_mmc.setncatts(file_ps.__dict__)
file_se.setncatts(file_ps.__dict__)

# copy dimensions from mse file
for name, dimension in file_mse.dimensions.items():
    if any(name in s for s in ['plev']):
        continue
    file_mmc.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_se.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_mse.variables.items():
    if any(name in s for s in ['mse' 'plev']):
        continue
    
    x = file_mmc.createVariable(name, variable.datatype, variable.dimensions)
    file_mmc[name].setncatts(file_mse[name].__dict__)
    file_mmc[name][:] = file_mse[name][:]

    x = file_se.createVariable(name, variable.datatype, variable.dimensions)
    file_se[name].setncatts(file_mse[name].__dict__)
    file_se[name][:] = file_mse[name][:]
    
mmc = file_mmc.createVariable('divmmc', 'f4', ("time","lat"))
mmc.units = "W m**-2"
mmc.long_name = "vertically integrated moist static energy flux divergence due to mean meridional circulation"
mmc[:,:] = div_vm_mmc_vint

se = file_se.createVariable('divse', 'f4', ("time","lat"))
se.units = "W m**-2"
se.long_name = "vertically integrated moist static energy flux divergence due to stationary eddies"
se[:,:] = div_vm_se_vint

file_mmc.close()
file_se.close()
