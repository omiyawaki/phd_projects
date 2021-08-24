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
path_mses = sys.argv[2]
path_mse = sys.argv[3]
path_mmc = sys.argv[4]
path_se = sys.argv[5]
path_te = sys.argv[6]
path_vas = sys.argv[7]
path_ps = sys.argv[8]

# open files
file_atm = Dataset(path_atm, 'r')
file_mses = Dataset(path_mses, 'r')
file_mse = Dataset(path_mse, 'r')
file_vas = Dataset(path_vas, 'r')
file_ps = Dataset(path_ps, 'r')

# read data
ps = file_ps.variables['aps'][:] # (mon x lat x lon)
va = file_atm.variables['v'][:] # (mon x lev x lat x lon)
mse = file_mse.variables['mse'][:] # (mon x lev x lat x lon)
mses = file_mses.variables['mses'][:] # (mon x lat x lon)
vas = file_vas.variables['vas'][:] # (mon x lat x lon)
lat = file_atm.variables['lat'][:] # (lat)
plev = file_atm.variables['lev'][:] # (lev)

# # for datasets that fill data below surface as missing data, fill with nans
ps = ps.filled(fill_value=np.nan)
va = va.filled(fill_value=np.nan)
mse = mse.filled(fill_value=np.nan)
vas = vas.filled(fill_value=np.nan)
mses = mses.filled(fill_value=np.nan)

# # check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
# lat2d = file_ps.variables['lat'][:] # (mon x lat x lon)
# lat3d = file_mse.variables['lat'][:] # (mon x lev x lat x lon)
# if not np.array_equal(lat2d,lat3d):
#     filledps = ps.filled(fill_value=np.nan)
#     filledlat2d = lat2d.filled(fill_value=np.nan)
#     filledlat3d = lat3d.filled(fill_value=np.nan)
#     f = interpolate.interp1d(filledlat2d, filledps, axis=1)
#     ps = f(filledlat3d)
#     filledps = None; f = None;

# Zonal and annual means of surface pressure
ps_z = np.mean(ps, axis=2)
ps_a = np.mean(ps, axis=(0))
ps_za = np.mean(ps, axis=(0,2))

ps_a_tile = np.tile(ps_a, [ps.shape[0], 1, 1])

# mask out data below surface
pa = np.transpose(np.tile(plev, [ps.shape[0], ps.shape[1], ps.shape[2], 1]), [0,3,1,2])
subsurf = pa > np.transpose(np.tile(ps_a_tile, [len(plev),1,1,1]),[1,0,2,3])
va[subsurf] = np.nan
mse[subsurf] = np.nan
pa[subsurf] = np.nan

# take monthly mean
va_t = np.nanmean(va, axis=0)
mse_t = np.nanmean(mse, axis=0)
vas_t = np.nanmean(vas, axis=0)
mses_t = np.nanmean(mses, axis=0)
pa_t = np.nanmean(pa, axis=0)

# take zonal mean
pa_zt = np.mean(pa_t, axis=2)
va_zt = np.mean(va_t, axis=2)
mse_zt = np.mean(mse_t, axis=2)
vas_zt = np.mean(vas_t, axis=1)
mses_zt = np.mean(mses_t, axis=1)

# compute stationary eddy transport
vms_se = np.squeeze(np.mean((vas_t - vas_zt[...,np.newaxis]) * (mses_t - mses_zt[...,np.newaxis]), axis=1))
vm_se = np.squeeze(np.mean((va_t - va_zt[...,np.newaxis]) * (mse_t - mse_zt[...,np.newaxis]), axis=2))

# deviation from monthly mean zonal mean
va_dt = va - va_t
mse_dt = mse - mse_t
vas_dt = vas - vas_t
mses_dt = mses - mses_t
va_dzt = va_dt - np.mean(va_dt, axis=3)[...,np.newaxis]
mse_dzt = mse_dt - np.mean(mse_dt, axis=3)[...,np.newaxis]
vas_dzt = vas_dt - np.mean(vas_dt, axis=2)[...,np.newaxis]
mses_dzt = mses_dt - np.mean(mses_dt, axis=2)[...,np.newaxis]

# compute transient eddy transport
vms_te = np.mean(vas_dzt * mses_dzt, axis=(0,2))
vm_te = np.mean(va_dzt * mse_dzt, axis=(0,3))

rlat = np.radians(lat)
clat = np.cos(rlat)

# compute vertical integral of heat transports
vm_mmc_vint = np.empty_like(lat)
vm_se_vint = np.empty_like(lat)
for ilat in range(len(lat)):
    ps_z_local = ps_za[ilat]
    vas_z_local = vas_zt[ilat]
    mses_z_local = mses_zt[ilat]
    vms_se_z_local = vms_se[ilat]

    # # remove subsurface data
    # abovesurf = (plev < ps_za[ilat])
    # plev_local = plev[abovesurf]
    # va_z_local = va_z[itime, abovesurf, ilat]
    # mse_z_local = mse_z[itime, abovesurf, ilat]
    # vm_se_local = vm_se[itime, abovesurf, ilat]

    plev_local = pa_zt[:, ilat]
    va_z_local = va_zt[ :, ilat]
    mse_z_local = mse_zt[:, ilat]
    vm_se_local = vm_se[:, ilat]

    plev_local = plev_local[~np.isnan(plev_local)]
    va_z_local = va_z_local[~np.isnan(va_z_local)]
    mse_z_local = mse_z_local[~np.isnan(mse_z_local)]
    vm_se_local = vm_se_local[~np.isnan(vm_se_local)]

    if plev[1]-plev[0]>0: # if pressure increases with index
        # plev_local = np.append(plev_local, ps_za[ilat]) 
        plev_local = np.append(plev_local, ps_z_local) 
        va_z_local = np.append(va_z_local, vas_z_local)
        mse_z_local = np.append(mse_z_local, mses_z_local)
        vm_se_local = np.append(vm_se_local, vms_se_z_local)
    else:
        # plev_local = np.insert(plev_local, 0, ps_za[ilat]) 
        plev_local = np.insert(plev_local, 0, ps_z_local) 
        va_z_local = np.insert(va_z_local, 0, vas_z_local)
        mse_z_local = np.insert(mse_z_local, 0, mses_z_local)
        vm_se_local = np.insert(vm_se_local, 0, vms_se_z_local)

    # va_zdv_local = va_z_local - np.trapz(va_z_local, plev_local)/np.trapz(np.ones_like(va_z_local), plev_local)
    # mse_zdv_local = mse_z_local - np.trapz(mse_z_local, plev_local)/np.trapz(np.ones_like(mse_z_local), plev_local)

    # vm_mmc_vint[itime, ilat] = 2*np.pi*a*clat[ilat]/g*np.trapz(va_zdv_local*mse_zdv_local, plev_local)
    # vm_se_vint[itime, ilat] = 2*np.pi*a*clat[ilat]/g*np.trapz(vm_se_local, plev_local)

    # conservation of mass correction terms
    va_zmc_local = np.trapz(va_z_local, plev_local)/np.trapz(np.ones_like(va_z_local), plev_local)
    mse_zmc_local = np.trapz(mse_z_local, plev_local)/np.trapz(np.ones_like(mse_z_local), plev_local)

    vm_mmc_vint[ilat] = 2*np.pi*a*clat[ilat]/g* (np.trapz(va_z_local*mse_z_local - va_zmc_local*mse_zmc_local, plev_local))
    vm_se_vint[ilat] = 2*np.pi*a*clat[ilat]/g*np.trapz(vm_se_local, plev_local)

    if plev[1]-plev[0]<0: # if pressure increases with index
        vm_mmc_vint[ilat] = -vm_mmc_vint[ilat]
        vm_se_vint[ilat] = -vm_se_vint[ilat]

# save file as netCDF
file_mmc = Dataset(path_mmc, "w", format='NETCDF4_CLASSIC')
file_se = Dataset(path_se, "w", format='NETCDF4_CLASSIC')
file_te = Dataset(path_te, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_mmc.setncatts(file_ps.__dict__)
file_se.setncatts(file_ps.__dict__)
file_te.setncatts(file_ps.__dict__)

# copy dimensions from mse file
for name, dimension in file_mse.dimensions.items():
    if any(name in s for s in ['plev']):
        continue
    file_mmc.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_se.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_te.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
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
    
    x = file_te.createVariable(name, variable.datatype, variable.dimensions)
    file_te[name].setncatts(file_mse[name].__dict__)
    file_te[name][:] = file_mse[name][:]

mmc = file_mmc.createVariable('vmmmc', 'f4', ("time","lat"))
mmc.units = "W"
mmc.long_name = "vertically integrated moist static energy flux transport due to mean meridional circulation"
mmc[:,:] = vm_mmc_vint

se = file_se.createVariable('vmse', 'f4', ("time","lat"))
se.units = "W"
se.long_name = "vertically integrated moist static energy flux transport due to stationary eddies"
se[:,:] = vm_se_vint

te = file_te.createVariable('vmte', 'f4', ("time","lat"))
te.units = "W"
te.long_name = "vertically integrated moist static energy flux transport due to transient eddies"
te[:,:] = vm_te_vint

file_mmc.close()
file_se.close()
file_te.close()
