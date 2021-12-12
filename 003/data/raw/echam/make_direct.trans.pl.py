import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# list of 2d vars
vars2d = ['ps', 'vas', 'mses']

# list of 3d vars
vars3d = ['v', 'mse']

# paths to ps and mse files
path_atm = sys.argv[1]
path_mses = sys.argv[2]
path_mse = sys.argv[3]
path_mmc = sys.argv[4]
path_se = sys.argv[5]
path_te = sys.argv[6]
path_vas = sys.argv[7]
path_ps = sys.argv[8]
path_aht = sys.argv[9]
path_mse_a = sys.argv[10]
path_atm_a = sys.argv[11]
path_ps_a = sys.argv[12]
path_mses_a = sys.argv[13]
path_vas_a = sys.argv[14]

# open files
print('Loading ATM file...')
file_atm = Dataset(path_atm, 'r')
print('Done.\n')

print('Loading mses file...')
file_mses = Dataset(path_mses, 'r')
print('Done.\n')

print('Loading mse file...')
file_mse = Dataset(path_mse, 'r')
print('Done.\n')

print('Loading vas file...')
file_vas = Dataset(path_vas, 'r')
print('Done.\n')

print('Loading ps file...')
file_ps = Dataset(path_ps, 'r')
print('Done.\n')

print('Loading aht file...')
file_aht = Dataset(path_aht, 'r')
print('Done.\n')

print('Loading mse_a file...')
file_mse_a = Dataset(path_mse_a, 'r')
print('Done.\n')

print('Loading mse_a file...')
file_atm_a = Dataset(path_atm_a, 'r')
print('Done.\n')

print('Loading ps file...')
file_ps_a = Dataset(path_ps_a, 'r')
print('Done.\n')

print('Loading mses_a file...')
file_mses_a = Dataset(path_mses_a, 'r')
print('Done.\n')

print('Loading vas_a file...')
file_vas_a = Dataset(path_vas_a, 'r')
print('Done.\n')

# read data
ps = file_ps.variables['aps'][:] # (mon x lat x lon)
ps_a = np.squeeze(file_ps_a.variables['aps'][:]) # (mon x lat x lon)
va = file_atm.variables['v'][:] # (mon x lev x lat x lon)
va_a = np.squeeze(file_atm_a.variables['v'][:]) # (mon x lev x lat x lon)
mse = file_mse.variables['mse'][:] # (mon x lev x lat x lon)
mses = file_mses.variables['mses'][:] # (mon x lat x lon)
mse_a = np.squeeze(file_mse_a.variables['mse'][:]) # (mon x lev x lat x lon)
vas = file_vas.variables['vas'][:] # (mon x lat x lon)
mses_a = np.squeeze(file_mses_a.variables['mses'][:]) # (mon x lev x lat x lon)
vas_a = np.squeeze(file_vas_a.variables['vas'][:]) # (mon x lev x lat x lon)
lat = file_atm.variables['lat'][:] # (lat)
plev = file_atm.variables['lev'][:] # (lev)

p_std = np.linspace(1.1e5,1e4,51) # make higher res pressure grid to interpolate to for the lower troposphere
p_std = np.concatenate((p_std, plev[plev<1e4]))

# # for datasets that fill data below surface as missing data, fill with nans
plev = plev.filled(fill_value=np.nan)
ps = ps.filled(fill_value=np.nan)
ps_a = ps_a.filled(fill_value=np.nan)
va = va.filled(fill_value=np.nan)
va_a = va_a.filled(fill_value=np.nan)
mse = mse.filled(fill_value=np.nan)
mse_a = mse_a.filled(fill_value=np.nan)
mses_a = mses_a.filled(fill_value=np.nan)
vas = vas.filled(fill_value=np.nan)
vas_a = vas_a.filled(fill_value=np.nan)
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

# Zonal and monthly means of surface pressure
ps_z = np.mean(ps, axis=2)
ps_t = np.mean(ps, axis=(0))
ps_za = np.mean(ps, axis=(0,2))

ps_t_tile = np.tile(ps_t, [ps.shape[0], 1, 1])
ps_a_tile = np.tile(ps_a, [ps.shape[0], 1, 1])

# mask out data below surface
pa = np.transpose(np.tile(plev, [ps.shape[0], ps.shape[1], ps.shape[2], 1]), [0,3,1,2])

# subsurf = pa > np.transpose(np.tile(ps_a_tile, [len(plev),1,1,1]),[1,0,2,3])
# subsurf = pa > np.transpose(np.tile(ps_t_tile, [len(plev),1,1,1]),[1,0,2,3])
subsurf = pa > np.transpose(np.tile(ps, [len(plev),1,1,1]),[1,0,2,3])

va[subsurf] = np.nan
mse[subsurf] = np.nan
# pa[subsurf] = np.nan

pa_a = np.transpose(np.tile(plev, [ps_a.shape[0], ps_a.shape[1], 1]), [2,0,1])
subsurf = pa_a > np.tile(ps_a, [len(plev),1,1])

va_a[subsurf] = np.nan
mse_a[subsurf] = np.nan
# pa_a[subsurf] = np.nan

# take monthly mean
va_t = np.nanmean(va, axis=0)
mse_t = np.nanmean(mse, axis=0)
vas_t = np.nanmean(vas, axis=0)
mses_t = np.nanmean(mses, axis=0)
pa_t = np.nanmean(pa, axis=0)

# take zonal mean
pa_zt = np.nanmean(pa_t, axis=2)
va_zt = np.nanmean(va_t, axis=2)
mse_zt = np.nanmean(mse_t, axis=2)
vas_zt = np.nanmean(vas_t, axis=1)
mses_zt = np.nanmean(mses_t, axis=1)

ps_az = np.nanmean(ps_a, axis=1)
vas_az = np.nanmean(vas_a, axis=1)
mses_az = np.nanmean(mses_a, axis=1)
pa_az = np.nanmean(pa_a, axis=2)
va_az = np.nanmean(va_a, axis=2)
mse_az = np.nanmean(mse_a, axis=2)

# compute stationary eddy transport
vms_se = np.squeeze(np.nanmean((vas_t - vas_zt[...,np.newaxis]) * (mses_t - mses_zt[...,np.newaxis]), axis=1))
vm_se = np.squeeze(np.nanmean((va_t - va_zt[...,np.newaxis]) * (mse_t - mse_zt[...,np.newaxis]), axis=2))

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

# re-mask subsurface data based on annual mean ps (to set the lower integration bound according to pa)
subsurf_a = pa_az > np.tile(ps_az, (pa_az.shape[0],1))
va_zt[subsurf_a] = np.nan
mse_zt[subsurf_a] = np.nan
va_az[subsurf_a] = np.nan
mse_az[subsurf_a] = np.nan

# replace nans with zeros in preparation for taking vertical integrals
vm_se[np.isnan(vm_se)] = 0

mask_zt = np.ones_like(va_zt)
mask_zt[np.isnan(va_zt)] = 0

va_zt[np.isnan(va_zt)] = 0
mse_zt[np.isnan(mse_zt)] = 0

mask_az = np.ones_like(va_az)
mask_az[np.isnan(va_az)] = 0

va_az[np.isnan(va_az)] = 0
mse_az[np.isnan(mse_az)] = 0

# mass correct MMC
# va_mc = va_zt - np.trapz(va_az, pa_az, axis=0)/np.trapz(mask_az, pa_az, axis=0)
# mse_mc = mse_zt - np.trapz(mse_az, pa_az, axis=0)/np.trapz(mask_az, pa_az, axis=0)
va_mc = va_zt - np.trapz(va_zt, pa_zt, axis=0)/np.trapz(mask_zt, pa_zt, axis=0)
mse_mc = mse_zt - np.trapz(mse_zt, pa_zt, axis=0)/np.trapz(mask_zt, pa_zt, axis=0)

vm_mmc = va_mc * mse_mc

# compute vertical integral of heat transports (simplified)
if plev[1]-plev[0]>0: # if pressure increases with index
    vm_se_vint = 2*np.pi*a*clat/g*np.trapz(vm_se, pa_az, axis=0)
    vm_mmc_vint = 2*np.pi*a*clat/g*np.trapz(vm_mmc, pa_az, axis=0)
else:
    vm_se_vint = -2*np.pi*a*clat/g*np.trapz(vm_se, pa_az, axis=0)
    vm_mmc_vint = -2*np.pi*a*clat/g*np.trapz(vm_mmc, pa_az, axis=0)

# vm_mmc_vint = np.empty_like(lat)
# vm_se_vint = np.empty_like(lat)
# for ilat in tqdm(range(len(lat))):
#     ps_z_local = ps_za[ilat]
#     ps_az_local = ps_az[ilat]
#     mses_z_local = mses_zt[ilat]
#     mses_az_local = mses_az[ilat]
#     vms_se_z_local = vms_se[ilat]

#     vas_z_local = vas_zt[ilat]
#     vas_az_local = vas_az[ilat]

#     plev_local = pa_zt[:, ilat]
#     plev_az_local = pa_az[:, ilat]
#     va_z_local = va_zt[ :, ilat]
#     va_az_local = va_az[ :, ilat]
#     mse_z_local = mse_zt[:, ilat]
#     mse_az_local = mse_az[:, ilat]
#     vm_se_local = vm_se[:, ilat]

#     # remove subsurface data
#     subsurf = plev_local > ps_z_local
#     plev_local[subsurf] = np.nan
#     va_z_local[subsurf] = np.nan
#     mse_z_local[subsurf] = np.nan
#     vm_se_local[subsurf] = np.nan

#     subsurf_a = plev_az_local > ps_az_local
#     plev_az_local[subsurf_a] = np.nan
#     va_az_local[subsurf_a] = np.nan
#     mse_az_local[subsurf_a] = np.nan

#     plev_local = plev_local[~np.isnan(plev_local)]
#     va_z_local = va_z_local[~np.isnan(va_z_local)]
#     mse_z_local = mse_z_local[~np.isnan(mse_z_local)]
#     plev_az_local = plev_az_local[~np.isnan(plev_az_local)]
#     va_az_local = va_az_local[~np.isnan(va_az_local)]
#     mse_az_local = mse_az_local[~np.isnan(mse_az_local)]
#     vm_se_local = vm_se_local[~np.isnan(vm_se_local)]

#     print(plev_az_local)
#     print(va_az_local)
#     print(mse_az_local)

#     # if plev[1]-plev[0]>0: # if pressure increases with index
#     #     # plev_local = np.append(plev_local, ps_za[ilat]) 
#     #     plev_local = np.append(plev_local, ps_z_local) 
#     #     va_z_local = np.append(va_z_local, vas_z_local)
#     #     mse_z_local = np.append(mse_z_local, mses_z_local)
#     #     vm_se_local = np.append(vm_se_local, vms_se_z_local)

#     #     plev_az_local = np.append(plev_az_local, ps_az_local) 
#     #     va_az_local = np.append(va_az_local, vas_az_local)
#     #     mse_az_local = np.append(mse_az_local, mses_az_local)

#     # else:
#     #     # plev_local = np.insert(plev_local, 0, ps_za[ilat]) 
#     #     plev_local = np.insert(plev_local, 0, ps_z_local) 
#     #     va_z_local = np.insert(va_z_local, 0, vas_z_local)
#     #     mse_z_local = np.insert(mse_z_local, 0, mses_z_local)
#     #     vm_se_local = np.insert(vm_se_local, 0, vms_se_z_local)

#     #     plev_az_local = np.insert(plev_az_local, 0, ps_z_local) 
#     #     va_az_local = np.insert(va_az_local, 0, vas_az_local)
#     #     mse_az_local = np.insert(mse_az_local, 0, mses_az_local)

#     plev_local_for_vmse = plev_local # retain this plev for SE use only

#     # # set lower boundary equal to annual mean ps
#     # if ps_az_local >= ps_z_local: # if annual mean ps is below the monthly mean ps, extrapolate
#     #     # mse_fint = interpolate.interp1d(plev_local, mse_z_local, fill_value='extrapolate')
#     #     # va_fint = interpolate.interp1d(plev_local, va_z_local, fill_value='extrapolate')

#     #     # if plev[1]-plev[0]>0: # if pressure increases with index
#     #     #     mse_z_local = np.append(mse_z_local, mse_fint(ps_az_local))
#     #     #     va_z_local = np.append(va_z_local, va_fint(ps_az_local))
#     #     #     plev_local = np.append(plev_local, ps_az_local)
#     #     # else:
#     #     #     mse_z_local = np.insert(mse_z_local, 0, mse_fint(ps_az_local))
#     #     #     va_z_local = np.insert(va_z_local, 0, va_fint(ps_az_local))
#     #     #     plev_local = np.insert(plev_local, 0, ps_az_local)

#     #     print('skip annual mean insert')

#     # else: # else, interpolate and replace the monthly ps datapoint
#     #     mse_fint = interpolate.interp1d(plev_local, mse_z_local)
#     #     va_fint = interpolate.interp1d(plev_local, va_z_local)

#     #     if plev[1]-plev[0]>0: # if pressure increases with index
#     #         mse_z_local[-1] = mse_fint(ps_az_local)
#     #         va_z_local[-1] = va_fint(ps_az_local)
#     #         plev_local[-1] = ps_az_local
#     #     else:
#     #         mse_z_local[0] = mse_fint(ps_az_local)
#     #         va_z_local[0] = va_fint(ps_az_local)
#     #         plev_local[0] = ps_az_local

#     # va_zdv_local = va_z_local - np.trapz(va_z_local, plev_local)/np.trapz(np.ones_like(va_z_local), plev_local)
#     # mse_zdv_local = mse_z_local - np.trapz(mse_z_local, plev_local)/np.trapz(np.ones_like(mse_z_local), plev_local)

#     # vm_mmc_vint[itime, ilat] = 2*np.pi*a*clat[ilat]/g*np.trapz(va_zdv_local*mse_zdv_local, plev_local)
#     # vm_se_vint[itime, ilat] = 2*np.pi*a*clat[ilat]/g*np.trapz(vm_se_local, plev_local)

#     # conservation of mass correction terms
#     # va_zmc_local = np.trapz(va_z_local, plev_local)/np.trapz(np.ones_like(va_z_local), plev_local)
#     # mse_zmc_local = np.trapz(mse_z_local, plev_local)/np.trapz(np.ones_like(mse_z_local), plev_local)
#     va_zmc_local = np.trapz(va_az_local, plev_az_local)/np.trapz(np.ones_like(va_az_local), plev_az_local)
#     mse_zmc_local = np.trapz(mse_az_local, plev_az_local)/np.trapz(np.ones_like(mse_az_local), plev_az_local)

#     vm_mmc_vint[ilat] = 2*np.pi*a*clat[ilat]/g* (np.trapz(va_z_local*mse_z_local - va_zmc_local*mse_zmc_local, plev_local))
#     vm_se_vint[ilat] = 2*np.pi*a*clat[ilat]/g*np.trapz(vm_se_local, plev_local_for_vmse)

#     if plev[1]-plev[0]<0: # if pressure increases with index
#         vm_mmc_vint[ilat] = -vm_mmc_vint[ilat]
#         vm_se_vint[ilat] = -vm_se_vint[ilat]

# save file as netCDF
file_mmc = Dataset(path_mmc, "w", format='NETCDF4_CLASSIC')
file_se = Dataset(path_se, "w", format='NETCDF4_CLASSIC')
file_te = Dataset(path_te, "w", format='NETCDF4_CLASSIC')

# copy attributes from aht file
file_mmc.setncatts(file_aht.__dict__)
file_se.setncatts(file_aht.__dict__)
file_te.setncatts(file_aht.__dict__)

# copy dimensions from mse file
for name, dimension in file_aht.dimensions.items():
    # if any(name in s for s in ['plev']):
    #     continue
    file_mmc.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_se.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_te.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_aht.variables.items():
    # if any(name in s for s in ['mse' 'plev']):
    #     continue
    
    x = file_mmc.createVariable(name, variable.datatype, variable.dimensions)
    file_mmc[name].setncatts(file_aht[name].__dict__)
    file_mmc[name][:] = file_aht[name][:]

    x = file_se.createVariable(name, variable.datatype, variable.dimensions)
    file_se[name].setncatts(file_aht[name].__dict__)
    file_se[name][:] = file_aht[name][:]
    
    x = file_te.createVariable(name, variable.datatype, variable.dimensions)
    file_te[name].setncatts(file_aht[name].__dict__)
    file_te[name][:] = file_aht[name][:]

mmc = file_mmc.createVariable('vmmmc', 'f4', ("time","lat"))
mmc.units = "W"
mmc.long_name = "vertically integrated moist static energy flux transport due to mean meridional circulation"
mmc[:,:] = vm_mmc_vint

se = file_se.createVariable('vmse', 'f4', ("time","lat"))
se.units = "W"
se.long_name = "vertically integrated moist static energy flux transport due to stationary eddies"
se[:,:] = vm_se_vint

# te = file_te.createVariable('vmte', 'f4', ("time","lat"))
# te.units = "W"
# te.long_name = "vertically integrated moist static energy flux transport due to transient eddies"
# te[:,:] = vm_te_vint

file_mmc.close()
file_se.close()
file_te.close()
