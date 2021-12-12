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

# open as regular numpy arrays
file_ps.set_auto_mask(False)
file_va.set_auto_mask(False)
file_vas.set_auto_mask(False)
file_mse.set_auto_mask(False)
file_mses.set_auto_mask(False)

# read data
ps = file_ps.variables['ps'][:] # (mon x lat x lon)
va = file_va.variables['va'][:] # (mon x lev x lat x lon)
vas = file_vas.variables['vas'][:] # (mon x lev x lat x lon)
mse = file_mse.variables['mse'][:] # (mon x lev x lat x lon)
mses = file_mses.variables['mses'][:] # (mon x lev x lat x lon)
plev = file_mse.variables['plev'][:]

# higher res plev to interpolate to
p_std = np.linspace(1.1e5,1e4,51) # make higher res pressure grid to interpolate to for the lower troposphere
p_std = np.concatenate((p_std, plev[plev<1e4]))

# # # for datasets that fill data below surface as missing data, fill with nans
# va = va.filled(fill_value=np.nan)
# mse = mse.filled(fill_value=np.nan)

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

# mask out data below surface
pa = np.transpose(np.tile(plev, [ps.shape[0], ps.shape[1], ps.shape[2], 1]), [0,3,1,2])
# subsurf = pa > np.transpose(np.tile(ps_a_tile, [len(plev),1,1,1]),[1,0,2,3])
subsurf = pa > np.transpose(np.tile(ps, [len(plev),1,1,1]),[1,0,2,3])
va[subsurf] = np.nan
mse[subsurf] = np.nan
# pa[subsurf] = np.nan

# # initialize array to store interpolated values
# va_i = np.empty([va.shape[0], len(p_std), va.shape[2], va.shape[3]])
# mse_i = np.empty([mse.shape[0], len(p_std), mse.shape[2], mse.shape[3]])

# # insert surface data and interpolate to higher res grid
# for imon in tqdm(range(va.shape[0])):
#     for ilat in range(va.shape[2]):
#         for ilon in range(va.shape[3]):
#             va_col = va[imon,:,ilat,ilon]
#             mse_col = mse[imon,:,ilat,ilon]

#             # truncate subsurf data
#             plev_col = plev[~np.isnan(va_col)] 
#             va_col = va_col[~np.isnan(va_col)]
#             mse_col = mse_col[~np.isnan(mse_col)]

#             # insert surface data
#             if plev[1]-plev[0]>0:
#                 va_col = np.append(va_col, vas[imon,ilat,ilon])
#                 mse_col = np.append(mse_col, mses[imon,ilat,ilon])
#                 plev_col = np.append(plev_col, ps[imon,ilat,ilon])
#             else:
#                 va_col = np.insert(va_col, 0, vas[imon,ilat,ilon])
#                 mse_col = np.insert(mse_col, 0, mses[imon,ilat,ilon])
#                 plev_col = np.insert(plev_col, 0, ps[imon,ilat,ilon])

#             fint_va = interpolate.interp1d(plev_col, va_col, bounds_error=False)
#             fint_mse = interpolate.interp1d(plev_col, mse_col, bounds_error=False)

#             va_i[imon,:,ilat,ilon] = fint_va(p_std)
#             mse_i[imon,:,ilat,ilon] = fint_mse(p_std)

# # replace original with interpolated grid data
# va = va_i
# mse = mse_i
# pa = np.transpose(np.tile(p_std, [ps.shape[0], ps.shape[1], ps.shape[2], 1]), [0,3,1,2])

# take zonal mean
va_z = np.nanmean(va, axis=3)
vas_z = np.nanmean(vas, axis=2)
mse_z = np.nanmean(mse, axis=3)
mses_z = np.nanmean(mses, axis=2)
pa_z = np.nanmean(pa, axis=3)

# take annual mean
va_az = np.nanmean(va, axis=(0,3))
mse_az = np.nanmean(mse, axis=(0,3))

# compute stationary eddy transport
# vms_se = np.squeeze(np.nanmean((vas - vas_z[...,np.newaxis]) * (mses - mses_z[...,np.newaxis]), axis=2))
vm_se = np.squeeze(np.nanmean((va - va_z[...,np.newaxis]) * (mse - mse_z[...,np.newaxis]), axis=3))

rlat = np.radians(lat3d)
clat = np.cos(rlat)

# re-mask subsurface data based on annual mean ps (to set the lower integration bound according to pa)
subsurf_a = pa_z > np.tile(ps_za, (pa_z.shape[0],pa_z.shape[1],1))
va_z[subsurf_a] = np.nan
mse_z[subsurf_a] = np.nan

pa_az = np.nanmean(pa, axis=(0,3))
subsurf_az = pa_az > np.tile(ps_za, (pa_z.shape[1],1))
va_az[subsurf_az] = np.nan
mse_az[subsurf_az] = np.nan

# replace nans with zeros in preparation for taking vertical integrals
vm_se[np.isnan(vm_se)] = 0

mask_z = np.ones_like(va_z)
mask_z[np.isnan(va_z)] = 0

va_z[np.isnan(va_z)] = 0
mse_z[np.isnan(mse_z)] = 0

mask_az = np.ones_like(va_az)
mask_az[np.isnan(va_az)] = 0

va_az[np.isnan(va_az)] = 0
mse_az[np.isnan(mse_az)] = 0

# mass correct MMC
va_mc = va_z - np.trapz(va_az, pa_az, axis=0)/np.trapz(mask_az, pa_az, axis=0)
mse_mc = mse_z - np.trapz(mse_az, pa_az, axis=0)/np.trapz(mask_az, pa_az, axis=0)
# va_mc = va_z - np.trapz(va_z, pa_z, axis=0)/np.trapz(mask_z, pa_z, axis=0)
# mse_mc = mse_z - np.trapz(mse_z, pa_z, axis=0)/np.trapz(mask_z, pa_z, axis=0)

vm_mmc = va_mc * mse_mc

# compute vertical integral of heat transports (simplified)
if plev[1]-plev[0]>0: # if pressure increases with index
    vm_se_vint = 2*np.pi*a*clat/g*np.trapz(vm_se, pa_z, axis=1)
    vm_mmc_vint = 2*np.pi*a*clat/g*np.trapz(vm_mmc, pa_z, axis=1)
else:
    vm_se_vint = -2*np.pi*a*clat/g*np.trapz(vm_se, pa_z, axis=1)
    vm_mmc_vint = -2*np.pi*a*clat/g*np.trapz(vm_mmc, pa_z, axis=1)

# # compute vertical integral of heat transports
# vm_mmc_vint = np.empty_like(ps_z)
# vm_se_vint = np.empty_like(ps_z)
# for itime in tqdm(range(ps_z.shape[0])):
#     for ilat in range(ps_z.shape[1]):
#         # ps_z_local = ps_za[ilat]
#         ps_z_local = ps_z[itime, ilat]
#         vas_z_local = vas_z[itime, ilat]
#         mses_z_local = mses_z[itime, ilat]
#         vms_se_z_local = vms_se[itime, ilat]

#         # remove subsurface data
#         abovesurf = (plev < ps_za[ilat])
#         plev_local = plev[abovesurf]
#         va_z_local = va_z[itime, abovesurf, ilat]
#         mse_z_local = mse_z[itime, abovesurf, ilat]
#         vm_se_local = vm_se[itime, abovesurf, ilat]

#         # plev_local = pa_z[itime, :, ilat]
#         # va_z_local = va_z[itime, :, ilat]
#         # mse_z_local = mse_z[itime, :, ilat]
#         # vm_se_local = vm_se[itime, :, ilat]

#         # plev_local = plev_local[~np.isnan(plev_local)]
#         # va_z_local = va_z_local[~np.isnan(va_z_local)]
#         # mse_z_local = mse_z_local[~np.isnan(mse_z_local)]
#         # vm_se_local = vm_se_local[~np.isnan(vm_se_local)]

#         if plev[1]-plev[0]>0: # if pressure increases with index
#             # plev_local = np.append(plev_local, ps_za[ilat]) 
#             plev_local = np.append(plev_local, ps_z_local) 
#             va_z_local = np.append(va_z_local, vas_z_local)
#             mse_z_local = np.append(mse_z_local, mses_z_local)
#             vm_se_local = np.append(vm_se_local, vms_se_z_local)
#         else:
#             # plev_local = np.insert(plev_local, 0, ps_za[ilat]) 
#             plev_local = np.insert(plev_local, 0, ps_z_local) 
#             va_z_local = np.insert(va_z_local, 0, vas_z_local)
#             mse_z_local = np.insert(mse_z_local, 0, mses_z_local)
#             vm_se_local = np.insert(vm_se_local, 0, vms_se_z_local)

#         # va_zdv_local = va_z_local - np.trapz(va_z_local, plev_local)/np.trapz(np.ones_like(va_z_local), plev_local)
#         # mse_zdv_local = mse_z_local - np.trapz(mse_z_local, plev_local)/np.trapz(np.ones_like(mse_z_local), plev_local)

#         # vm_mmc_vint[itime, ilat] = 2*np.pi*a*clat[ilat]/g*np.trapz(va_zdv_local*mse_zdv_local, plev_local)
#         # vm_se_vint[itime, ilat] = 2*np.pi*a*clat[ilat]/g*np.trapz(vm_se_local, plev_local)

#         # conservation of mass correction terms
#         va_zmc_local = np.trapz(va_z_local, plev_local)/np.trapz(np.ones_like(va_z_local), plev_local)
#         mse_zmc_local = np.trapz(mse_z_local, plev_local)/np.trapz(np.ones_like(mse_z_local), plev_local)

#         vm_mmc_vint[itime, ilat] = 2*np.pi*a*clat[ilat]/g* (np.trapz(va_z_local*mse_z_local - va_zmc_local*mse_zmc_local, plev_local))
#         vm_se_vint[itime, ilat] = 2*np.pi*a*clat[ilat]/g*np.trapz(vm_se_local, plev_local)

#         if plev[1]-plev[0]<0: # if pressure increases with index
#             vm_mmc_vint[itime, ilat] = -vm_mmc_vint[itime, ilat]
#             vm_se_vint[itime, ilat] = -vm_se_vint[itime, ilat]

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
    
mmc = file_mmc.createVariable('vmmmc', 'f4', ("time","lat"))
mmc.units = "W"
mmc.long_name = "vertically integrated moist static energy flux transport due to mean meridional circulation"
mmc[:,:] = vm_mmc_vint

se = file_se.createVariable('vmse', 'f4', ("time","lat"))
se.units = "W"
se.long_name = "vertically integrated moist static energy flux transport due to stationary eddies"
se[:,:] = vm_se_vint

file_mmc.close()
file_se.close()































# # compute vertically integrated V and MSE for annual mean (for mass conservation correction)
# va_za = np.nanmean(va_z, axis=0)
# vas_za = np.nanmean(vas_z, axis=0)
# mse_za = np.nanmean(mse_z, axis=0)
# mses_za = np.nanmean(mses_z, axis=0)
# pa_za = np.nanmean(pa_z, axis=0)

# va_za_vint = np.empty([ps_z.shape[1]])
# mse_za_vint = np.empty([ps_z.shape[1]])
# for ilat in range(ps_z.shape[1]):
#     ps_za_local = ps_za[ilat]
#     vas_za_local = vas_za[ilat]
#     mses_za_local = mses_za[ilat]

#     plev_local = pa_za[:, ilat]
#     va_za_local = va_za[:, ilat]
#     mse_za_local = mse_za[:, ilat]

#     plev_local = plev_local[~np.isnan(plev_local)]
#     va_za_local = va_za_local[~np.isnan(va_za_local)]
#     mse_za_local = mse_za_local[~np.isnan(mse_za_local)]

#     if plev[1]-plev[0]>0: # if pressure increases with index
#         plev_local = np.append(plev_local, ps_za_local) 
#         va_za_local = np.append(va_za_local, vas_za_local)
#         mse_za_local = np.append(mse_za_local, mses_za_local)
#     else:
#         plev_local = np.insert(plev_local, 0, ps_za_local) 
#         va_za_local = np.insert(va_za_local, 0, vas_za_local)
#         mse_za_local = np.insert(mse_za_local, 0, mses_za_local)

#     va_za_vint[ilat] = np.trapz(va_za_local, plev_local)/np.trapz(np.ones_like(va_za_local), plev_local)
#     mse_za_vint[ilat] = np.trapz(mse_za_local, plev_local)/np.trapz(np.ones_like(mse_za_local), plev_local)

