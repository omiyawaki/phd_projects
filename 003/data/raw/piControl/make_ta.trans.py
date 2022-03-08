import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps and ta files
path_ps = sys.argv[1]
path_va = sys.argv[2]
path_ta = sys.argv[3]
path_mmc = sys.argv[4]
path_se = sys.argv[5]

# open files
file_ps = Dataset(path_ps, 'r')
file_va = Dataset(path_va, 'r')
file_ta = Dataset(path_ta, 'r')

# read data
ps = file_ps.variables['ps'][:] # (mon x lat x lon)
va = file_va.variables['va'][:] # (mon x lev x lat x lon)
ta = file_ta.variables['ta'][:] # (mon x lev x lat x lon)
plev = file_ta.variables['plev'][:]

# for datasets that fill data below surface as missing data, fill with nans
va = va.filled(fill_value=np.nan)
ta = ta.filled(fill_value=np.nan)

# check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
lat2d = file_ps.variables['lat'][:] # (mon x lat x lon)
lat3d = file_ta.variables['lat'][:] # (mon x lev x lat x lon)
if not np.array_equal(lat2d,lat3d):
    print('\nInterpolating ps to 3d lat grid...\n')
    filledps = ps.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filledps, axis=1)
    ps = f(filledlat3d)
    filledps = None; f = None;
    print('\nDone.\n')

# Zonal and annual means of surface pressure
ps_z = np.mean(ps, axis=2)
ps_a = np.mean(ps, axis=(0))
ps_za = np.mean(ps, axis=(0,2))

ps_a_tile = np.tile(ps_a, [ps.shape[0], 1, 1])

# mask out data below surface
pa = np.transpose(np.tile(plev, [ps.shape[0], ps.shape[1], ps.shape[2], 1]), [0,3,1,2])
# subsurf = pa > np.transpose(np.tile(ps_a_tile, [len(plev),1,1,1]),[1,0,2,3])
subsurf = pa > np.transpose(np.tile(ps, [len(plev),1,1,1]),[1,0,2,3])
# va[subsurf] = np.nan
# ta[subsurf] = np.nan
# pa[subsurf] = np.nan

# take zonal mean
va_z = np.nanmean(va, axis=3)
ta_z = np.nanmean(ta, axis=3)
pa_z = np.nanmean(pa, axis=3)

# take annual mean
va_az = np.nanmean(va, axis=(0,3))
ta_az = np.nanmean(ta, axis=(0,3))

# compute stationary eddy transport
vta_se = cpd * np.squeeze(np.nanmean((va - va_z[...,np.newaxis]) * (ta - ta_z[...,np.newaxis]), axis=3))
# print(vta_se[0,:,0])

# 21 22 23
# print(va[7,3,21,:])
# print(va[7,3,22,:])
# print(va[7,3,23,:])
# sys.exit()

rlat = np.radians(lat3d)
clat = np.cos(rlat)

# re-mask subsurface data based on annual mean ps (to set the lower integration bound according to pa)
subsurf_a = pa_z > np.tile(ps_za, (pa_z.shape[0],pa_z.shape[1],1))
va_z[subsurf_a] = np.nan
ta_z[subsurf_a] = np.nan

pa_az = np.nanmean(pa, axis=(0,3))
subsurf_az = pa_az > np.tile(ps_za, (pa_z.shape[1],1))
va_az[subsurf_az] = np.nan
ta_az[subsurf_az] = np.nan

# replace nans with zeros in preparation for taking vertical integrals
vta_se[np.isnan(vta_se)] = 0

mask_z = np.ones_like(va_z)
mask_z[np.isnan(va_z)] = 0

va_z[np.isnan(va_z)] = 0
ta_z[np.isnan(ta_z)] = 0

mask_az = np.ones_like(va_az)
mask_az[np.isnan(va_az)] = 0

va_az[np.isnan(va_az)] = 0
ta_az[np.isnan(ta_az)] = 0

# mass correct MMC
va_mc = va_z - np.trapz(va_az, pa_az, axis=0)/np.trapz(mask_az, pa_az, axis=0)
ta_mc = ta_z - np.trapz(ta_az, pa_az, axis=0)/np.trapz(mask_az, pa_az, axis=0)

vta_mmc = cpd * va_mc * ta_mc

# # compute vertical integral of heat transports (simplified)
# if plev[1]-plev[0]>0: # if pressure increases with index
#     vta_se_vint = 2*np.pi*a*clat/g*np.trapz(vta_se, pa_z, axis=1)
#     vta_mmc_vint = 2*np.pi*a*clat/g*np.trapz(vta_mmc, pa_z, axis=1)
# else:
#     vta_se_vint = -2*np.pi*a*clat/g*np.trapz(vta_se, pa_z, axis=1)
#     vta_mmc_vint = -2*np.pi*a*clat/g*np.trapz(vta_mmc, pa_z, axis=1)

vta_se_vint = np.empty([vta_se.shape[0], vta_se.shape[2]])
vta_mmc_vint = np.empty([vta_mmc.shape[0], vta_mmc.shape[2]])
# compute vertical integral of heat transports (lower integration bound is the annual surface pressure)
for ilat in range(ps_z.shape[1]):
    clat_local = clat[ilat]
    ps_local = ps_za[ilat]
    abovesurf = (plev < ps_local)
    plev_local = plev[abovesurf]

    vta_se_local = vta_se[:,abovesurf,ilat]
    vta_se_int = interpolate.interp1d(plev_local, vta_se_local, axis=1, fill_value="extrapolate")
    vta_se_surf = vta_se_int(ps_local)[..., np.newaxis]

    vta_mmc_local = vta_mmc[:,abovesurf,ilat]
    vta_mmc_int = interpolate.interp1d(plev_local, vta_mmc_local, axis=1, fill_value="extrapolate")
    vta_mmc_surf = vta_mmc_int(ps_local)[..., np.newaxis]

    if plev[1]-plev[0]>0: # if pressure increases with index
        plev_local = np.append(plev_local, ps_local) 
        vta_se_local = np.append(vta_se_local, vta_se_surf, axis=1)
        vta_mmc_local = np.append(vta_mmc_local, vta_mmc_surf, axis=1)

        vta_se_vint[:,ilat] = 2*np.pi*a*clat_local/g*np.trapz(vta_se_local, plev_local, axis=1)
        vta_mmc_vint[:,ilat] = 2*np.pi*a*clat_local/g*np.trapz(vta_mmc_local, plev_local, axis=1)
    else:
        plev_local = np.insert(plev_local, 0, ps_local) 
        vta_se_local = np.concatenate((vta_se_surf, vta_se_local), axis=1)
        vta_mmc_local = np.concatenate((vta_mmc_surf, vta_mmc_local), axis=1)

        vta_se_vint[:,ilat] = -2*np.pi*a*clat_local/g*np.trapz(vta_se_local, plev_local, axis=1)
        vta_mmc_vint[:,ilat] = -2*np.pi*a*clat_local/g*np.trapz(vta_mmc_local, plev_local, axis=1)

# 21 22 23

# print(vta_se[7,:,21])
# print(vta_se[7,:,22])
# print(vta_se[7,:,23])

# print(vta_se_vint[7,21])
# print(vta_se_vint[7,22])
# print(vta_se_vint[7,23])
# print(np.around(1e-15*vta_se_vint[7,21:24],3))

# sys.exit()

# save file as netCDF
file_mmc = Dataset(path_mmc, "w", format='NETCDF4_CLASSIC')
file_se = Dataset(path_se, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_mmc.setncatts(file_ps.__dict__)
file_se.setncatts(file_ps.__dict__)

# copy dimensions from ta file
for name, dimension in file_ta.dimensions.items():
    if any(name in s for s in ['plev']):
        continue
    file_mmc.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_se.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_ta.variables.items():
    if any(name in s for s in ['ta' 'plev']):
        continue
    
    x = file_mmc.createVariable(name, variable.datatype, variable.dimensions)
    file_mmc[name].setncatts(file_ta[name].__dict__)
    file_mmc[name][:] = file_ta[name][:]

    x = file_se.createVariable(name, variable.datatype, variable.dimensions)
    file_se[name].setncatts(file_ta[name].__dict__)
    file_se[name][:] = file_ta[name][:]
    
mmc = file_mmc.createVariable('vtammc', 'f4', ("time","lat"))
mmc.units = "W"
mmc.long_name = "vertically integrated moist static energy flux transport due to mean meridional circulation"
mmc[:,:] = vta_mmc_vint

se = file_se.createVariable('vtase', 'f4', ("time","lat"))
se.units = "W"
se.long_name = "vertically integrated moist static energy flux transport due to stationary eddies"
se[:,:] = vta_se_vint

file_mmc.close()
file_se.close()
