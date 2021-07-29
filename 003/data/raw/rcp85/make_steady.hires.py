import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;
pref = np.linspace(0,1e5,50) # higher resolution pressure grid

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

# Zonal and annual means of surface pressure
ps_z = np.mean(ps, axis=2)
ps_a = np.mean(ps, axis=(0))
ps_za = np.mean(ps, axis=(0,2))

ps_a_tile = np.tile(ps_a, [ps.shape[0], 1, 1])

va_hi = np.empty([va.shape[0],len(pref),va.shape[2], va.shape[3]])
mse_hi = np.empty([mse.shape[0],len(pref),mse.shape[2], mse.shape[3]])
# interpolate to higher resolution pressure grid
for itime in tqdm(range(ps.shape[0])):
    for ilat in range(ps.shape[1]):
        for ilon in range(ps.shape[2]):
            # remove subsurface data and insert surface data
            abovesurf = (plev < ps[itime,ilat,ilon])
            plev_local = plev[abovesurf]
            va_local = va[itime, abovesurf, ilat, ilon]
            mse_local = mse[itime, abovesurf, ilat, ilon]
            if plev[1]-plev[0]>0: # if pressure increases with index
                plev_local = np.append(plev_local, ps[itime,ilat,ilon])
                va_local = np.append(va_local, vas[itime,ilat,ilon])
                mse_local = np.append(mse_local, mses[itime,ilat,ilon])
            else:
                plev_local = np.insert(plev_local, 0, ps[itime,ilat,ilon])
                va_local = np.insert(va_local, 0, vas[itime,ilat,ilon])
                mse_local = np.insert(mse_local, 0, mses[itime,ilat,ilon])

            f_va = interpolate.interp1d(plev_local, va_local, kind='cubic', bounds_error=False)
            f_mse = interpolate.interp1d(plev_local, mse_local, kind='cubic', bounds_error=False)

            va_hi[itime,:,ilat,ilon] = f_va(pref)
            mse_hi[itime,:,ilat,ilon] = f_mse(pref)

va = None; mse = None; plev=None;
va = va_hi; mse = mse_hi; plev=pref;
        
# take zonal mean
va_z = np.nanmean(va, axis=3)
mse_z = np.nanmean(mse, axis=3)

# compute stationary eddy transport
vm_se = np.squeeze(np.mean((va - va_z[...,np.newaxis]) * (mse - mse_z[...,np.newaxis]), axis=3))

# # compute MSE flux divergences
rlat = np.radians(lat3d)
clat = np.cos(rlat)

# compute vertical integral
div_mmc_vint = np.empty_like(ps_z)
div_se_vint = np.empty_like(ps_z)
for itime in tqdm(range(ps_z.shape[0])):
    for ilat in range(ps_z.shape[1]):

        # remove subsurface data
        # abovesurf = (plev < ps_z[itime, ilat])

        abovesurf = (plev < ps_za[ilat])
        plev_local = plev[abovesurf]
        va_z_local = va_z[itime, abovesurf, ilat]
        mse_z_local = mse_z[itime, abovesurf, ilat]
        vm_se_local = vm_se[itime, abovesurf, ilat]
        # plev_local = plev
        # va_z_local = va_z[itime, :, ilat]
        # mse_z_local = mse_z[itime, :, ilat]
        # vm_se_local = vm_se[itime, :, ilat]

        # insert surface value
        # plev_local = np.insert(plev_local, 0, ps_z[itime, ilat]) 

        # plev_local = np.insert(plev_local, 0, ps_za[ilat]) 
        # va_z_local = np.insert(va_z_local, 0, vas_z_local)
        # mse_z_local = np.insert(mse_z_local, 0, mses_z_local)
        # vm_se_local = np.insert(vm_se_local, 0, vms_se_local)

        va_zdv_local = va_z_local + np.trapz(va_z_local, plev_local)/np.trapz(np.ones_like(va_z_local), plev_local)
        mse_zdv_local = mse_z_local + np.trapz(mse_z_local, plev_local)/np.trapz(np.ones_like(mse_z_local), plev_local)

        div_mmc_vint[itime, ilat] = 2*np.pi*a*clat[ilat]/g*np.trapz(va_zdv_local*mse_zdv_local, plev_local)
        div_se_vint[itime, ilat] = 2*np.pi*a*clat[ilat]/g*np.trapz(vm_se_local, plev_local)

# div_mmc = None; div_se = None;

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
mmc.units = "W"
mmc.long_name = "vertically integrated moist static energy flux divergence due to mean meridional circulation"
mmc[:,:] = div_mmc_vint

se = file_se.createVariable('divse', 'f4', ("time","lat"))
se.units = "W"
se.long_name = "vertically integrated moist static energy flux divergence due to stationary eddies"
se[:,:] = div_se_vint

file_mmc.close()
file_se.close()
