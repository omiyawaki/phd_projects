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
path_va = sys.argv[9]

# open files
file_atm = Dataset(path_atm, 'r')
file_mses = Dataset(path_mses, 'r')
file_mse = Dataset(path_mse, 'r')
file_va = Dataset(path_va, 'r')

# read data
ps = file_atm.variables['aps'][:] # (mon x lat x lon)
va = file_va.variables['vasi'][:] # (mon x lev x lat x lon)
mse = file_mse.variables['msesi'][:] # (mon x lev x lat x lon)
lat = file_atm.variables['lat'][:] # (lat)
si =1e-5*file_mse.variables['lev'][:] # (lev)

# compute pressure field
ps = np.tile(ps, [len(si),1,1,1])
pai = si[:,np.newaxis,np.newaxis,np.newaxis]*ps
pa = np.transpose(pai, [1,0,2,3])

# # for datasets that fill data below surface as missing data, fill with nans
va = va.filled(fill_value=np.nan)
mse = mse.filled(fill_value=np.nan)

# Zonal and annual means of surface pressure
ps_z = np.mean(ps, axis=2)
ps_a = np.mean(ps, axis=(0))
ps_za = np.mean(ps, axis=(0,2))

ps_a_tile = np.tile(ps_a, [ps.shape[0], 1, 1])

# take monthly mean
va_t = np.nanmean(va, axis=0)
mse_t = np.nanmean(mse, axis=0)
pa_t = np.nanmean(pa, axis=0)

# take zonal mean of monthly mean
va_zt = np.nanmean(va_t, axis=2)
mse_zt = np.nanmean(mse_t, axis=2)
pa_zt = np.nanmean(pa_t, axis=2)

# compute stationary eddy transport
vm_se = np.nanmean((va_t - va_zt[...,np.newaxis]) * (mse_t - mse_zt[...,np.newaxis]), axis=2)

# deviation from monthly mean zonal mean
va_dt = va - va_t
mse_dt = mse - mse_t
va_dzt = va_dt - np.nanmean(va_dt, axis=3)[...,np.newaxis]
mse_dzt = mse_dt - np.nanmean(mse_dt, axis=3)[...,np.newaxis]

# compute transient eddy transport
vm_te = np.nanmean(va_dzt * mse_dzt, axis=(0,3))

# monthly mean deviation of zonally averaged circulation
va_z = np.nanmean(va, axis=3)
mse_z = np.nanmean(mse, axis=3)
va_zdt = va_zt - va_z
mse_zdt = mse_zt - mse_z

# compute transient overturning circulation transport
vm_toc = np.nanmean(va_zdt * mse_zdt, axis=(0))

# combine TOC with TE
# vm_te = vm_te + vm_toc

rlat = np.radians(lat)
clat = np.cos(rlat)

# compute vertical integral of heat transports
vm_mmc_vint = np.empty_like(lat)
vm_se_vint = np.empty_like(lat)
vm_te_vint = np.empty_like(lat)
for ilat in range(len(lat)):
    plev_local = pa_zt[:, ilat]
    va_zt_local = va_zt[:, ilat]
    mse_zt_local = mse_zt[:, ilat]
    vm_se_local = vm_se[:, ilat]
    vm_te_local = vm_te[:, ilat]

    plev_local = plev_local[~np.isnan(plev_local)]
    va_zt_local = va_zt_local[~np.isnan(va_zt_local)]
    mse_zt_local = mse_zt_local[~np.isnan(mse_zt_local)]
    vm_se_local = vm_se_local[~np.isnan(vm_se_local)]
    vm_te_local = vm_te_local[~np.isnan(vm_te_local)]

    # conservation of mass correction terms
    va_zmc_local = np.trapz(va_zt_local, plev_local)/np.trapz(np.ones_like(va_zt_local), plev_local)
    mse_zmc_local = np.trapz(mse_zt_local, plev_local)/np.trapz(np.ones_like(mse_zt_local), plev_local)

    vm_mmc_vint[ilat] = 2*np.pi*a*clat[ilat]/g* (np.trapz(va_zt_local*mse_zt_local - va_zmc_local*mse_zmc_local, plev_local))
    vm_se_vint[ilat] = 2*np.pi*a*clat[ilat]/g*np.trapz(vm_se_local, plev_local)
    vm_te_vint[ilat] = 2*np.pi*a*clat[ilat]/g*np.trapz(vm_te_local, plev_local)

    if plev_local[1]-plev_local[0]<0: # if pressure increases with index
        vm_mmc_vint[ilat] = -vm_mmc_vint[ilat]
        vm_se_vint[ilat] = -vm_se_vint[ilat]
        vm_te_vint[ilat] = -vm_te_vint[ilat]


# save file as netCDF
file_mmc = Dataset(path_mmc, "w", format='NETCDF4_CLASSIC')
file_se = Dataset(path_se, "w", format='NETCDF4_CLASSIC')
file_te = Dataset(path_te, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_mmc.setncatts(file_mses.__dict__)
file_se.setncatts(file_mses.__dict__)
file_te.setncatts(file_mses.__dict__)

# copy dimensions from mse file
for name, dimension in file_mse.dimensions.items():
    if any(name in s for s in ['lev']):
        continue

    file_mmc.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_se.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_te.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_mse.variables.items():
    if any(name in s for s in ['msesi' 'lev']):
        continue
    
    print(name)
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
