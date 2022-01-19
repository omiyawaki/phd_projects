import sys
import numpy as np
from scipy import interpolate, integrate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps and mse files
path_ps = sys.argv[1]
path_ua = sys.argv[2]
path_va = sys.argv[3]
path_eke = sys.argv[4]
path_veke = sys.argv[5]

# open files
print('Loading ps file...')
file_ps = Dataset(path_ps, 'r')
print('Done.\n')

print('Loading ua file...')
file_ua = Dataset(path_ua, 'r')
print('Done.\n')

print('Loading va file...')
file_va = Dataset(path_va, 'r')
print('Done.\n')

# read data
ps = file_ps.variables['ps'][:] # (mon x lat x lon)
ua = file_ua.variables['ua'][:] # (day x lev x lat x lon)
va = file_va.variables['va'][:] # (day x lev x lat x lon)

latps = file_ps.variables['lat'][:] # (lat)
lat = file_va.variables['lat'][:] # (lat)
plev = file_va.variables['plev'][:] # (lev)

ps = ps.filled(fill_value=np.nan)
ua = ua.filled(fill_value=np.nan)
va = va.filled(fill_value=np.nan)

# check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
if not np.array_equal(latps,lat):
    latps = latps.filled(fill_value=np.nan)
    lat = lat.filled(fill_value=np.nan)
    f = interpolate.interp1d(latps, ps, axis=1)
    ps = f(lat)

rlat = np.radians(lat)
clat = np.cos(rlat)

time = file_ps.variables['time'][:] 
time_units = file_ps.variables['time'].units
time_cal = file_ps.variables['time'].calendar
ps_date = num2date(time, units=time_units, calendar=time_cal)

time = file_ua.variables['time'][:] 
time_units = file_ua.variables['time'].units
time_cal = file_ua.variables['time'].calendar
daily_date = num2date(time, units=time_units, calendar=time_cal)

# for datasets that fill data below surface as missing data, fill with nans
psd = np.empty([ua.shape[0], ps.shape[1], ps.shape[2]])
# extend ps to daily
for itime in tqdm(range(psd.shape[0])):
    year=daily_date[itime].year
    mon=daily_date[itime].month
    idx_pstime = np.where([ (d.year == year) & (d.month == mon) for d in ps_date])
    psd[itime,...] = ps[idx_pstime[0],...]

# mask out data below surface
pa = np.transpose(np.tile(plev, [psd.shape[0], psd.shape[1], psd.shape[2], 1]), [0,3,1,2])
subsurf = pa > np.transpose(np.tile(psd, [len(plev),1,1,1]),[1,0,2,3])
ua[subsurf] = 0
va[subsurf] = 0

# take monthly mean
# pa_t = np.empty([ps.shape[0], pa.shape[1], pa.shape[2], pa.shape[3]])
# ua_t = np.empty([ps.shape[0], ua.shape[1], ua.shape[2], ua.shape[3]])
# va_t = np.empty([ps.shape[0], va.shape[1], va.shape[2], va.shape[3]])
# ua2_t = np.empty([ps.shape[0], ua.shape[1], ua.shape[2], ua.shape[3]])
# va2_t = np.empty([ps.shape[0], va.shape[1], va.shape[2], va.shape[3]])
res_raw = np.empty([ps.shape[0], ua.shape[1], ua.shape[2]])
tke_raw = np.empty([ps.shape[0], ua.shape[1], ua.shape[2]])
mke_raw = np.empty([ps.shape[0], ua.shape[1], ua.shape[2]])
eke_raw = np.empty([ps.shape[0], ua.shape[1], ua.shape[2]])
teke_raw = np.empty([ps.shape[0], ua.shape[1], ua.shape[2]])
seke_raw = np.empty([ps.shape[0], ua.shape[1], ua.shape[2]])
reke_raw = np.empty([ps.shape[0], ua.shape[1], ua.shape[2]])
for itime in tqdm(range(ps.shape[0])):
    year=ps_date[itime].year
    mon=ps_date[itime].month
    idx_dailytime = np.where([ (d.year == year) & (d.month == mon) for d in daily_date])

    # take monthly mean
    pa_t = np.nanmean(np.squeeze(pa[idx_dailytime[0],...]), axis=0)
    ua_t = np.nanmean(np.squeeze(ua[idx_dailytime[0],...]), axis=0)
    va_t = np.nanmean(np.squeeze(va[idx_dailytime[0],...]), axis=0)
    ua2_t = np.nanmean(np.squeeze(ua[idx_dailytime[0],...]**2), axis=0)
    va2_t = np.nanmean(np.squeeze(va[idx_dailytime[0],...]**2), axis=0)

    # pa_t[itime,...] = np.nanmean(np.squeeze(pa[idx_dailytime[0],...]), axis=0)
    # ua_t[itime,...] = np.nanmean(np.squeeze(ua[idx_dailytime[0],...]), axis=0)
    # va_t[itime,...] = np.nanmean(np.squeeze(va[idx_dailytime[0],...]), axis=0)
    # ua2_t[itime,...] = np.nanmean(np.squeeze(ua[idx_dailytime[0],...]**2), axis=0)
    # va2_t[itime,...] = np.nanmean(np.squeeze(va[idx_dailytime[0],...]**2), axis=0)

    # compute total kinetic energy
    tke_raw[itime,...] = 1/2 * np.nanmean( ua2_t + va2_t, axis=2 )

    # compute mean meridional kinetic energy
    # mke_raw[itime,...] = 1/2 * np.nanmean( ua_t**2 + va_t**2, axis=2 )
    mke_raw[itime,...] = 1/2 * np.nanmean( ua_t, axis=2)**2 + np.nanmean(va_t, axis=2 )**2

    # infer eddy kinetic energy as the residual
    reke_raw[itime,...] = tke_raw[itime,...] - mke_raw[itime,...]

    # deviation from monthly mean
    ua_dt = ua[idx_dailytime[0],...] - np.expand_dims(ua_t, axis=0)
    va_dt = va[idx_dailytime[0],...] - np.expand_dims(va_t, axis=0)

    # compute transient eddy kinetic energy
    teke_raw[itime,...] = 1/2 * np.nanmean( np.nanmean(ua_dt**2, axis=0) + np.nanmean(va_dt**2, axis=0), axis=2)

    # zonal mean of monthly mean
    ua_tz = np.nanmean(ua_t, axis=2)
    va_tz = np.nanmean(va_t, axis=2)

    # zonal deviation of monthly mean
    ua_tdz = ua_t - np.expand_dims(ua_tz, axis=2)
    va_tdz = va_t - np.expand_dims(va_tz, axis=2)

    # compute stationary eddy kinetic energy
    seke_raw[itime,...] = 1/2 * np.nanmean( ua_tdz**2 + va_tdz**2, axis=2)

    # compute transient + stationary eddy kinetic energy
    eke_raw[itime,...] = teke_raw[itime,...] + seke_raw[itime,...]

    # compute residual
    res_raw[itime,...] = tke_raw[itime,...] - ( teke_raw[itime,...] + seke_raw[itime,...] + mke_raw[itime,...] )

    # # deviation from monthly mean zonal mean
    # va_dzt = va_dt - np.nanmean(va_dt, axis=3)[...,np.newaxis]

# pa_tz = np.nanmean(pa_t, axis=3)

# compute vertical integral of kinetic energies (simplified)
if plev[1]-plev[0]>0: # if pressure increases with index
    eke_vint = 1/g*np.trapz(eke_raw, plev, axis=1)
    teke_vint = 1/g*np.trapz(teke_raw, plev, axis=1)
    seke_vint = 1/g*np.trapz(seke_raw, plev, axis=1)
    reke_vint = 1/g*np.trapz(reke_raw, plev, axis=1)
    mke_vint = 1/g*np.trapz(mke_raw, plev, axis=1)
    tke_vint = 1/g*np.trapz(tke_raw, plev, axis=1)
    res_vint = 1/g*np.trapz(res_raw, plev, axis=1)
else:
    eke_vint = -1/g*np.trapz(eke_raw, plev, axis=1)
    teke_vint = -1/g*np.trapz(teke_raw, plev, axis=1)
    seke_vint = -1/g*np.trapz(seke_raw, plev, axis=1)
    reke_vint = -1/g*np.trapz(reke_raw, plev, axis=1)
    mke_vint = -1/g*np.trapz(mke_raw, plev, axis=1)
    tke_vint = -1/g*np.trapz(tke_raw, plev, axis=1)
    res_vint = -1/g*np.trapz(res_raw, plev, axis=1)

# mean
eke_mean =  np.trapz(eke_raw, plev, axis=1)   /np.trapz(np.ones_like(eke_raw), plev, axis=1)
teke_mean = np.trapz(teke_raw, plev, axis=1) /np.trapz(np.ones_like(eke_raw), plev, axis=1)
seke_mean = np.trapz(seke_raw, plev, axis=1) /np.trapz(np.ones_like(eke_raw), plev, axis=1)
reke_mean = np.trapz(reke_raw, plev, axis=1) /np.trapz(np.ones_like(eke_raw), plev, axis=1)
mke_mean =  np.trapz(mke_raw, plev, axis=1)   /np.trapz(np.ones_like(eke_raw), plev, axis=1)
tke_mean =  np.trapz(tke_raw, plev, axis=1)   /np.trapz(np.ones_like(eke_raw), plev, axis=1)
res_mean =  np.trapz(res_raw, plev, axis=1)   /np.trapz(np.ones_like(eke_raw), plev, axis=1)

# save file as netCDF
file_eke = Dataset(path_eke, "w", format='NETCDF4_CLASSIC')
file_veke = Dataset(path_veke, "w", format='NETCDF4_CLASSIC')

# copy attributes from psfile
file_eke.setncatts(file_ps.__dict__)
file_veke.setncatts(file_ps.__dict__)

# copy dimensions from ps file
for name, dimension in file_ps.dimensions.items():
    # if any(name in s for s in ['lon']):
    #     continue
    file_eke.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_veke.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_ps.variables.items():
    if any(name in s for s in ['ps']):
        continue
    
    x = file_eke.createVariable(name, variable.datatype, variable.dimensions)
    file_eke[name].setncatts(file_ps[name].__dict__)
    file_eke[name][:] = file_ps[name][:]
    
    x = file_veke.createVariable(name, variable.datatype, variable.dimensions)
    file_veke[name].setncatts(file_ps[name].__dict__)
    file_veke[name][:] = file_ps[name][:]
    
eke = file_eke.createVariable('eke', 'f4', ("time","lat"))
eke.units = "m**2 s**-2"
eke.long_name = "zonal mean vertical eddy kinetic energy"
eke[:,:] = eke_mean

teke = file_eke.createVariable('teke', 'f4', ("time","lat"))
teke.units = "m**2 s**-2"
teke.long_name = "zonal mean vertical mean transient eddy kinetic energy"
teke[:,:] = teke_mean

seke = file_eke.createVariable('seke', 'f4', ("time","lat"))
seke.units = "m**2 s**-2"
seke.long_name = "zonal mean vertical mean stationary eddy kinetic energy"
seke[:,:] = seke_mean

reke = file_eke.createVariable('reke', 'f4', ("time","lat"))
reke.units = "m**2 s**-2"
reke.long_name = "zonal mean vertical mean eddy kinetic energy (inferred)"
reke[:,:] = reke_mean

mke = file_eke.createVariable('mke', 'f4', ("time","lat"))
mke.units = "m**2 s**-2"
mke.long_name = "zonal mean vertical mean, mean meridional kinetic energy"
mke[:,:] = mke_mean

tke = file_eke.createVariable('tke', 'f4', ("time","lat"))
tke.units = "m**2 s**-2"
tke.long_name = "zonal mean vertical mean total kinetic energy"
tke[:,:] = tke_mean

res = file_eke.createVariable('res', 'f4', ("time","lat"))
res.units = "m**2 s**-2"
res.long_name = "zonal mean vertical mean residual kinetic energy"
res[:,:] = res_mean

file_eke.close()

eke = file_veke.createVariable('eke', 'f4', ("time","lat"))
eke.units = "J m**-2"
eke.long_name = "vertically integrated eddy kinetic energy"
eke[:,:] = eke_vint

teke = file_veke.createVariable('teke', 'f4', ("time","lat"))
teke.units = "J m**-2"
teke.long_name = "vertically integrated transient eddy kinetic energy"
teke[:,:] = teke_vint

seke = file_veke.createVariable('seke', 'f4', ("time","lat"))
seke.units = "J m**-2"
seke.long_name = "vertically integrated stationary eddy kinetic energy"
seke[:,:] = seke_vint

reke = file_veke.createVariable('reke', 'f4', ("time","lat"))
reke.units = "J m**-2"
reke.long_name = "vertically integrated eddy kinetic energy (inferred)"
reke[:,:] = reke_vint

mke = file_veke.createVariable('mke', 'f4', ("time","lat"))
mke.units = "J m**-2"
mke.long_name = "vertically integrated mean meridional kinetic energy"
mke[:,:] = mke_vint

tke = file_veke.createVariable('tke', 'f4', ("time","lat"))
tke.units = "J m**-2"
tke.long_name = "vertically integrated total kinetic energy"
tke[:,:] = tke_vint

res = file_veke.createVariable('res', 'f4', ("time","lat"))
res.units = "J m**-2"
res.long_name = "vertically integrated residual kinetic energy"
res[:,:] = res_vint

file_veke.close()
