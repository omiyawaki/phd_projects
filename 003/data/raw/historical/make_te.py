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
path_va = sys.argv[2]
path_ta = sys.argv[3]
path_hus = sys.argv[4]
path_zg = sys.argv[5]
path_vmte = sys.argv[6]

# open files
print('Loading ps file...')
file_ps = Dataset(path_ps, 'r')
print('Done.\n')

print('Loading va file...')
file_va = Dataset(path_va, 'r')
print('Done.\n')

print('Loading ta file...')
file_ta = Dataset(path_ta, 'r')
print('Done.\n')

print('Loading hus file...')
file_hus = Dataset(path_hus, 'r')
print('Done.\n')

print('Loading zg file...')
file_zg = Dataset(path_zg, 'r')
print('Done.\n')

# read data
ps = file_ps.variables['ps'][:] # (mon x lat x lon)
va = file_va.variables['va'][:] # (mon x lev x lat x lon)
ta = file_ta.variables['ta'][:] # (mon x lev x lat x lon)
hus = file_hus.variables['hus'][:] # (mon x lev x lat x lon)
zg = file_zg.variables['zg'][:] # (mon x lev x lat x lon)

latps = file_ps.variables['lat'][:] # (lat)
lat = file_va.variables['lat'][:] # (lat)
plev = file_va.variables['plev'][:] # (lev)

ps = ps.filled(fill_value=np.nan)
ta = ta.filled(fill_value=np.nan)
hus = hus.filled(fill_value=np.nan)
zg = zg.filled(fill_value=np.nan)

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

time = file_hus.variables['time'][:] 
time_units = file_hus.variables['time'].units
time_cal = file_hus.variables['time'].calendar
daily_date = num2date(time, units=time_units, calendar=time_cal)

# for datasets that fill data below surface as missing data, fill with nans
psd = np.empty([ta.shape[0], ps.shape[1], ps.shape[2]])
# extend ps to daily
for itime in tqdm(range(psd.shape[0])):
    year=daily_date[itime].year
    mon=daily_date[itime].month
    idx_pstime = np.where([ (d.year == year) & (d.month == mon) for d in ps_date])
    psd[itime,...] = ps[idx_pstime[0],...]

# mask out data below surface
pa = np.transpose(np.tile(plev, [psd.shape[0], psd.shape[1], psd.shape[2], 1]), [0,3,1,2])
subsurf = pa > np.transpose(np.tile(psd, [len(plev),1,1,1]),[1,0,2,3])
va[subsurf] = np.nan
ta[subsurf] = np.nan
zg[subsurf] = np.nan
hus[subsurf] = np.nan

# compute mse
mse = cpd*ta + L*hus + g*zg

# take monthly mean
pa_t = np.empty([ps.shape[0], pa.shape[1], pa.shape[2], pa.shape[3]])
va_t = np.empty([ps.shape[0], va.shape[1], va.shape[2], va.shape[3]])
mse_t = np.empty([ps.shape[0], mse.shape[1], mse.shape[2], mse.shape[3]])
vm_te = np.empty([ps.shape[0], mse.shape[1], mse.shape[2]])
for itime in tqdm(range(ps.shape[0])):
    year=ps_date[itime].year
    mon=ps_date[itime].month
    idx_dailytime = np.where([ (d.year == year) & (d.month == mon) for d in daily_date])
    # take monthly mean
    pa_t[itime,...] = np.nanmean(np.squeeze(pa[idx_dailytime[0],...]), axis=0)
    va_t[itime,...] = np.nanmean(np.squeeze(va[idx_dailytime[0],...]), axis=0)
    mse_t[itime,...] = np.nanmean(np.squeeze(mse[idx_dailytime[0],...]), axis=0)
    # deviation from monthly mean
    va_dt = np.nanmean(np.squeeze(va[idx_dailytime[0],...]), axis=0) - va[idx_dailytime[0],...]
    mse_dt = np.nanmean(np.squeeze(mse[idx_dailytime[0],...]), axis=0) - mse[idx_dailytime[0],...]
    # deviation from monthly mean zonal mean
    va_dzt = va_dt - np.nanmean(va_dt, axis=3)[...,np.newaxis]
    mse_dzt = mse_dt - np.nanmean(mse_dt, axis=3)[...,np.newaxis]
    # compute zonal mean transient eddy transport
    vm_te[itime,...] = np.nanmean(va_dzt * mse_dzt, axis=(0,3))

pa_tz = np.nanmean(pa_t, axis=3)

# compute vertical integral of heat transports (simplified)
if plev[1]-plev[0]>0: # if pressure increases with index
    vm_te_vint = 2*np.pi*a*clat/g*np.trapz(vm_te, plev, axis=1)
else:
    vm_te_vint = -2*np.pi*a*clat/g*np.trapz(vm_te, plev, axis=1)

# save file as netCDF
file_vmte = Dataset(path_vmte, "w", format='NETCDF4_CLASSIC')

# copy attributes from psfile
file_vmte.setncatts(file_ps.__dict__)

# copy dimensions from ps file
for name, dimension in file_ps.dimensions.items():
    # if any(name in s for s in ['lon']):
    #     continue
    file_vmte.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_ps.variables.items():
    if any(name in s for s in ['ps']):
        continue
    
    x = file_vmte.createVariable(name, variable.datatype, variable.dimensions)
    file_vmte[name].setncatts(file_ps[name].__dict__)
    file_vmte[name][:] = file_ps[name][:]
    
vmte = file_vmte.createVariable('vmte', 'f4', ("time","lat"))
vmte.units = "W"
vmte.long_name = "vertically integrated moist static energy flux transport due to transient eddies"
vmte[:,:] = vm_te_vint

file_vmte.close()
