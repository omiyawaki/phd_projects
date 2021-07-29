import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;
pref = np.linspace(1e3,1e5,50) # higher resolution pressure grid

# paths to ps and mse files
path_ps = sys.argv[1]
path_va = sys.argv[2]
path_vas = sys.argv[3]
path_mse = sys.argv[4]
path_mses = sys.argv[5]
path_vah = sys.argv[6]
path_mseh = sys.argv[7]

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

            f_va = interpolate.interp1d(plev_local, va_local, kind='linear', bounds_error=False)
            f_mse = interpolate.interp1d(plev_local, mse_local, kind='linear', bounds_error=False)

            va_hi[itime,:,ilat,ilon] = f_va(pref)
            mse_hi[itime,:,ilat,ilon] = f_mse(pref)

# save file as netCDF
file_vah = Dataset(path_vah, "w", format='NETCDF4_CLASSIC')
file_mseh = Dataset(path_mseh, "w", format='NETCDF4_CLASSIC')

# copy dimensions from mse file
for name, dimension in file_mse.dimensions.items():
    if any(name in s for s in ['plev']):
        file_vah.createDimension(name, len(pref))
        file_mseh.createDimension(name, len(pref))
    else:
        file_vah.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
        file_mseh.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables from mse file
for name, variable in file_mse.variables.items():
    if any(name in s for s in ['mse']):
        continue

    x = file_vah.createVariable(name, variable.datatype, variable.dimensions)
    file_vah[name].setncatts(file_mse[name].__dict__)
    if name=='plev':
        file_vah[name][:] = pref
    else:
        file_vah[name][:] = file_mse[name][:]

    x = file_mseh.createVariable(name, variable.datatype, variable.dimensions)
    file_mseh[name].setncatts(file_mse[name].__dict__)
    if name=='plev':
        file_mseh[name][:] = pref
    else:
        file_mseh[name][:] = file_mse[name][:]
    
vah = file_vah.createVariable('vah', 'f4', ("time","plev","lat","lon"))
vah.units = "m/s"
vah.long_name = "northward velocity"
vah[:,:,:,:] = va_hi

mseh = file_mseh.createVariable('mseh', 'f4', ("time","plev","lat","lon"))
mseh.units = "J/kg"
mseh.long_name = "moist static energy"
mseh[:,:,:,:] = mse_hi

file_vah.close()
file_mseh.close()
