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
path_vas = sys.argv[4]
path_ps = sys.argv[5]
path_vasi = sys.argv[6]
path_msesi = sys.argv[7]

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
plev = plev.filled(fill_value=np.nan)
si = 1e-5*plev
ps = ps.filled(fill_value=np.nan)
va = va.filled(fill_value=np.nan)
mse = mse.filled(fill_value=np.nan)
vas = vas.filled(fill_value=np.nan)
mses = mses.filled(fill_value=np.nan)

# mask out data below surface
pa = np.transpose(np.tile(plev, [ps.shape[0], ps.shape[1], ps.shape[2], 1]), [0,3,1,2])
subsurf = pa > np.transpose(np.tile(ps, [len(plev),1,1,1]),[1,0,2,3])
va[subsurf] = np.nan
mse[subsurf] = np.nan
pa[subsurf] = np.nan

# interpolate to sigma coord
vasi = np.empty_like(va)
msesi = np.empty_like(mse)
for itime in tqdm(range(mse.shape[0])):
    for ilat in range(mse.shape[2]):
        for ilon in range(mse.shape[3]):
            ps_local = ps[itime, ilat, ilon]
            vas_local = vas[itime, ilat, ilon]
            mses_local = mses[itime, ilat, ilon]

            plev_local = pa[itime, :, ilat, ilon]
            va_local = va[itime, :, ilat, ilon]
            mse_local = mse[itime, :, ilat, ilon]

            plev_local = plev_local[~np.isnan(plev_local)]
            va_local = va_local[~np.isnan(va_local)]
            mse_local = mse_local[~np.isnan(mse_local)]

            if plev[1]-plev[0]>0: # if pressure increases with index
                if not plev_local[-1] == ps_local:
                    plev_local = np.append(plev_local, ps_local) 
                    va_local = np.append(va_local, vas_local)
                    mse_local = np.append(mse_local, mses_local)
            else:
                if not plev_local[0] == ps_local:
                    plev_local = np.insert(plev_local, 0, ps_local) 
                    va_local = np.insert(va_local, 0, vas_local)
                    mse_local = np.insert(mse_local, 0, mses_local)

            si_local = plev_local / ps_local

            f_va = interpolate.interp1d(si_local, va_local, kind='linear', fill_value="extrapolate")
            vasi[itime, :, ilat, ilon] = f_va(si)

            f_mse = interpolate.interp1d(si_local, mse_local, kind='linear', fill_value="extrapolate")
            msesi[itime, :, ilat, ilon] = f_mse(si)

# save file as netCDF
file_vasi = Dataset(path_vasi, "w", format='NETCDF4_CLASSIC')
file_msesi = Dataset(path_msesi, "w", format='NETCDF4_CLASSIC')

# # copy attributes from atm file
# file_vasi.setncatts(file_atm.__dict__)

# copy dimensions from atm file
for name, dimension in file_atm.dimensions.items():
    file_vasi.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_msesi.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables from atm file
for name, variable in file_atm.variables.items():
    if any(name in s for s in ['t', 'u', 'v', 'q', 'aps', 'omega', 'geopoth']):
        continue
    
    x = file_vasi.createVariable(name, variable.datatype, variable.dimensions)
    file_vasi[name].setncatts(file_atm[name].__dict__)
    file_vasi[name][:] = file_atm[name][:]

    x = file_msesi.createVariable(name, variable.datatype, variable.dimensions)
    file_msesi[name].setncatts(file_atm[name].__dict__)
    file_msesi[name][:] = file_atm[name][:]
    
va_si = file_vasi.createVariable('vasi', 'f4', ("time","lev","lat","lon"))
va_si.units = "m/s"
va_si.long_name = "v component velocity"
va_si[:,:,:,:] = vasi

mse_si = file_msesi.createVariable('msesi', 'f4', ("time","lev","lat","lon"))
mse_si.units = "m/s"
mse_si.long_name = "v component velocity"
mse_si[:,:,:,:] = msesi

file_vasi.close()
file_msesi.close()
