import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps, ta, zg, and hus files
path_ps = sys.argv[1]
path_mse = sys.argv[2]
path_mses = sys.argv[3]

# open files
file_ps = Dataset(path_ps, 'r')
file_mse = Dataset(path_mse, 'r')

# read data
ps = file_ps.variables['ps'][:] # (day x lev x lat x lon)
mse = file_mse.variables['mse'][:] # (day x lev x lat x lon)
plev = file_mse.variables['plev'][:]

# for datasets that fill data below surface as missing data, fill with nans
ps = ps.filled(fill_value=np.nan)
mse = mse.filled(fill_value=np.nan)

# Zonal and annual means of surface pressure
ps_z = np.mean(ps, axis=2)
ps_a = np.mean(ps, axis=(0))
ps_za = np.mean(ps, axis=(0,2))

ps_a_tile = np.tile(ps_a, [ps.shape[0], 1, 1])

# mask out data below surface
pa = np.transpose(np.tile(plev, [ps.shape[0], ps.shape[1], ps.shape[2], 1]), [0,3,1,2])
subsurf = pa > np.transpose(np.tile(ps_a_tile, [len(plev),1,1,1]),[1,0,2,3])
mse[subsurf] = np.nan
pa[subsurf] = np.nan

# interpolate MSE at ps
ms = np.empty([mse.shape[0], mse.shape[2], mse.shape[3]])
for itime in tqdm(range(mse.shape[0])):
    for ilat in range(mse.shape[2]):
        for ilon in range(mse.shape[3]):
            plev_local = pa[itime,:,ilat,ilon]
            mse_local = mse[itime,:,ilat,ilon]

            plev_local = plev_local[~np.isnan(plev_local)]
            mse_local = mse_local[~np.isnan(mse_local)]

            f = interpolate.interp1d(plev_local, mse_local, fill_value='extrapolate')
            # f = interpolate.interp1d(plev, mse[itime,:,ilat,ilon], fill_value='extrapolate')
            ms[itime,ilat,ilon] = f(ps[itime,ilat,ilon])

# save file as netCDF
file_mses = Dataset(path_mses, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_mses.setncatts(file_ps.__dict__)

# copy dimensions from mse file
for name, dimension in file_mse.dimensions.items():
    if any(name in s for s in ['plev']):
        continue
    file_mses.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_mse.variables.items():
    if any(name in s for s in ['mse' 'plev']):
        continue
    
    x = file_mses.createVariable(name, variable.datatype, variable.dimensions)
    file_mses[name].setncatts(file_mse[name].__dict__)
    file_mses[name][:] = file_mse[name][:]
    
mses = file_mses.createVariable('mses', 'f4', ("time","lat","lon"))
mses.units = "J/kg"
mses.long_name = "moist static energy"
mses[:,:,:] = ms

file_mses.close()
