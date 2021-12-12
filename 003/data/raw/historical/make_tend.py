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
path_mse = sys.argv[2]
path_tend = sys.argv[3]

# open files
file_ps = Dataset(path_ps, 'r')
file_mse = Dataset(path_mse, 'r')

# read data
ps = file_ps.variables['ps'][:] # (mon x lat x lon)
mse = file_mse.variables['mse'][:] # (day x lev x lat x lon)
plev = file_mse.variables['plev'][:]

# for datasets that fill data below surface as missing data, fill with nans
mse = mse.filled(fill_value=np.nan)

# check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
lat2d = file_ps.variables['lat'][:] # (mon x lat x lon)
lat3d = file_mse.variables['lat'][:] # (day x lev x lat x lon)
if not np.array_equal(lat2d,lat3d):
    filledps = ps.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filledps, axis=1)
    ps = f(filledlat3d)
    filledps = None; f = None;

# replace subsurface data in general with nans (important for computing dhdt near the surface)
ps3d = np.tile(ps, [plev.size, 1, 1, 1])
ps3d = np.transpose( ps3d, [1, 0, 2, 3] )
pa3d = np.tile(plev, [ps.shape[0], ps.shape[1], ps.shape[2], 1])
pa3d = np.transpose( pa3d, [0, 3, 1, 2] )
idx_below = pa3d > ps3d
ps3d = None;

pa3d[idx_below]=np.nan
mse[idx_below]=np.nan

# take time tendency
# 3d mse tendency
dmsedt = np.empty(mse.shape) 
dmsedt[1:-1,:,:,:] = (mse[2:,:,:,:]-mse[0:-2,:,:,:])/(2*365*86400/12)
dmsedt[0] = (mse[1,:,:,:]-mse[0,:,:,:])/(365*86400/12)
dmsedt[-1] = (mse[-1,:,:,:]-mse[-2,:,:,:])/(365*86400/12)

# compute vertical integral
if plev[1]-plev[0]>0: # if pressure increases with index
    dvmsedt = 1/g*np.trapz(dmsedt, pa3d, axis=1)
else:
    dvmsedt = -1/g*np.trapz(dmsedt, pa3d, axis=1)

# save file as netCDF
file_tend = Dataset(path_tend, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_tend.setncatts(file_ps.__dict__)

# copy dimensions from mse file
for name, dimension in file_mse.dimensions.items():
    if any(name in s for s in ['plev']):
        continue
    file_tend.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_mse.variables.items():
    if any(name in s for s in ['mse' 'plev']):
        continue
    
    x = file_tend.createVariable(name, variable.datatype, variable.dimensions)
    file_tend[name].setncatts(file_mse[name].__dict__)
    file_tend[name][:] = file_mse[name][:]
    
tend = file_tend.createVariable('tend', 'f4', ("time","lat","lon"))
tend.units = "W/m^2"
tend.long_name = "time tendency of vertically integrated moist static energy"
tend[:,:,:] = dvmsedt

file_tend.close()
