import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to t
path_t = sys.argv[1]
path_ps = sys.argv[2]
path_tz = sys.argv[3]

# open files
file_t = Dataset(path_t, 'r')
file_ps = Dataset(path_ps, 'r')

# read data
ps = np.squeeze(file_ps.variables['aps'][:]) # (lat x lon)
t = file_t.variables['t'][:] # (mon x plev x lat x lon)
try:
    plev = file_t.variables['plev'][:]
except:
    plev = file_t.variables['lev'][:]

# for datsets that fill data below surface as missing data, fill with nans
t = t.filled(fill_value=np.nan)

# check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
lat2d = file_ps.variables['lat'][:] 
lat3d = file_t.variables['lat'][:] 
if not np.array_equal(lat2d,lat3d):
    filledps = ps.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filledps, axis=0)
    ps = f(filledlat3d)
    filledps = None; f = None;

# replace subsurface data in general with nans (importnt for computing dhdt near the surface)
ps_tile = np.tile(ps, [t.shape[0], t.shape[1], 1, 1])
pa = np.transpose(np.tile(plev, [t.shape[0], t.shape[2], t.shape[3], 1]), [0,3,1,2])
subsurf = pa > ps_tile
t[subsurf] = np.nan

# tke zonal mean
# mz = np.transpose( np.tile(np.nanmean(t,axis=2), [t.shape[2],1,1]), [1,2,0] )
mz = np.nanmean(t,axis=3)

# save file as netCDF
file_tz = Dataset(path_tz, "w", format='NETCDF4_CLASSIC')

# # copy attributes from zg file
# file_tz.setncatts(file_zg.__dict__)

# copy dimensions from t file
for name, dimension in file_t.dimensions.items():
    file_tz.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables
for name, variable in file_t.variables.items():
    if any(name in s for s in ['t']):
        continue
    
    x = file_tz.createVariable(name, variable.datatype, variable.dimensions)
    file_tz[name].setncatts(file_t[name].__dict__)
    file_tz[name][:] = file_t[name][:]
    
try:
    tz = file_tz.createVariable('t', 'f4', ("time","plev","lat"))
except:
    tz = file_tz.createVariable('t', 'f4', ("time","lev","lat"))
tz.units = "K"
tz.long_name = "air temperature"
tz[:] = mz

file_tz.close()
