import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ta
path_ta = sys.argv[1]
path_ps = sys.argv[2]
path_taz = sys.argv[3]

# open files
file_ta = Dataset(path_ta, 'r')
file_ps = Dataset(path_ps, 'r')

# read data
ps = np.squeeze(file_ps.variables['ps'][:]) # (lat x lon)
ta = file_ta.variables['ta'][:] # (mon x lev x lat x lon)
plev = file_ta.variables['plev'][:]

# for datasets that fill data below surface as missing data, fill with nans
ta = ta.filled(fill_value=np.nan)

# check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
lat2d = file_ps.variables['lat'][:] 
lat3d = file_ta.variables['lat'][:] 
if not np.array_equal(lat2d,lat3d):
    filledps = ps.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filledps, axis=0)
    ps = f(filledlat3d)
    filledps = None; f = None;

# replace subsurface data in general with nans (important for computing dhdt near the surface)
ps_tile = np.tile(ps, [ta.shape[0], ta.shape[1], 1, 1])
pa = np.transpose(np.tile(plev, [ta.shape[0], ta.shape[2], ta.shape[3], 1]), [0,3,1,2])
subsurf = pa > ps_tile
ta[subsurf] = np.nan

# take zonal mean
# mz = np.transpose( np.tile(np.nanmean(ta,axis=2), [ta.shape[2],1,1]), [1,2,0] )
mz = np.nanmean(ta,axis=3)

# save file as netCDF
file_taz = Dataset(path_taz, "w", format='NETCDF4_CLASSIC')

# # copy attributes from zg file
# file_taz.setncatts(file_zg.__dict__)

# copy dimensions from ta file
for name, dimension in file_ta.dimensions.items():
    file_taz.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables
for name, variable in file_ta.variables.items():
    if any(name in s for s in ['ta']):
        continue
    
    x = file_taz.createVariable(name, variable.datatype, variable.dimensions)
    file_taz[name].setncatts(file_ta[name].__dict__)
    file_taz[name][:] = file_ta[name][:]
    
taz = file_taz.createVariable('ta', 'f4', ("time","plev","lat"))
taz.units = "K"
taz.long_name = "air temperature"
taz[:] = mz

file_taz.close()
