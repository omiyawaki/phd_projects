import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;
p0 = 1100.e2 # bottom integral bound [hPa]
pt = 0. # top integral bound [hPa]

# paths to ps and ta files
path_ps = sys.argv[1]
path_ta = sys.argv[2]
path_tmax = sys.argv[3]

# open files
file_ps = Dataset(path_ps, 'r')
file_ta = Dataset(path_ta, 'r')

# read data
ps = np.squeeze(file_ps.variables['aps'][:]) # (lat x lon)
ta = file_ta.variables['t'][:] # (mon x lev x lat x lon)
try:
    plev = file_ta.variables['plev'][:]
except:
    plev = file_ta.variables['lev'][:]

# for datasets that fill data below surface as missing data, fill with nans
ta = ta.filled(fill_value=np.nan)

# check if 2d and 3d lat grids are the same and interpolate 3d data to 2d grid if different
lat2d = file_ps.variables['lat'][:] 
lat3d = file_ta.variables['lat'][:] 
if not np.array_equal(lat2d,lat3d):
    # filledps = ps.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat3d, ta, axis=2, bounds_error=False)
    ta = f(filledlat2d)
    f = None;

# replace subsurface data in general with nans (important for computing dhdt near the surface)
ps_tile = np.transpose(np.tile(ps, [ta.shape[1], 1, 1, 1]), [1,0,2,3])
pa = np.transpose(np.tile(plev, [ta.shape[0], ta.shape[2], ta.shape[3], 1]), [0,3,1,2])
subsurf = pa > ps_tile
ta[subsurf] = np.nan

# mask out data above 500 hPa (focus on lower troposphere only)
upper = pa < 500e2 
ta[upper] = np.nan

# find local maximum temperature (in height)
tm = np.nanmax(ta, axis=1)

# save file as netCDF
file_tmax = Dataset(path_tmax, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_tmax.setncatts(file_ps.__dict__)

# copy dimensions from ps file
for name, dimension in file_ps.dimensions.items():
    # if any(name in s for s in ['plev']):
    #     continue
    file_tmax.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_ps.variables.items():
    if any(name in s for s in ['aps']):
        continue
    
    x = file_tmax.createVariable(name, variable.datatype, variable.dimensions)
    file_tmax[name].setncatts(file_ps[name].__dict__)
    file_tmax[name][:] = file_ps[name][:]
    
tmax = file_tmax.createVariable('tmax', 'f4', ("time","lat","lon"))
tmax.units = "K"
tmax.long_name = "local maximum temperature in the lower troposphere (below 500 hPa)"
tmax[:,:,:] = tm

file_tmax.close()
