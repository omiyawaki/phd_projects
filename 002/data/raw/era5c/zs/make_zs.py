import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for ERA5

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps, ta, zg, and hus files
path_ps = ../srfc/era5c_srfc_1980_1980.nc
path_zg = ../zg/era5c_zg_1980_1980.nc
path_zs = ./era5c_zs_1980_1980.nc

# open files
file_ps = Dataset(path_ps, 'r')
file_zg = Dataset(path_zg, 'r')

# read data
ps = file_ps.variables['sp'][:] # (mon x lat x lon)
zg = file_zg.variables['z'][:] # (day x lev x lat x lon)
# zg = zg/g # era5 geopot height is gz
plev = 100 * file_zg.variables['level'][:]  # convert hPa to Pa

# for datasets that fill data below surface as missing data, fill with nans
zg = zg.filled(fill_value=np.nan)

# check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
lat2d = file_ps.variables['lat'][:] # (mon x lat x lon)
lat3d = file_zg.variables['lat'][:] # (day x lev x lat x lon)
lat3d_zg = file_zg.variables['lat'][:] # (day x lev x lat x lon)
if not np.array_equal(lat2d,lat3d_zg):
    filledps = ps.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filledps, axis=1)
    ps = f(filledlat3d)
    # if not np.array_equal(lat3d,lat3d_zg):
    #     f = interpolate.interp1d(filledlat2d, ta, axis=2)
    #     ta = f(filledlat3d)
    #     f = interpolate.interp1d(filledlat2d, hus, axis=2)
    #     hus = f(filledlat3d)
    filledps = None; f = None;

# evaluate zs by interpolating zg at ps

# replace subsurface data in general with nans (important for computing dhdt near the surface)
ps3d = np.tile(ps, [plev.size, 1, 1, 1])
ps3d = np.transpose( ps3d, [1, 0, 2, 3] )
pa3d = np.tile(plev, [ps.shape[0], ps.shape[1], ps.shape[2], 1])
pa3d = np.transpose( pa3d, [0, 3, 1, 2] )
idx_below = pa3d > ps3d
ps3d = None;

pa3d[idx_below]=np.nan
zg[idx_below]=np.nan


# save file as netCDF
file_tend = Dataset(path_tend, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_tend.setncatts(file_ps.__dict__)

# copy dimensions from zg file
for name, dimension in file_zg.dimensions.items():
    if any(name in s for s in ['level']):
        continue
    file_tend.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_zg.variables.items():
    if any(name in s for s in ['z' 'level']):
        continue
    
    x = file_tend.createVariable(name, variable.datatype, variable.dimensions)
    file_tend[name].setncatts(file_zg[name].__dict__)
    file_tend[name][:] = file_zg[name][:]
    
tend = file_tend.createVariable('tend', 'f4', ("time","lat","lon"))
tend.units = "W/m^2"
tend.long_name = "time tendency of vertically integrated moist static energy"
tend[:,:,:] = dvmsedt

file_tend.close()
