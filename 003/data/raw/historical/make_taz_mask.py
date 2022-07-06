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
path_tas = sys.argv[4]

# open files
file_ta = Dataset(path_ta, 'r')
file_tas = Dataset(path_tas, 'r')
file_ps = Dataset(path_ps, 'r')

# read data
ps = np.squeeze(file_ps.variables['ps'][:]) # (mon x lat x lon)
ta = file_ta.variables['ta'][:] # (mon x lev x lat x lon)
tas = file_tas.variables['tas'][:] # (mon x lev x lat x lon)
plev = file_ta.variables['plev'][:]

# for datasets that fill data below surface as missing data, fill with nans
ta = ta.filled(fill_value=np.nan)
plev = plev.filled(fill_value=np.nan)

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
ps_tile = np.transpose( np.tile(ps, [ta.shape[1], 1, 1, 1]), [1,0,2,3])
pa = np.transpose(np.tile(plev, [ta.shape[0], ta.shape[2], ta.shape[3], 1]), [0,3,1,2])
subsurf = pa > ps_tile
ta[subsurf] = np.nan

# interpolate data after inserting surface value
tafill = np.empty_like(ta)
for itime in tqdm(np.arange(ta.shape[0])):
    for ilat in np.arange(ta.shape[2]):
        for ilon in np.arange(ta.shape[3]):
            # remove subsurface data
            qs0 = tas[itime,ilat,ilon]
            q0 = ta[itime,:,ilat,ilon]
            ps0 = ps[itime,ilat,ilon]
            p0 = pa[itime,:,ilat,ilon]
            q0 = q0[ p0 < ps0 ]
            p0 = p0[ p0 < ps0 ]
            # insert surface data
            q0 = np.insert(q0, 0, qs0)
            p0 = np.insert(p0, 0, ps0)
            # reinterpolate to original grid
            fint = interpolate.interp1d(p0, q0, bounds_error=False)
            tafill[itime,:,ilat,ilon] = fint(plev)

ta = tafill

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
taz.long_name = "temperature"
taz[:] = mz

file_taz.close()
