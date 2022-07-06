import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to hus
path_hus = sys.argv[1]
path_ps = sys.argv[2]
path_husz = sys.argv[3]
path_huss = sys.argv[4]

# open files
file_hus = Dataset(path_hus, 'r')
file_huss = Dataset(path_huss, 'r')
file_ps = Dataset(path_ps, 'r')

# read data
ps = np.squeeze(file_ps.variables['ps'][:]) # (mon x lat x lon)
hus = file_hus.variables['hus'][:] # (mon x lev x lat x lon)
huss = file_huss.variables['huss'][:] # (mon x lev x lat x lon)
plev = file_hus.variables['plev'][:]

# for datasets that fill data below surface as missing data, fill with nans
hus = hus.filled(fill_value=np.nan)
plev = plev.filled(fill_value=np.nan)

# check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
lat2d = file_ps.variables['lat'][:] 
lat3d = file_hus.variables['lat'][:] 
if not np.array_equal(lat2d,lat3d):
    filledps = ps.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filledps, axis=0)
    ps = f(filledlat3d)
    filledps = None; f = None;

# replace subsurface data in general with nans (important for computing dhdt near the surface)
ps_tile = np.transpose( np.tile(ps, [hus.shape[1], 1, 1, 1]), [1,0,2,3])
pa = np.transpose(np.tile(plev, [hus.shape[0], hus.shape[2], hus.shape[3], 1]), [0,3,1,2])
subsurf = pa > ps_tile
hus[subsurf] = np.nan

# interpolate data after inserting surface value
husfill = np.empty_like(hus)
for itime in tqdm(np.arange(hus.shape[0])):
    for ilat in np.arange(hus.shape[2]):
        for ilon in np.arange(hus.shape[3]):
            # remove subsurface data
            qs0 = huss[itime,ilat,ilon]
            q0 = hus[itime,:,ilat,ilon]
            ps0 = ps[itime,ilat,ilon]
            p0 = pa[itime,:,ilat,ilon]
            q0 = q0[ p0 < ps0 ]
            p0 = p0[ p0 < ps0 ]
            # insert surface data
            q0 = np.insert(q0, 0, qs0)
            p0 = np.insert(p0, 0, ps0)
            # reinterpolate to original grid
            fint = interpolate.interp1d(p0, q0, bounds_error=False)
            husfill[itime,:,ilat,ilon] = fint(plev)

hus = husfill

# take zonal mean
# mz = np.transpose( np.tile(np.nanmean(hus,axis=2), [hus.shape[2],1,1]), [1,2,0] )
mz = np.nanmean(hus,axis=3)

# save file as netCDF
file_husz = Dataset(path_husz, "w", format='NETCDF4_CLASSIC')

# # copy attributes from zg file
# file_husz.setncatts(file_zg.__dict__)

# copy dimensions from hus file
for name, dimension in file_hus.dimensions.items():
    file_husz.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables
for name, variable in file_hus.variables.items():
    if any(name in s for s in ['hus']):
        continue
    
    x = file_husz.createVariable(name, variable.datatype, variable.dimensions)
    file_husz[name].setncatts(file_hus[name].__dict__)
    file_husz[name][:] = file_hus[name][:]
    
husz = file_husz.createVariable('hus', 'f4', ("time","plev","lat"))
husz.units = "kg/kg"
husz.long_name = "specific humidity"
husz[:] = mz

file_husz.close()
