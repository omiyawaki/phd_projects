import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to hur
path_hur = sys.argv[1]
path_ps = sys.argv[2]
path_hurz = sys.argv[3]
path_hurs = sys.argv[4]

# open files
file_hur = Dataset(path_hur, 'r')
file_hurs = Dataset(path_hurs, 'r')
file_ps = Dataset(path_ps, 'r')

# read data
ps = np.squeeze(file_ps.variables['ps'][:]) # (mon x lat x lon)
hur = file_hur.variables['hur'][:] # (mon x lev x lat x lon)
hurs = file_hurs.variables['hurs'][:] # (mon x lev x lat x lon)
plev = file_hur.variables['plev'][:]

# for datasets that fill data below surface as missing data, fill with nans
hur = hur.filled(fill_value=np.nan)
plev = plev.filled(fill_value=np.nan)

# check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
lat2d = file_ps.variables['lat'][:] 
lat3d = file_hur.variables['lat'][:] 
if not np.array_equal(lat2d,lat3d):
    filledps = ps.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filledps, axis=0)
    ps = f(filledlat3d)
    filledps = None; f = None;

# replace subsurface data in general with nans (important for computing dhdt near the surface)
ps_tile = np.transpose( np.tile(ps, [hur.shape[1], 1, 1, 1]), [1,0,2,3])
pa = np.transpose(np.tile(plev, [hur.shape[0], hur.shape[2], hur.shape[3], 1]), [0,3,1,2])
subsurf = pa > ps_tile
hur[subsurf] = np.nan

# interpolate data after inserting surface value
hurfill = np.empty_like(hur)
for itime in tqdm(np.arange(hur.shape[0])):
    for ilat in np.arange(hur.shape[2]):
        for ilon in np.arange(hur.shape[3]):
            # remove subsurface data
            qs0 = hurs[itime,ilat,ilon]
            q0 = hur[itime,:,ilat,ilon]
            ps0 = ps[itime,ilat,ilon]
            p0 = pa[itime,:,ilat,ilon]
            q0 = q0[ p0 < ps0 ]
            p0 = p0[ p0 < ps0 ]
            # insert surface data
            q0 = np.insert(q0, 0, qs0)
            p0 = np.insert(p0, 0, ps0)
            # reinterpolate to original grid
            fint = interpolate.interp1d(p0, q0, bounds_error=False)
            hurfill[itime,:,ilat,ilon] = fint(plev)

hur = hurfill

# take zonal mean
# mz = np.transpose( np.tile(np.nanmean(hur,axis=2), [hur.shape[2],1,1]), [1,2,0] )
mz = np.nanmean(hur,axis=3)

# save file as netCDF
file_hurz = Dataset(path_hurz, "w", format='NETCDF4_CLASSIC')

# # copy attributes from zg file
# file_hurz.setncatts(file_zg.__dict__)

# copy dimensions from hur file
for name, dimension in file_hur.dimensions.items():
    file_hurz.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables
for name, variable in file_hur.variables.items():
    if any(name in s for s in ['hur']):
        continue
    
    x = file_hurz.createVariable(name, variable.datatype, variable.dimensions)
    file_hurz[name].setncatts(file_hur[name].__dict__)
    file_hurz[name][:] = file_hur[name][:]
    
hurz = file_hurz.createVariable('hur', 'f4', ("time","plev","lat"))
hurz.units = "%"
hurz.long_name = "relative humidity"
hurz[:] = mz

file_hurz.close()
