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
path_va = sys.argv[2]
path_vas = sys.argv[3]

# open files
file_ps = Dataset(path_ps, 'r')
file_va = Dataset(path_va, 'r')

# read data
ps = file_ps.variables['aps'][:] # (day x lev x lat x lon)
va = file_va.variables['v'][:] # (day x lev x lat x lon)
plev = file_va.variables['lev'][:]

# for datasets that fill data below surface as missing data, fill with nans
ps = ps.filled(fill_value=np.nan)
va = va.filled(fill_value=np.nan)

# Zonal and annual means of surface pressure
ps_z = np.mean(ps, axis=2)
ps_a = np.mean(ps, axis=(0))
ps_za = np.mean(ps, axis=(0,2))

ps_a_tile = np.tile(ps_a, [ps.shape[0], 1, 1])

# mask out data below surface
pa = np.transpose(np.tile(plev, [ps.shape[0], ps.shape[1], ps.shape[2], 1]), [0,3,1,2])
subsurf = pa > np.transpose(np.tile(ps_a_tile, [len(plev),1,1,1]),[1,0,2,3])
va[subsurf] = np.nan
pa[subsurf] = np.nan

# interpolate va at ps
vs = np.empty([va.shape[0], va.shape[2], va.shape[3]])
for itime in tqdm(range(va.shape[0])):
    for ilat in range(va.shape[2]):
        for ilon in range(va.shape[3]):
            plev_local = pa[itime,:,ilat,ilon]
            va_local = va[itime,:,ilat,ilon]

            plev_local = plev_local[~np.isnan(plev_local)]
            va_local = va_local[~np.isnan(va_local)]

            f = interpolate.interp1d(plev_local, va_local, fill_value='extrapolate')
            # f = interpolate.interp1d(plev, va[itime,:,ilat,ilon], fill_value='extrapolate')
            vs[itime,ilat,ilon] = f(ps[itime,ilat,ilon])

# save file as netCDF
file_vas = Dataset(path_vas, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_vas.setncatts(file_ps.__dict__)

# copy dimensions from va file
for name, dimension in file_va.dimensions.items():
    if any(name in s for s in ['lev']):
        continue
    file_vas.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_va.variables.items():
    if any(name in s for s in ['u', 'v', 't', 'q', 'omega', 'geopoth', 'lev']):
        continue
    
    x = file_vas.createVariable(name, variable.datatype, variable.dimensions)
    file_vas[name].setncatts(file_va[name].__dict__)
    file_vas[name][:] = file_va[name][:]
    
vas = file_vas.createVariable('vas', 'f4', ("time","lat","lon"))
vas.units = "m/s"
vas.long_name = "surface northward velocity"
vas[:,:,:] = vs

file_vas.close()
