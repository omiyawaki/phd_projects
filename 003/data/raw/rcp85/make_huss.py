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
path_hus = sys.argv[2]
path_huss = sys.argv[3]

# open files
file_ps = Dataset(path_ps, 'r')
file_hus = Dataset(path_hus, 'r')

# read data
ps = file_ps.variables['ps'][:] # (day x lev x lat x lon)
hus = file_hus.variables['hus'][:] # (day x lev x lat x lon)
plev = file_hus.variables['plev'][:]

# for datasets that fill data below surface as missing data, fill with nans
ps = ps.filled(fill_value=np.nan)
hus = hus.filled(fill_value=np.nan)

# Zonal and annual means of surface pressure
ps_z = np.mean(ps, axis=2)
ps_a = np.mean(ps, axis=(0))
ps_za = np.mean(ps, axis=(0,2))

ps_a_tile = np.tile(ps_a, [ps.shape[0], 1, 1])

# mask out data below surface
pa = np.transpose(np.tile(plev, [ps.shape[0], ps.shape[1], ps.shape[2], 1]), [0,3,1,2])
subsurf = pa > np.transpose(np.tile(ps_a_tile, [len(plev),1,1,1]),[1,0,2,3])
hus[subsurf] = np.nan
pa[subsurf] = np.nan

# interpolate hus at ps
ms = np.empty([hus.shape[0], hus.shape[2], hus.shape[3]])
for itime in tqdm(range(hus.shape[0])):
    for ilat in range(hus.shape[2]):
        for ilon in range(hus.shape[3]):
            plev_local = pa[itime,:,ilat,ilon]
            hus_local = hus[itime,:,ilat,ilon]

            plev_local = plev_local[~np.isnan(plev_local)]
            hus_local = hus_local[~np.isnan(hus_local)]

            f = interpolate.interp1d(plev_local, hus_local, fill_value='extrapolate')
            # f = interpolate.interp1d(plev, hus[itime,:,ilat,ilon], fill_value='extrapolate')
            ms[itime,ilat,ilon] = f(ps[itime,ilat,ilon])

# save file as netCDF
file_huss = Dataset(path_huss, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_huss.setncatts(file_ps.__dict__)

# copy dimensions from hus file
for name, dimension in file_hus.dimensions.items():
    if any(name in s for s in ['plev']):
        continue
    file_huss.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_hus.variables.items():
    if any(name in s for s in ['hus' 'plev']):
        continue
    
    x = file_huss.createVariable(name, variable.datatype, variable.dimensions)
    file_huss[name].setncatts(file_hus[name].__dict__)
    file_huss[name][:] = file_hus[name][:]
    
huss = file_huss.createVariable('huss', 'f4', ("time","lat","lon"))
huss.units = "kg/kg"
huss.long_name = "surface specific humidity"
huss[:,:,:] = ms

file_huss.close()
