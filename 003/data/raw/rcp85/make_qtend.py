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

# paths to ps and hus files
path_ps = sys.argv[1]
path_hus = sys.argv[2]
path_beta = sys.argv[3]
path_qtend = sys.argv[4]

# open files
file_ps = Dataset(path_ps, 'r')
file_hus = Dataset(path_hus, 'r')
file_beta = Dataset(path_beta, 'r')

# read data
ps = np.squeeze(file_ps.variables['ps'][:]) # (lat x lon)
hus = file_hus.variables['hus'][:] # (mon x lev x lat x lon)
beta = np.squeeze(file_beta.variables['beta'][:]) # (lev x lat x lon)

plev = file_hus.variables['plev'][:]
plev_half = 1/2 * (plev[1:] + plev[:-1])
plev_half = np.sort(np.append(plev_half, [pt, p0]))
plev_full = np.sort(np.concatenate((plev, plev_half)))
if plev[1]-plev[0]<0:
    plev_half = plev_half[::-1]
    plev_full = plev_full[::-1]
dplev = plev_half[1:] - plev_half[:-1]

# for datasets that fill data below surface as missing data, fill with nans
hus = hus.filled(fill_value=np.nan)

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
ps_tile = np.tile(ps, [hus.shape[0], hus.shape[1], 1, 1])
pa = np.transpose(np.tile(plev, [hus.shape[0], hus.shape[2], hus.shape[3], 1]), [0,3,1,2])
subsurf = pa > ps_tile
hus[subsurf] = np.nan

# take time qtendency
# 3d hus qtendency
dhusdt = np.empty(hus.shape) 
dhusdt[1:-1,:,:,:] = (hus[2:,:,:,:]-hus[0:-2,:,:,:])/(2*365*86400/12)
dhusdt[0] = (hus[1,:,:,:]-hus[0,:,:,:])/(365*86400/12)
dhusdt[-1] = (hus[-1,:,:,:]-hus[-2,:,:,:])/(365*86400/12)

# compute vertical integral
if plev[1]-plev[0]>0: # if pressure increases with index
    # dvhusdt = 1/g*np.trapz(dhusdt, pa3d, axis=1)
    dvhusdt = 1/g * np.nansum( beta[np.newaxis,...] * dhusdt * dplev[np.newaxis,...,np.newaxis,np.newaxis], axis=1)
else:
    # dvhusdt = -1/g*np.trapz(dhusdt, pa3d, axis=1)
    dvhusdt = -1/g * np.nansum( beta[np.newaxis,...] * dhusdt * dplev[np.newaxis,...,np.newaxis,np.newaxis], axis=1)

# save file as netCDF
file_qtend = Dataset(path_qtend, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_qtend.setncatts(file_ps.__dict__)

# copy dimensions from hus file
for name, dimension in file_hus.dimensions.items():
    if any(name in s for s in ['plev']):
        continue
    file_qtend.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_hus.variables.items():
    if any(name in s for s in ['hus' 'plev', 'plev_bnds']):
        continue
    
    x = file_qtend.createVariable(name, variable.datatype, variable.dimensions)
    file_qtend[name].setncatts(file_hus[name].__dict__)
    file_qtend[name][:] = file_hus[name][:]
    
qtend = file_qtend.createVariable('qtend', 'f4', ("time","lat","lon"))
qtend.units = "kg m**-2 s**-1"
qtend.long_name = "time qtendency of vertically integrated specific humidity"
qtend[:,:,:] = dvhusdt

file_qtend.close()
