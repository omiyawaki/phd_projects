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

# paths to ps and mse files
path_ps = sys.argv[1]
path_mse = sys.argv[2]
path_beta = sys.argv[3]
path_tendv = sys.argv[4]

# open files
file_ps = Dataset(path_ps, 'r')
file_mse = Dataset(path_mse, 'r')
file_beta = Dataset(path_beta, 'r')

# read data
ps = np.squeeze(file_ps.variables['aps'][:]) # (lat x lon)
mse = file_mse.variables['mse'][:] # (mon x lev x lat x lon)
beta = np.squeeze(file_beta.variables['beta'][:]) # (lev x lat x lon)

try:
    plev = file_mse.variables['plev'][:]
except:
    plev = file_mse.variables['lev'][:]
plev_half = 1/2 * (plev[1:] + plev[:-1])
plev_half = np.sort(np.append(plev_half, [pt, p0]))
plev_full = np.sort(np.concatenate((plev, plev_half)))
if plev[1]-plev[0]<0:
    plev_half = plev_half[::-1]
    plev_full = plev_full[::-1]
dplev = plev_half[1:] - plev_half[:-1]

# for datasets that fill data below surface as missing data, fill with nans
mse = mse.filled(fill_value=np.nan)

# check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
lat2d = file_ps.variables['lat'][:] 
lat3d = file_mse.variables['lat'][:] 
if not np.array_equal(lat2d,lat3d):
    filledps = ps.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filledps, axis=0)
    ps = f(filledlat3d)
    filledps = None; f = None;

# replace subsurface data in general with nans (important for computing dhdt near the surface)
ps_tile = np.tile(ps, [mse.shape[0], mse.shape[1], 1, 1])
pa = np.transpose(np.tile(plev, [mse.shape[0], mse.shape[2], mse.shape[3], 1]), [0,3,1,2])
subsurf = pa > ps_tile
mse[subsurf] = np.nan

# take time tendency
# 3d mse tendency
dmsedt = np.empty(mse.shape) 
dmsedt[1:-1,:,:,:] = (mse[2:,:,:,:]-mse[0:-2,:,:,:])/(2*365*86400/12)
dmsedt[0] = (mse[1,:,:,:]-mse[0,:,:,:])/(365*86400/12)
dmsedt[-1] = (mse[-1,:,:,:]-mse[-2,:,:,:])/(365*86400/12)

# save file as netCDF
file_tendv = Dataset(path_tendv, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_tendv.setncatts(file_ps.__dict__)

# copy dimensions from mse file
for name, dimension in file_mse.dimensions.items():
    file_tendv.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_mse.variables.items():
    if any(name in s for s in ['mse']):
        continue
    
    x = file_tendv.createVariable(name, variable.datatype, variable.dimensions)
    file_tendv[name].setncatts(file_mse[name].__dict__)
    file_tendv[name][:] = file_mse[name][:]
    
try:
    tendv = file_tendv.createVariable('tendv', 'f4', ("time","plev","lat","lon"))
except:
    tendv = file_tendv.createVariable('tendv', 'f4', ("time","lev","lat","lon"))
tendv.units = "W/kg"
tendv.long_name = "time tendency of moist static energy"
tendv[:,:,:] = dmsedt

file_tendv.close()
