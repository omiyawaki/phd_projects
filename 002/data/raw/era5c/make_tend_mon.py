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
path_ps = sys.argv[1]
path_ta = sys.argv[2]
path_zg = sys.argv[3]
path_hus = sys.argv[4]
path_tend = sys.argv[5]
path_tas = sys.argv[6]
path_orog = sys.argv[7]
path_huss = sys.argv[8]

# open files
file_ps = Dataset(path_ps, 'r')
file_ta = Dataset(path_ta, 'r')
file_zg = Dataset(path_zg, 'r')
file_hus = Dataset(path_hus, 'r')
# file_tas = Dataset(path_tas, 'r')
# file_orog = Dataset(path_orog, 'r')

# read data
# orog = file_orog.variables['orog'][:] # (lat x lon)
ps = file_ps.variables['sp'][:] # (mon x lat x lon)
# tas = file_tas.variables['tas'][:] # (mon x lev x lat x lon)
ta = file_ta.variables['t'][:] # (day x lev x lat x lon)
zg = file_zg.variables['z'][:] # (day x lev x lat x lon)
zg = zg/g # era5 geopot height is gz
hus = file_hus.variables['q'][:] # (day x lev x lat x lon)
plev = 100 * file_ta.variables['level'][:]  # convert hPa to Pa

# for datasets that fill data below surface as missing data, fill with nans
ta = ta.filled(fill_value=np.nan)
zg = zg.filled(fill_value=np.nan)
hus = hus.filled(fill_value=np.nan)

# check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
lat2d = file_ps.variables['lat'][:] # (mon x lat x lon)
lat3d = file_zg.variables['lat'][:] # (day x lev x lat x lon)
lat3d_zg = file_zg.variables['lat'][:] # (day x lev x lat x lon)
if not np.array_equal(lat2d,lat3d_zg):
    filledps = ps.filled(fill_value=np.nan)
    # filledtas = tas.filled(fill_value=np.nan)
    # filledorog = orog.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filledps, axis=1)
    ps = f(filledlat3d)
    if not np.array_equal(lat3d,lat3d_zg):
        f = interpolate.interp1d(filledlat2d, ta, axis=2)
        ta = f(filledlat3d)
        f = interpolate.interp1d(filledlat2d, hus, axis=2)
        hus = f(filledlat3d)
    # f = interpolate.interp1d(filledlat2d, filledtas, axis=1)
    # tas = f(filledlat3d)
    # f = interpolate.interp1d(filledlat2d, filledorog, axis=0)
    # orog = f(filledlat3d)
    filledps = None; filledtas = None; filledorog = None; f = None;

# replace subsurface data in general with nans (important for computing dhdt near the surface)
ps3d = np.tile(ps, [plev.size, 1, 1, 1])
ps3d = np.transpose( ps3d, [1, 0, 2, 3] )
pa3d = np.tile(plev, [ps.shape[0], ps.shape[1], ps.shape[2], 1])
pa3d = np.transpose( pa3d, [0, 3, 1, 2] )
idx_below = pa3d > ps3d
ps3d = None;

pa3d[idx_below]=np.nan
ta[idx_below]=np.nan
zg[idx_below]=np.nan
hus[idx_below]=np.nan

# calculate MSE
mse = cpd*ta + g*zg + L*hus
ta = None; zg = None; # collect garbage

# take time tendency
# 3d mse tendency
dmsedt = np.empty(mse.shape) 
dmsedt[1:-1,:,:,:] = (mse[2:,:,:,:]-mse[0:-2,:,:,:])/(2*365*86400/12)
dmsedt[0] = (mse[1,:,:,:]-mse[0,:,:,:])/(365*86400/12)
dmsedt[-1] = (mse[-1,:,:,:]-mse[-2,:,:,:])/(365*86400/12)

# compute vertical integral
if plev[1]-plev[0]>0: # if pressure increases with index
    dvmsedt = 1/g*np.trapz(dmsedt, pa3d, axis=1)
else:
    dvmsedt = -1/g*np.trapz(dmsedt, pa3d, axis=1)

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
