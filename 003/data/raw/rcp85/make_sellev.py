import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to mse
path_in = sys.argv[1]
path_out = sys.argv[2]
path_ps = sys.argv[3] 
varname = sys.argv[4] 
slev = float(sys.argv[5]) # selected level

# open files
file_in = Dataset(path_in, 'r')
file_ps = Dataset(path_ps, 'r')

# read data
vin = np.squeeze(file_in.variables[varname][:]) # (mon x lev x lat x lon)
plev = np.squeeze(file_in.variables['plev'][:]) # (lev)
lat3d = np.squeeze(file_in.variables['lat'][:]) # (lat)
lat2d = np.squeeze(file_ps.variables['lat'][:]) # (lat)
ps = np.squeeze(file_ps.variables['ps'][:]) # (mon x lat x lon)

# for datasets that fill data below surface as missing data, fill with nans
vin = vin.filled(fill_value=np.nan)

# Check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
if not np.array_equal(lat2d,lat3d):
    print('\nInterpolating ps to 3d lat grid...\n')
    filledps = ps.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filledps, axis=1)
    ps = f(filledlat3d)
    filledps = None; f = None;
    print('\nDone.\n')

# Mask out data below surface
ps_tile = np.transpose( np.tile(ps, [vin.shape[1], 1, 1, 1]), [1,0,2,3] )
pa = np.transpose(np.tile(plev, [vin.shape[0], vin.shape[2], vin.shape[3], 1]), [0,3,1,2])
subsurf = pa > ps_tile
vin[subsurf] = np.nan
pa = None; ps_tile = None

# select level
fint = interpolate.interp1d(plev, vin, axis=1, bounds_error=False)
svin = fint(slev)

# save file as netCDF
file_out = Dataset(path_out, "w", format='NETCDF4_CLASSIC')

# copy dimensions 
for name, dimension in file_in.dimensions.items():
    file_out.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables
for name, variable in file_in.variables.items():
    if any(name in s for s in [varname]):
        continue
    
    x = file_out.createVariable(name, variable.datatype, variable.dimensions)
    file_out[name].setncatts(file_in[name].__dict__)
    file_out[name][:] = file_in[name][:]
    
outa = file_out.createVariable(varname, 'f4', ("time", "lat", "lon"))
outa.units = file_in[varname].units
outa.long_name = file_in[varname].long_name
outa[:] = svin

file_out.close()
