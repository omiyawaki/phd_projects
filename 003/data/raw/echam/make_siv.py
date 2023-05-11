import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps, si, zg, and sit files
path_si = sys.argv[1]
path_sit = sys.argv[2]
path_siv = sys.argv[3]

# open files
file_si = Dataset(path_si, 'r')
file_sit = Dataset(path_sit, 'r')

# read data
si = file_si.variables['seaice'][:] # (mon x lat x lon)
sit = file_sit.variables['siced'][:] # (mon x lat x lon)
lat = file_sit.variables['lat'][:] # (lat)
lon = file_sit.variables['lon'][:] # (lon)

# for Datasets that fill data below surface as missing data, fill with nans
si = si.filled(fill_value=np.nan)
sit = sit.filled(fill_value=np.nan)

# calculate grid cell area
r = 6371e3
dlat = -np.deg2rad(lat[1]-lat[0])
dlon = np.deg2rad(lon[1]-lon[0])
area = r**2 * np.cos(np.deg2rad(lat[:,None])) * dlat * dlon

# calculate siv
sicev = si * sit * area[None,:,:]

# save file as netCDF
file_siv = Dataset(path_siv, "w", format='NETCDF4_CLASSIC')

# # copy attributes from atm file
# file_siv.setncatts(file_atm.__dict__)

# copy dimensions from atm file
for name, dimension in file_si.dimensions.items():
    file_siv.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables from atm file
for name, variable in file_si.variables.items():
    if any(name in s for s in ['seaice']):
        continue
    
    x = file_siv.createVariable(name, variable.datatype, variable.dimensions)
    file_siv[name].setncatts(file_si[name].__dict__)
    file_siv[name][:] = file_si[name][:]
    
try:
    siv = file_siv.createVariable('siv', 'f4', ("time","lat","lon"))
except:
    siv = file_siv.createVariable('siv', 'f4', ("time","lat","lon"))
siv.units = "m**3"
siv.long_name = "sea ice volume"
siv[:,:,:] = sicev

file_siv.close()
