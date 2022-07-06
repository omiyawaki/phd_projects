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
path_t = sys.argv[1]
path_q = sys.argv[2]
path_z = sys.argv[3]
path_mse = sys.argv[4]

# open files
file_t = Dataset(path_t, 'r')
file_q = Dataset(path_q, 'r')
file_z = Dataset(path_z, 'r')

# read data
ta = file_t.variables['t'][:] # (mon x lev x lat x lon)
hus = file_q.variables['q'][:] # (mon x lev x lat x lon)
zg = file_z.variables['geopoth'][:] # (mon x lev x lat x lon)

# for datasets that fill data below surface as missing data, fill with nans
ta = ta.filled(fill_value=np.nan)
zg = zg.filled(fill_value=np.nan)
hus = hus.filled(fill_value=np.nan)

# calculate MSE
m = cpd*ta + g*zg + L*hus
ta = None; zg = None; hus = None; # collect garbage

# save file as netCDF
file_mse = Dataset(path_mse, "w", format='NETCDF4_CLASSIC')

# # copy attributes from atm file
# file_mse.setncatts(file_atm.__dict__)

# copy dimensions from atm file
for name, dimension in file_t.dimensions.items():
    file_mse.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables from atm file
for name, variable in file_t.variables.items():
    if any(name in s for s in ['t']):
        continue
    
    x = file_mse.createVariable(name, variable.datatype, variable.dimensions)
    file_mse[name].setncatts(file_t[name].__dict__)
    file_mse[name][:] = file_t[name][:]
    
try:
    mse = file_mse.createVariable('mse', 'f4', ("time","plev","lat","lon"))
except:
    mse = file_mse.createVariable('mse', 'f4', ("time","lev","lat","lon"))
mse.units = "J/kg"
mse.long_name = "moist static energy"
mse[:,:,:,:] = m

file_mse.close()
