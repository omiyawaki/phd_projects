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
path_ta = sys.argv[1]
path_zg = sys.argv[2]
path_hus = sys.argv[3]
path_mse = sys.argv[4]

# open files
file_ta = Dataset(path_ta, 'r')
file_zg = Dataset(path_zg, 'r')
file_hus = Dataset(path_hus, 'r')

# read data
ta = file_ta.variables['ts'][:] # (day x lat x lon)
zg = file_zg.variables['orog'][:] # (lat x lon)
zg = zg[np.newaxis,...]
hus = file_hus.variables['huss'][:] # (day x lat x lon)

# for datasets that fill data below surface as missing data, fill with nans
ta = ta.filled(fill_value=np.nan)
zg = zg.filled(fill_value=np.nan)
hus = hus.filled(fill_value=np.nan)

# calculate MSE
m = cpd*ta + g*zg + L*hus
ta = None; zg = None; hus = None; # collect garbage

# save file as netCDF
file_mse = Dataset(path_mse, "w", format='NETCDF4_CLASSIC')

# # copy attributes from ta file
# file_mse.setncatts(file_ta.__dict__)

# copy dimensions from ta file
for name, dimension in file_ta.dimensions.items():
    file_mse.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables from ta file
for name, variable in file_ta.variables.items():
    if any(name in s for s in ['ts']):
        continue
    
    x = file_mse.createVariable(name, variable.datatype, variable.dimensions)
    file_mse[name].setncatts(file_ta[name].__dict__)
    file_mse[name][:] = file_ta[name][:]
    
mse = file_mse.createVariable('mses', 'f4', ("time","lat","lon"))
mse.units = "J/kg"
mse.long_name = "surface moist static energy"
mse[:,:,:] = m

file_mse.close()
