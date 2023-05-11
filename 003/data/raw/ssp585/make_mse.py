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
print(file_hus)
hus = file_hus.variables['hus'][:] # (day x lev x lat x lon)
ta = file_ta.variables['ta'][:] # (day x lev x lat x lon)
zg = file_zg.variables['zg'][:] # (day x lev x lat x lon)

# for datasets that fill data below surface as missing data, fill with nans
ta = ta.filled(fill_value=np.nan)
zg = zg.filled(fill_value=np.nan)
hus = hus.filled(fill_value=np.nan)
hus[np.abs(hus)>1e20]=np.nan

# calculate MSE
m = cpd*ta + g*zg + L*hus
ta = None; zg = None; hus = None; # collect garbage

# save file as netCDF
file_mse = Dataset(path_mse, "w", format='NETCDF4_CLASSIC')

# # copy attributes from zg file
# file_mse.setncatts(file_zg.__dict__)

# copy dimensions from zg file
for name, dimension in file_zg.dimensions.items():
    file_mse.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables from zg file
for name, variable in file_zg.variables.items():
    if any(name in s for s in ['zg']):
        continue
    
    x = file_mse.createVariable(name, variable.datatype, variable.dimensions)
    file_mse[name].setncatts(file_zg[name].__dict__)
    file_mse[name][:] = file_zg[name][:]
    
mse = file_mse.createVariable('mse', 'f4', ("time","plev","lat","lon"))
mse.units = "J/kg"
mse.long_name = "moist static energy"
mse[:,:,:,:] = m

file_mse.close()
