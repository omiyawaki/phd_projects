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
path_dse = sys.argv[3]

# open files
file_ta = Dataset(path_ta, 'r')
file_zg = Dataset(path_zg, 'r')

# read data
ta = file_ta.variables['ta'][:] # (day x lev x lat x lon)
zg = file_zg.variables['zg'][:] # (day x lev x lat x lon)

# for datasets that fill data below surface as missing data, fill with nans
ta = ta.filled(fill_value=np.nan)
zg = zg.filled(fill_value=np.nan)

# calculate dse
m = cpd*ta + g*zg
ta = None; zg = None; # collect garbage

# save file as netCDF
file_dse = Dataset(path_dse, "w", format='NETCDF4_CLASSIC')

# # copy attributes from zg file
# file_dse.setncatts(file_zg.__dict__)

# copy dimensions from zg file
for name, dimension in file_zg.dimensions.items():
    file_dse.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables from zg file
for name, variable in file_zg.variables.items():
    if any(name in s for s in ['zg']):
        continue
    
    x = file_dse.createVariable(name, variable.datatype, variable.dimensions)
    file_dse[name].setncatts(file_zg[name].__dict__)
    file_dse[name][:] = file_zg[name][:]
    
dse = file_dse.createVariable('dse', 'f4', ("time","plev","lat","lon"))
dse.units = "J/kg"
dse.long_name = "dry static energy"
dse[:,:,:,:] = m

file_dse.close()
