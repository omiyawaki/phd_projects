import sys
import numpy as np
from scipy import interpolate, integrate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps and mse files
path_srafs = sys.argv[1]
path_albedo = sys.argv[2]
path_srafsd = sys.argv[3]

# open files
file_srafs = Dataset(path_srafs, 'r')
file_albedo = Dataset(path_albedo, 'r')

# read data
srafs = file_srafs.variables['srafs'][:] # (mon x lat x lon)
albedo = file_albedo.variables['albedo'][:] # (mon x lat x lon)

# compute downwelling sw surface flux
sw_sfc_down = srafs/(1 - albedo)

# save file as netCDF
file_srafsd = Dataset(path_srafsd, "w", format='NETCDF4_CLASSIC')

# copy attributes from srafs file
file_srafsd.setncatts(file_srafs.__dict__)

# copy dimensions from srafs file
for name, dimension in file_srafs.dimensions.items():
    file_srafsd.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from srafs file
for name, variable in file_srafs.variables.items():
    if any(name in s for s in ['srafs']):
        continue
    
    x = file_srafsd.createVariable(name, variable.datatype, variable.dimensions)
    file_srafsd[name].setncatts(file_srafs[name].__dict__)
    file_srafsd[name][:] = file_srafs[name][:]
    
srafsd = file_srafsd.createVariable('srafsd', 'f4', ("time","lat","lon"))
srafsd.units = "W m**-2"
srafsd.long_name = "downwelling surface shortwave flux (clear sky)"
srafsd[:] = sw_sfc_down

file_srafsd.close()
