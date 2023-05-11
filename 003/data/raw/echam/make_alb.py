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
path_sradsu = sys.argv[1]
path_sradsd = sys.argv[2]
path_alb = sys.argv[3]

# open files
file_sradsu = Dataset(path_sradsu, 'r')
file_sradsd = Dataset(path_sradsd, 'r')

# read data
sradsu = file_sradsu.variables['sradsu'][:] # (mon x lat x lon)
sradsd = file_sradsd.variables['sradsd'][:] # (mon x lat x lon)

# compute albedo
albedo = -sradsu/sradsd
albedo[sradsd<=0.1]=np.nan

# save file as netCDF
file_alb = Dataset(path_alb, "w", format='NETCDF4_CLASSIC')

# copy attributes from sradsu file
file_alb.setncatts(file_sradsu.__dict__)

# copy dimensions from sradsu file
for name, dimension in file_sradsu.dimensions.items():
    file_alb.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from sradsu file
for name, variable in file_sradsu.variables.items():
    if any(name in s for s in ['sradsu']):
        continue
    
    x = file_alb.createVariable(name, variable.datatype, variable.dimensions)
    file_alb[name].setncatts(file_sradsu[name].__dict__)
    file_alb[name][:] = file_sradsu[name][:]
    
alb = file_alb.createVariable('alb', 'f4', ("time","lat","lon"))
alb.units = "unitless"
alb.long_name = "surface albedo (diagnosed)"
alb[:] = albedo

file_alb.close()
