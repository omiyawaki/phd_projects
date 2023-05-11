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
path_sradsd = sys.argv[1]
path_srad0d = sys.argv[2]
path_swtau = sys.argv[3]

# open files
file_sradsd = Dataset(path_sradsd, 'r')
file_srad0d = Dataset(path_srad0d, 'r')

# read data
sradsd = file_sradsd.variables['sradsd'][:] # (mon x lat x lon)
srad0d = file_srad0d.variables['srad0d'][:] # (mon x lat x lon)

# compute atmospheric sw optical thickness
tau = -np.log(sradsd/srad0d)
tau[sradsd<=0.1]=np.nan

# save file as netCDF
file_swtau = Dataset(path_swtau, "w", format='NETCDF4_CLASSIC')

# copy attributes from sradsd file
file_swtau.setncatts(file_sradsd.__dict__)

# copy dimensions from sradsd file
for name, dimension in file_sradsd.dimensions.items():
    file_swtau.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from sradsd file
for name, variable in file_sradsd.variables.items():
    if any(name in s for s in ['sradsd']):
        continue
    
    x = file_swtau.createVariable(name, variable.datatype, variable.dimensions)
    file_swtau[name].setncatts(file_sradsd[name].__dict__)
    file_swtau[name][:] = file_sradsd[name][:]
    
swtau = file_swtau.createVariable('swtau', 'f4', ("time","lat","lon"))
swtau.units = "unitless"
swtau.long_name = "atmospheric sw optical thickness"
swtau[:] = tau

file_swtau.close()
