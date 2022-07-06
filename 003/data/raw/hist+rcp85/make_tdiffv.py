import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts/misc')
from tvregdiff import TVRegDiff
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP
cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps and mse files
path_aht = sys.argv[1]
path_gmse = sys.argv[2]
path_diffv = sys.argv[3]

# open files
file_aht = Dataset(path_aht, 'r')
file_gmse = Dataset(path_gmse, 'r')

# read data
aht = np.squeeze(file_aht.variables['aht'][:]) # (mon x lat)
gmse = np.squeeze(file_gmse.variables['gmse'][:]) # (mon x lat)

# compute diffusivity
if len(aht.shape)==1:
    diffv = np.empty([1,len(gmse)])
    diffv[0,:] = - aht / gmse
else:
    diffv = - aht / (2*np.pi*gmse)

# save file as netCDF
file_diffv = Dataset(path_diffv, "w", format='NETCDF4_CLASSIC')

# copy attributes
file_diffv.setncatts(file_aht.__dict__)

# copy dimensions
for name, dimension in file_aht.dimensions.items():
    file_diffv.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables
for name, variable in file_aht.variables.items():
    if any(name in s for s in ['aht']):
        continue
    
    x = file_diffv.createVariable(name, variable.datatype, variable.dimensions)
    file_diffv[name].setncatts(file_aht[name].__dict__)
    file_diffv[name][:] = file_aht[name][:]
    
var_diffv = file_diffv.createVariable('tdiffv', 'f4', ("time","lat"))
var_diffv.units = "kg s**-1"
var_diffv.long_name = "92500 Pa Atmospheric MSE diffusivity"
var_diffv[:] = diffv
file_diffv.close()
