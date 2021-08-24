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
path_mse = sys.argv[1]
path_mses = sys.argv[2]

# open files
file_mse = Dataset(path_mse, 'r')

# read data
mse = file_mse.variables['mse'][:] # (day x lev x lat x lon)
lev = file_mse.variables['lev'][:] # (lev)

# ms = np.empty([mse.shape[0],mse.shape[2],mse.shape[3]])
# extract surface MSE
if lev[1]-lev[0]>0:
    ms = mse[:,-1,:,:]
else:
    ms = mse[:,0,:,:]

# save file as netCDF
file_mses = Dataset(path_mses, "w", format='NETCDF4_CLASSIC')

# copy dimensions from mse file
for name, dimension in file_mse.dimensions.items():
    if any(name in s for s in ['lev', 'hyai', 'hybi', 'hyam', 'hybm']):
        continue

    file_mses.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables from mse file
for name, variable in file_mse.variables.items():
    if any(name in s for s in ['mse', 'lev', 'hyai', 'hybi', 'hyam', 'hybm']):
        continue
    
    x = file_mses.createVariable(name, variable.datatype, variable.dimensions)
    file_mses[name].setncatts(file_mse[name].__dict__)
    file_mses[name][:] = file_mse[name][:]
    
mses = file_mses.createVariable('mses', 'f4', ("time","lat","lon"))
mses.units = "J/kg"
mses.long_name = "surface moist static energy"
mses[:,:,:] = ms

file_mse.close()
