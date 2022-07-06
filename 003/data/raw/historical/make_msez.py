import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to mse
path_mse = sys.argv[1]
path_msez = sys.argv[2]

# open files
file_mse = Dataset(path_mse, 'r')

# read data
mse = np.squeeze(file_mse.variables['mse'][:]) # (mon x lat x lon)

# for datasets that fill data below surface as missing data, fill with nans
mse = mse.filled(fill_value=np.nan)

# take zonal mean
mz = np.transpose( np.tile(np.nanmean(mse,axis=2), [mse.shape[2],1,1]), [1,2,0] )

# save file as netCDF
file_msez = Dataset(path_msez, "w", format='NETCDF4_CLASSIC')

# # copy attributes from zg file
# file_msez.setncatts(file_zg.__dict__)

# copy dimensions from mse file
for name, dimension in file_mse.dimensions.items():
    file_msez.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables
for name, variable in file_mse.variables.items():
    if any(name in s for s in ['mse']):
        continue
    
    x = file_msez.createVariable(name, variable.datatype, variable.dimensions)
    file_msez[name].setncatts(file_mse[name].__dict__)
    file_msez[name][:] = file_mse[name][:]
    
msez = file_msez.createVariable('mse', 'f4', ("time","lat","lon"))
msez.units = "J/kg"
msez.long_name = "moist static energy"
msez[:] = mz

file_msez.close()
