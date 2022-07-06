import sys
import numpy as np
from scipy import interpolate, integrate
from scipy.ndimage import uniform_filter1d
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# parameters
n_roll = 80 # rolling mean over n years around a given point in time

# paths to ps and mse files
path_noisy = sys.argv[1]
path_smooth = sys.argv[2]
varname = sys.argv[3]

# open files
print('Loading data...')
file_noisy = Dataset(path_noisy, 'r')
print('Done.\n')

# read data
varname = ''.join([i for i in varname if not i.isdigit()])
nd = file_noisy.variables[varname][:]
nd = nd.filled(fill_value=np.nan)

# smoothen data by taking the rolling mean
sd = uniform_filter1d(nd, n_roll, axis=0)

# save file as netCDF
file_smooth = Dataset(path_smooth, "w", format='NETCDF4_CLASSIC')

# copy attributes from original file
file_smooth.setncatts(file_noisy.__dict__)

# copy dimensions from original file
for name, dimension in file_noisy.dimensions.items():
    # if any(name in s for s in ['lon']):
    #     continue
    file_smooth.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

# copy all variables except time from ps file
for name, variable in file_noisy.variables.items():
    if any(name in s for s in [varname]):
        continue
    
    x = file_smooth.createVariable(name, variable.datatype, variable.dimensions)
    file_smooth[name].setncatts(file_noisy[name].__dict__)
    file_smooth[name][:] = file_noisy[name][:]
    
noisy = file_noisy.variables[varname]
smooth = file_smooth.createVariable('%s_sm' % (varname), 'f4', noisy.dimensions)
smooth.units = noisy.units
smooth.long_name = ("smoothed %s") % (noisy.long_name)
smooth[:] = sd
file_smooth.close()

