import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps and siced files
path_siced = sys.argv[1]
path_stend = sys.argv[2]

# open files
file_siced = Dataset(path_siced, 'r')

# read data
siced = file_siced.variables['siced'][:] # (time x lev x lat x lon)

# take time stendency
dsiceddt = np.empty(siced.shape) 
dsiceddt[0] = (siced[1,:,:]-siced[0,:,:])/(360*86400/12)
dsiceddt[-1] = (siced[-1,:,:]-siced[-2,:,:])/(360*86400/12)

# -- center difference
dsiceddt[1:-1,:,:] = (siced[2:,:,:]-siced[:-2,:,:])/(2*360*86400/12)
# -- forward difference
# dsiceddt[1:-1,:,:] = (siced[2:,:,:]-siced[1:-1,:,:])/(360*86400/12)

# save file as netCDF
file_stend = Dataset(path_stend, "w", format='NETCDF4_CLASSIC')

# copy attributes from siced file
file_stend.setncatts(file_siced.__dict__)

# copy dimensions from siced file
for name, dimension in file_siced.dimensions.items():
    file_stend.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_siced.variables.items():
    if any(name in s for s in ['siced']):
        continue
    
    x = file_stend.createVariable(name, variable.datatype, variable.dimensions)
    file_stend[name].setncatts(file_siced[name].__dict__)
    file_stend[name][:] = file_siced[name][:]
    
stend = file_stend.createVariable('stend', 'f4', ("time","lat","lon"))
stend.units = "m s**-1"
stend.long_name = "time stendency of sea ice thickness"
stend[:,:,:] = dsiceddt

file_stend.close()
