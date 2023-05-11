import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps and tsi files
path_tsi = sys.argv[1]
path_titend = sys.argv[2]

# open files
file_tsi = Dataset(path_tsi, 'r')

# read data
tsi = file_tsi.variables['tsi'][:] # (time x lev x lat x lon)

# take time titendency
dtsidt = np.empty(tsi.shape) 
dtsidt[1:-1,:,:] = (tsi[2:,:,:]-tsi[:-2,:,:])/(2*360*86400/12)
dtsidt[0] = (tsi[1,:,:]-tsi[0,:,:])/(360*86400/12)
dtsidt[-1] = (tsi[-1,:,:]-tsi[-2,:,:])/(360*86400/12)

# save file as netCDF
file_titend = Dataset(path_titend, "w", format='NETCDF4_CLASSIC')

# copy attributes from tsi file
file_titend.setncatts(file_tsi.__dict__)

# copy dimensions from tsi file
for name, dimension in file_tsi.dimensions.items():
    file_titend.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_tsi.variables.items():
    if any(name in s for s in ['tsi']):
        continue
    
    x = file_titend.createVariable(name, variable.datatype, variable.dimensions)
    file_titend[name].setncatts(file_tsi[name].__dict__)
    file_titend[name][:] = file_tsi[name][:]
    
titend = file_titend.createVariable('titend', 'f4', ("time","lat","lon"))
titend.units = "K s**-1"
titend.long_name = "time tendency of ice temperature"
titend[:,:,:] = dtsidt

file_titend.close()
