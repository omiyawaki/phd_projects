import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps and tsw files
path_tsw = sys.argv[1]
path_twtend = sys.argv[2]

# open files
file_tsw = Dataset(path_tsw, 'r')

# read data
tsw = file_tsw.variables['tsw'][:] # (time x lev x lat x lon)

# take time twtendency
dtswdt = np.empty(tsw.shape) 
dtswdt[1:-1,:,:] = (tsw[2:,:,:]-tsw[:-2,:,:])/(2*365*86400/12)
dtswdt[0] = (tsw[1,:,:]-tsw[0,:,:])/(365*86400/12)
dtswdt[-1] = (tsw[-1,:,:]-tsw[-2,:,:])/(365*86400/12)

# save file as netCDF
file_twtend = Dataset(path_twtend, "w", format='NETCDF4_CLASSIC')

# copy attributes from tsw file
file_twtend.setncatts(file_tsw.__dict__)

# copy dimensions from tsw file
for name, dimension in file_tsw.dimensions.items():
    file_twtend.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_tsw.variables.items():
    if any(name in s for s in ['tsw']):
        continue
    
    x = file_twtend.createVariable(name, variable.datatype, variable.dimensions)
    file_twtend[name].setncatts(file_tsw[name].__dict__)
    file_twtend[name][:] = file_tsw[name][:]
    
twtend = file_twtend.createVariable('twtend', 'f4', ("time","lat","lon"))
twtend.units = "K s**-1"
twtend.long_name = "time tendency of water temperature"
twtend[:,:,:] = dtswdt

file_twtend.close()
