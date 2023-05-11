import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps and tsurf files
path_tsurf = sys.argv[1]
path_ttend = sys.argv[2]

# open files
file_tsurf = Dataset(path_tsurf, 'r')

# read data
tsurf = file_tsurf.variables['tsurf'][:] # (time x lev x lat x lon)

# take time ttendency
dtsurfdt = np.empty(tsurf.shape) 
dtsurfdt[1:-1,:,:] = (tsurf[2:,:,:]-tsurf[:-2,:,:])/(2*360*86400/12)
dtsurfdt[0] = (tsurf[1,:,:]-tsurf[0,:,:])/(360*86400/12)
dtsurfdt[-1] = (tsurf[-1,:,:]-tsurf[-2,:,:])/(360*86400/12)

# save file as netCDF
file_ttend = Dataset(path_ttend, "w", format='NETCDF4_CLASSIC')

# copy attributes from tsurf file
file_ttend.setncatts(file_tsurf.__dict__)

# copy dimensions from tsurf file
for name, dimension in file_tsurf.dimensions.items():
    file_ttend.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_tsurf.variables.items():
    if any(name in s for s in ['tsurf']):
        continue
    
    x = file_ttend.createVariable(name, variable.datatype, variable.dimensions)
    file_ttend[name].setncatts(file_tsurf[name].__dict__)
    file_ttend[name][:] = file_tsurf[name][:]
    
ttend = file_ttend.createVariable('ttend', 'f4', ("time","lat","lon"))
ttend.units = "K s**-1"
ttend.long_name = "time tendency of surface temperature"
ttend[:,:,:] = dtsurfdt

file_ttend.close()
