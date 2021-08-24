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
path_va = sys.argv[1]
path_vas = sys.argv[2]

# open files
file_va = Dataset(path_va, 'r')

# read data
va = file_va.variables['v'][:] # (day x lev x lat x lon)
lev = file_va.variables['lev'][:] # (lev)

# vs = np.empty([va.shape[0],va.shape[2],va.shape[3]])
# extract surface va
if lev[1]-lev[0]>0:
    vs = va[:,-1,:,:]
else:
    vs = va[:,0,:,:]

# save file as netCDF
file_vas = Dataset(path_vas, "w", format='NETCDF4_CLASSIC')

# copy dimensions from va file
for name, dimension in file_va.dimensions.items():
    if any(name in s for s in ['lev', 'hyai', 'hybi', 'hyam', 'hybm']):
        continue

    file_vas.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables from va file
for name, variable in file_va.variables.items():
    if any(name in s for s in ['aps', 'u', 'v', 't', 'q', 'omega', 'geopoth',  'lev', 'hyai', 'hybi', 'hyam', 'hybm']):
        continue
    
    x = file_vas.createVariable(name, variable.datatype, variable.dimensions)
    file_vas[name].setncatts(file_va[name].__dict__)
    file_vas[name][:] = file_va[name][:]
    
vas = file_vas.createVariable('vas', 'f4', ("time","lat","lon"))
vas.units = "m/s"
vas.long_name = "surface northward velocity"
vas[:,:,:] = vs

file_va.close()
