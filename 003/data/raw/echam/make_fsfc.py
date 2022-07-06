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
path_ahfl = sys.argv[1]
path_ahfs = sys.argv[2]
path_trads = sys.argv[3]
path_srads = sys.argv[4]
path_fsfc = sys.argv[5]

# open files
file_ahfl = Dataset(path_ahfl, 'r')
file_ahfs = Dataset(path_ahfs, 'r')
file_srads = Dataset(path_srads, 'r')
file_trads = Dataset(path_trads, 'r')

# read data
ahfl = file_ahfl.variables['ahfl'][:] # (mon x lat x lon)
ahfs = file_ahfs.variables['ahfs'][:] # (mon x lat x lon)
srads = file_srads.variables['srads'][:] # (mon x lat x lon)
trads = file_trads.variables['trads'][:] # (mon x lat x lon)

# compute surface turbulent fluxes
sfc = ahfl + ahfs + trads + srads

# save file as netCDF
file_fsfc = Dataset(path_fsfc, "w", format='NETCDF4_CLASSIC')

# copy attributes from ahfl file
file_fsfc.setncatts(file_ahfl.__dict__)

# copy dimensions from ahfl file
for name, dimension in file_ahfl.dimensions.items():
    file_fsfc.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from ahfl file
for name, variable in file_ahfl.variables.items():
    if any(name in s for s in ['ahfl']):
        continue
    
    x = file_fsfc.createVariable(name, variable.datatype, variable.dimensions)
    file_fsfc[name].setncatts(file_ahfl[name].__dict__)
    file_fsfc[name][:] = file_ahfl[name][:]
    
fsfc = file_fsfc.createVariable('fsfc', 'f4', ("time","lat","lon"))
fsfc.units = "W m**-2"
fsfc.long_name = "net surface heat flux"
fsfc[:] = sfc

file_fsfc.close()
