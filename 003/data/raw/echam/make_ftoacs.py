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
path_traf0 = sys.argv[1]
path_sraf0 = sys.argv[2]
path_ftoacs = sys.argv[3]

# open files
file_sraf0 = Dataset(path_sraf0, 'r')
file_traf0 = Dataset(path_traf0, 'r')

# read data
sraf0 = file_sraf0.variables['sraf0'][:] # (mon x lat x lon)
traf0 = file_traf0.variables['traf0'][:] # (mon x lat x lon)

# compute net TOA rafiative flux
toacs = traf0 + sraf0

# save file as netCDF
file_ftoacs = Dataset(path_ftoacs, "w", format='NETCDF4_CLASSIC')

# copy attributes from traf0 file
file_ftoacs.setncatts(file_traf0.__dict__)

# copy dimensions from traf0 file
for name, dimension in file_traf0.dimensions.items():
    file_ftoacs.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from traf0 file
for name, variable in file_traf0.variables.items():
    if any(name in s for s in ['traf0']):
        continue
    
    x = file_ftoacs.createVariable(name, variable.datatype, variable.dimensions)
    file_ftoacs[name].setncatts(file_traf0[name].__dict__)
    file_ftoacs[name][:] = file_traf0[name][:]
    
ftoacs = file_ftoacs.createVariable('ftoacs', 'f4', ("time","lat","lon"))
ftoacs.units = "W m**-2"
ftoacs.long_name = "net TOA rafiative flux (clear sky)"
ftoacs[:] = toacs

file_ftoacs.close()
