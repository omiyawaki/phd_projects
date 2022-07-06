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
path_ref = sys.argv[1]
path_fice = sys.argv[2]
path_sic = sys.argv[3]

# open files
file_ref = Dataset(path_ref, 'r')
file_fice = Dataset(path_fice, 'r')

# read data
ref = file_ref.variables['sic'][:] # (mon x lat x lon)
fice = file_fice.variables['seaice'][:] # (mon x lat x lon)

# save file as netCDF
file_sic = Dataset(path_sic, "w", format='NETCDF4_CLASSIC')

# copy attributes from ref file
file_sic.setncatts(file_ref.__dict__)

# copy dimensions from ref file
for name, dimension in file_ref.dimensions.items():
    file_sic.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from ref file
for name, variable in file_ref.variables.items():
    if any(name in s for s in ['sic']):
        continue
    
    x = file_sic.createVariable(name, variable.datatype, variable.dimensions)
    file_sic[name].setncatts(file_ref[name].__dict__)
    file_sic[name][:] = file_ref[name][:]
    
sic = file_sic.createVariable('sic', 'f4', ("time","lat","lon"))
sic.units = 'unitless'
sic.long_name = 'sea ice concentration'
sic[:] = fice

file_sic.close()
