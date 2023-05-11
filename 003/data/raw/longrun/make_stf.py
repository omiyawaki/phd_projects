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
path_hfls = sys.argv[1]
path_hfss = sys.argv[2]
path_stf = sys.argv[3]

# open files
print('Loading hfls file...')
file_hfls = Dataset(path_hfls, 'r')
print('Done.\n')

print('Loading hfss file...')
file_hfss = Dataset(path_hfss, 'r')
print('Done.\n')

# read data
hfls = file_hfls.variables['hfls'][:] # (mon x lat x lon)
hfss = file_hfss.variables['hfss'][:] # (mon x lat x lon)

# compute surface turbulent fluxes
lhsh = hfls + hfss

# save file as netCDF
file_stf = Dataset(path_stf, "w", format='NETCDF4_CLASSIC')

# copy attributes from hfls file
file_stf.setncatts(file_hfls.__dict__)

# copy dimensions from hfls file
for name, dimension in file_hfls.dimensions.items():
    file_stf.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from hfls file
for name, variable in file_hfls.variables.items():
    if any(name in s for s in ['hfls']):
        continue
    
    x = file_stf.createVariable(name, variable.datatype, variable.dimensions)
    file_stf[name].setncatts(file_hfls[name].__dict__)
    file_stf[name][:] = file_hfls[name][:]
    
stf = file_stf.createVariable('stf', 'f4', ("time","lat","lon"))
stf.units = "W m**-2"
stf.long_name = "surface latent plus sensible heat flux"
stf[:] = lhsh

file_stf.close()
