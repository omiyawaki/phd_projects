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
path_out = sys.argv[2]

# open files
file_ref = Dataset(path_ref, 'r')

# read data
ref = file_ref.variables['CO2'][:] # (time)

pic = 8*284.725*np.ones_like(ref)

# save file as netCDF
file_out = Dataset(path_out, "w", format='NETCDF4_CLASSIC')

# copy attributes from ref file
file_out.setncatts(file_ref.__dict__)

# copy dimensions from ref file
for name, dimension in file_ref.dimensions.items():
    file_out.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from ref file
for name, variable in file_ref.variables.items():
    if any(name in s for s in ['CO2']):
        continue
    
    x = file_out.createVariable(name, variable.datatype, variable.dimensions)
    file_out[name].setncatts(file_ref[name].__dict__)
    file_out[name][:] = file_ref[name][:]
    
out = file_out.createVariable('CO2', 'f4', ("time"))
out.units = "1e-6 mol/mol"
out.long_name = "CO2"
out[:] = pic

file_out.close()
