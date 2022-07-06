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
path_a32 = sys.argv[2]
path_qflux = sys.argv[3]

# open files
file_ref = Dataset(path_ref, 'r')
file_a32 = Dataset(path_a32, 'r')

# read data
ref = file_ref.variables['aflux'][:] # (mon x lat x lon)
a32 = np.flip(np.squeeze(file_a32.variables['q_flux_32'][:])) # (mon x lat x lon) qflux

a32 = np.transpose( np.tile(a32, [ref.shape[0],ref.shape[2],1]), [0,2,1] )

# save file as netCDF
file_qflux = Dataset(path_qflux, "w", format='NETCDF4_CLASSIC')

# copy attributes from ref file
file_qflux.setncatts(file_ref.__dict__)

# copy dimensions from ref file
for name, dimension in file_ref.dimensions.items():
    file_qflux.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from ref file
for name, variable in file_ref.variables.items():
    if any(name in s for s in ['aflux']):
        continue
    
    x = file_qflux.createVariable(name, variable.datatype, variable.dimensions)
    file_qflux[name].setncatts(file_ref[name].__dict__)
    file_qflux[name][:] = file_ref[name][:]
    
qflux = file_qflux.createVariable('aflux', 'f4', ("time","lat","lon"))
qflux.units = "W m**-2"
qflux.long_name = "LW flux over water"
qflux[:] = a32

file_qflux.close()
