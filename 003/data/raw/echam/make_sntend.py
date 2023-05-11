import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps and sni files
path_sni = sys.argv[1]
path_sntend = sys.argv[2]

# open files
file_sni = Dataset(path_sni, 'r')

# read data
sni = file_sni.variables['sni'][:] # (time x lev x lat x lon)

# take time sntendency
dsnidt = np.empty(sni.shape) 
dsnidt[0] = (sni[1,:,:]-sni[0,:,:])/(360*86400/12)
dsnidt[-1] = (sni[-1,:,:]-sni[-2,:,:])/(360*86400/12)

# -- center difference
dsnidt[1:-1,:,:] = (sni[2:,:,:]-sni[:-2,:,:])/(2*360*86400/12)
# -- forward difference
# dsnidt[1:-1,:,:] = (sni[2:,:,:]-sni[1:-1,:,:])/(360*86400/12)

# save file as netCDF
file_sntend = Dataset(path_sntend, "w", format='NETCDF4_CLASSIC')

# copy attributes from sni file
file_sntend.setncatts(file_sni.__dict__)

# copy dimensions from sni file
for name, dimension in file_sni.dimensions.items():
    file_sntend.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_sni.variables.items():
    if any(name in s for s in ['sni']):
        continue
    
    x = file_sntend.createVariable(name, variable.datatype, variable.dimensions)
    file_sntend[name].setncatts(file_sni[name].__dict__)
    file_sntend[name][:] = file_sni[name][:]
    
sntend = file_sntend.createVariable('sntend', 'f4', ("time","lat","lon"))
sntend.units = "m s**-1"
sntend.long_name = "time sntendency of snow equivalent water depth"
sntend[:,:,:] = dsnidt

file_sntend.close()
