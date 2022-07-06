import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to mse
path_in = sys.argv[1]
path_out = sys.argv[2]
llo = float(sys.argv[3]) # lower lat bound of area mean
lhi = float(sys.argv[4]) # upper lat bound of area mean
varname = sys.argv[5] # variable name

# create lat vector for area averaging
nlat = 101
lati = np.linspace(llo, lhi, nlat)
rlati = np.radians(lati)
clati = np.cos(rlati)

# open input file
file_in = Dataset(path_in, 'r')

# read data
varname = ''.join([i for i in varname if not i.isdigit()])
vin = np.squeeze(file_in.variables[varname][:]) # (mon x lat)
lat = np.squeeze(file_in.variables['lat'][:]) # (mon x lat)

# for datasets that fill data below surface as missing data, fill with nans
vin = vin.filled(fill_value=np.nan)

# take area mean
fint = interpolate.interp1d(lat, vin, axis=1, bounds_error=False)
vini = fint(lati)
vina = np.nansum(clati*vini,axis=1) / np.nansum(clati)

# save file as netCDF
file_out = Dataset(path_out, "w", format='NETCDF4_CLASSIC')

# copy dimensions 
for name, dimension in file_in.dimensions.items():
    file_out.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables
for name, variable in file_in.variables.items():
    if any(name in s for s in [varname]):
        continue
    
    x = file_out.createVariable(name, variable.datatype, variable.dimensions)
    file_out[name].setncatts(file_in[name].__dict__)
    file_out[name][:] = file_in[name][:]
    
outa = file_out.createVariable(varname, 'f4', ("time"))
outa.units = file_in[varname].units
outa.long_name = file_in[varname].long_name
outa[:] = vina

file_out.close()
