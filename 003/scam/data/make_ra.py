import sys
import numpy as np
from scipy import interpolate, integrate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to input and output files
path_in = sys.argv[1]
path_out = sys.argv[2]

# open files
print('Loading input file...')
file_in = Dataset(path_in, 'r')
print('Done.\n')

# read data
FSNT = file_in.variables['FSNT'][:] # (mon)
FLNT = file_in.variables['FLNT'][:] # (mon)
FSNS = file_in.variables['FSNS'][:] # (mon)
FLNS = file_in.variables['FLNS'][:] # (mon)

# compute radiative cooling
racool = FSNT - FLNT - FSNS + FLNS

# save file as netCDF
file_out = Dataset(path_out, "w", format='NETCDF4_CLASSIC')

# copy attributes from in file
file_out.setncatts(file_in.__dict__)

# copy dimensions from in file
for name, dimension in file_in.dimensions.items():
    file_out.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

# copy all variables corresponding to dimensions from in file
for name, variable in file_in.variables.items():
    if any(name in s for s in ['time', 'bnds', 'lat', 'lon', 'lev', 'ilev']):
        x = file_out.createVariable(name, variable.datatype, variable.dimensions)
        file_out[name].setncatts(file_in[name].__dict__)
        file_out[name][:] = file_in[name][:]
    
ra = file_out.createVariable('ra', 'f4', ("time", "lat", "lon"))
ra.units = "W m**-2"
ra.long_name = "radiative cooling"
ra[:] = racool

file_out.close()
