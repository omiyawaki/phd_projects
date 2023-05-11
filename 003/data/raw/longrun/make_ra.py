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
path_rsdt = sys.argv[1]
path_rsut = sys.argv[2]
path_rlut = sys.argv[3]
path_rsds = sys.argv[4]
path_rsus = sys.argv[5]
path_rlds = sys.argv[6]
path_rlus = sys.argv[7]
path_ra = sys.argv[8]

# open files
file_rsdt = Dataset(path_rsdt, 'r')
file_rsut = Dataset(path_rsut, 'r')
file_rlut = Dataset(path_rlut, 'r')
file_rsds = Dataset(path_rsds, 'r')
file_rsus = Dataset(path_rsus, 'r')
file_rlds = Dataset(path_rlds, 'r')
file_rlus = Dataset(path_rlus, 'r')

# read data
rsdt = file_rsdt.variables['rsdt'][:] # (mon x lat x lon)
rsut = file_rsut.variables['rsut'][:] # (mon x lat x lon)
rlut = file_rlut.variables['rlut'][:] # (mon x lat x lon)
rsds = file_rsds.variables['rsds'][:] # (mon x lat x lon)
rsus = file_rsus.variables['rsus'][:] # (mon x lat x lon)
rlds = file_rlds.variables['rlds'][:] # (mon x lat x lon)
rlus = file_rlus.variables['rlus'][:] # (mon x lat x lon)

# compute radiative cooling
racool = rsdt - rsut - rlut - rsds + rsus - rlds + rlus 

# save file as netCDF
file_ra = Dataset(path_ra, "w", format='NETCDF4_CLASSIC')

# copy attributes from rsdt file
file_ra.setncatts(file_rsdt.__dict__)

# copy dimensions from rsdt file
for name, dimension in file_rsdt.dimensions.items():
    file_ra.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from rsdt file
for name, variable in file_rsdt.variables.items():
    if any(name in s for s in ['rsdt']):
        continue
    
    x = file_ra.createVariable(name, variable.datatype, variable.dimensions)
    file_ra[name].setncatts(file_rsdt[name].__dict__)
    file_ra[name][:] = file_rsdt[name][:]
    
ra = file_ra.createVariable('ra', 'f4', ("time","lat","lon"))
ra.units = "W m**-2"
ra.long_name = "radiative cooling"
ra[:] = racool

file_ra.close()
