import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;
p0 = 1100.e2 # bottom integral bound [hPa]
pt = 0. # top integral bound [hPa]

# Paths to ps and mse files
path_mmc = sys.argv[1]
path_lmmc = sys.argv[2]
path_ummc = sys.argv[3]
path_res = sys.argv[4]

# Open files
file_mmc = Dataset(path_mmc, 'r')
file_lmmc = Dataset(path_lmmc, 'r')
file_ummc = Dataset(path_ummc, 'r')

# Read data
mmc = file_mmc.variables['vmmmc'][:] # (mon x lat)
lmmc = file_lmmc.variables['lvmmmc'][:] # (mon x lat)
ummc = file_ummc.variables['uvmmmc'][:] # (mon x lat)

# residual (difference between full vertical integral and sum of split integrals)
r = mmc - (lmmc + ummc)

# save file as netCDF
file_res = Dataset(path_res, "w", format='NETCDF4_CLASSIC')

# copy attributes from mmc file
file_res.setncatts(file_mmc.__dict__)

# copy dimensions from mmc file
for name, dimension in file_mmc.dimensions.items():
    if any(name in s for s in ['plev']):
        continue

    file_res.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from mmc file
for name, variable in file_mmc.variables.items():
    if any(name in s for s in ['vmmmc']):
        continue
    
    x = file_res.createVariable(name, variable.datatype, variable.dimensions)
    file_res[name].setncatts(file_mmc[name].__dict__)
    file_res[name][:] = file_mmc[name][:]

res = file_res.createVariable('res', 'f4', ("time","lat"))
res.units = "W"
res.long_name = "difference between fully vertically integrated MMC transport and sum of split (lower and upper) integrals"
res[:,:] = r

file_res.close()
