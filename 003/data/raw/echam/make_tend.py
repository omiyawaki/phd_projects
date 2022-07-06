import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps and mse files
path_mse = sys.argv[1]
path_ps = sys.argv[2]
path_tend = sys.argv[3]

# open files
file_mse = Dataset(path_mse, 'r')
file_ps = Dataset(path_ps, 'r')

# read data
mse = file_mse.variables['mse'][:] # (time x lev x lat x lon)
ps = file_ps.variables['aps'][:] # (time x lat x lon)
try:
    lev = file_mse.variables['plev'][:] # (lev)
except:
    lev = file_mse.variables['lev'][:] # (lev)

# compute pressure field
ps = np.transpose(np.tile(ps, [len(lev),1,1,1]), [1,0,2,3])
pa = np.transpose(np.tile(lev, [mse.shape[0],mse.shape[2], mse.shape[3],1]), [0,3,1,2])

# take time tendency
# 3d mse tendency
dmsedt = np.empty(mse.shape) 
dmsedt[1:-1,:,:,:] = (mse[2:,:,:,:]-mse[0:-2,:,:,:])/(2*365*86400/12)
dmsedt[0] = (mse[1,:,:,:]-mse[0,:,:,:])/(365*86400/12)
dmsedt[-1] = (mse[-1,:,:,:]-mse[-2,:,:,:])/(365*86400/12)

# set data below surface to 0 (will not be counted in vertical integral)
below = pa > ps
dmsedt[below] = 0

# compute vertical integral
if pa[0,1,0,0]-pa[0,0,0,0]>0: # if pressure increases with index
    dvmsedt = 1/g*np.trapz(dmsedt, pa, axis=1)
else:
    dvmsedt = -1/g*np.trapz(dmsedt, pa, axis=1)

# save file as netCDF
file_tend = Dataset(path_tend, "w", format='NETCDF4_CLASSIC')

# copy attributes from mse file
file_tend.setncatts(file_mse.__dict__)

# copy dimensions from mse file
for name, dimension in file_mse.dimensions.items():
    if any(name in s for s in ['plev']):
        continue
    file_tend.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_mse.variables.items():
    if any(name in s for s in ['mse' 'plev']):
        continue
    
    x = file_tend.createVariable(name, variable.datatype, variable.dimensions)
    file_tend[name].setncatts(file_mse[name].__dict__)
    file_tend[name][:] = file_mse[name][:]
    
tend = file_tend.createVariable('tend', 'f4', ("time","lat","lon"))
tend.units = "W m**-2"
tend.long_name = "time tendency of vertically integrated moist static energy"
tend[:,:,:] = dvmsedt

file_tend.close()
