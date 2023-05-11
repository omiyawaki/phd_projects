import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to tendv
path_tendv = sys.argv[1]
path_tendvz = sys.argv[2]

# open files
file_tendv = Dataset(path_tendv, 'r')

# read data
tendv = np.squeeze(file_tendv.variables['tendv'][:]) # (mon x lev x lat x lon)

# for datasets that fill data below surface as missing data, fill with nans
tendv = tendv.filled(fill_value=np.nan)

# take zonal mean
# mz = np.transpose( np.tile(np.nanmean(tendv,axis=2), [tendv.shape[2],1,1]), [1,2,0] )
mz = np.nanmean(tendv,axis=3)

# save file as netCDF
file_tendvz = Dataset(path_tendvz, "w", format='NETCDF4_CLASSIC')

# # copy attributes from zg file
# file_tendvz.setncatts(file_zg.__dict__)

# copy dimensions from tendv file
for name, dimension in file_tendv.dimensions.items():
    file_tendvz.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables
for name, variable in file_tendv.variables.items():
    if any(name in s for s in ['tendv']):
        continue
    
    x = file_tendvz.createVariable(name, variable.datatype, variable.dimensions)
    file_tendvz[name].setncatts(file_tendv[name].__dict__)
    file_tendvz[name][:] = file_tendv[name][:]
    
try:
    tendvz = file_tendvz.createVariable('tendv', 'f4', ("time","plev","lat"))
except:
    tendvz = file_tendvz.createVariable('tendv', 'f4', ("time","lev","lat"))
tendvz.units = "W/kg"
tendvz.long_name = "time tendency of moist static energy"
tendvz[:] = mz

file_tendvz.close()
