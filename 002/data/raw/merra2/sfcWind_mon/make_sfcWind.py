import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for MERRA2

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps, ta, zg, and hus files
path_uas = sys.argv[1]
path_vas = sys.argv[2]
path_sfcWind = sys.argv[3]

# open files
file_uas = Dataset(path_uas, 'r')
file_vas = Dataset(path_vas, 'r')

# read data
uas = file_uas.variables['U10M'][:] # (mon x lat x lon)
vas = file_vas.variables['V10M'][:] # (day x lev x lat x lon)

# calculate horizontal wind speed
uv10m = np.sqrt(uas**2 + vas**2)

# save file as netCDF
file_sfcWind = Dataset(path_sfcWind, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_sfcWind.setncatts(file_uas.__dict__)

# copy dimensions from uas file
for name, dimension in file_uas.dimensions.items():
    file_sfcWind.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_uas.variables.items():
    if any(name in s for s in ['U10M', 'V10M']):
        continue
    
    x = file_sfcWind.createVariable(name, variable.datatype, variable.dimensions)
    file_sfcWind[name].setncatts(file_uas[name].__dict__)
    file_sfcWind[name][:] = file_uas[name][:]
    
sfcWind = file_sfcWind.createVariable('sfcWind', 'f4', ("time","lat","lon"))
sfcWind.units = "m/s"
sfcWind.long_name = "10 m wind speed"
sfcWind[:,:,:] = uv10m

file_sfcWind.close()