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
path_bot = sys.argv[1]
path_tend = sys.argv[2]
path_aht = sys.argv[3]

# open files
file_bot = Dataset(path_bot, 'r')
file_tend = Dataset(path_tend, 'r')

# read data
srad0 = file_bot.variables['srad0'][:] # (mon x lat x lon)
trad0 = file_bot.variables['trad0'][:] # (mon x lat x lon)
srads = file_bot.variables['srads'][:] # (mon x lat x lon)
trads = file_bot.variables['trads'][:] # (mon x lat x lon)
ahfl = file_bot.variables['ahfl'][:] # (mon x lat x lon)
ahfs = file_bot.variables['ahfs'][:] # (mon x lat x lon)
tend = file_tend.variables['tend'][:] # (mon x lat x lon)
lat = file_tend.variables['lat'][:] # (lat)
rlat = np.radians(lat)
clat = np.cos(rlat)

# infer total energy flux divergence as the residual
fluxdiv = srad0 + trad0 - srads - trads - ahfl - ahfs - tend

# zonal mean
fluxdiv_z = np.nanmean(fluxdiv, axis=2)

# subtract global mean
fluxdiv_g = np.trapz(clat*fluxdiv_z, rlat, axis=1) / np.trapz(clat, rlat)
fluxdiv_z = fluxdiv_z - np.transpose(np.tile(fluxdiv_g, [len(lat), 1]))

# meridionally integrate
vE = 2*np.pi*a**2 * integrate.cumtrapz(clat*fluxdiv_z, rlat, axis=1, initial=0)
# if lat[1]-lat[0]<0:
#     vE = -vE

# save file as netCDF
file_aht = Dataset(path_aht, "w", format='NETCDF4_CLASSIC')

# copy attributes from tend file
file_aht.setncatts(file_tend.__dict__)

# copy dimensions from tend file
for name, dimension in file_tend.dimensions.items():
    file_aht.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_tend.variables.items():
    if any(name in s for s in ['tend']):
        continue
    
    x = file_aht.createVariable(name, variable.datatype, variable.dimensions)
    file_aht[name].setncatts(file_tend[name].__dict__)
    file_aht[name][:] = file_tend[name][:]
    
aht = file_aht.createVariable('aht', 'f4', ("time","lat"))
aht.units = "W"
aht.long_name = "vertically integrated total energy flux transport"
aht[:,:] = vE

file_aht.close()
