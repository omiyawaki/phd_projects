import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;
d = 40
di = 0.05
ci = 2106
cs = 2090
rhoi = 917 # kg m**-3 density of ice
rhos = 300 # density of snow
rhow = 1000 # density of liquid water
ki = 2.1656 # ice thermal conductivity [W/K/m]
ks = 0.31 # snow thermal conductivity [W/K/m]
T0 = 273.15-1.8 # sea ice freezing temperature

# paths to ps and tsurf files
path_tsurf = sys.argv[1]
path_tsi = sys.argv[2]
path_siced = sys.argv[3]
path_sni = sys.argv[4]
path_cond = sys.argv[5]

# open files
file_tsurf = Dataset(path_tsurf, 'r')
file_tsi = Dataset(path_tsi, 'r')
file_siced = Dataset(path_siced, 'r')
file_sni = Dataset(path_sni, 'r')

# read data
tsurf = file_tsurf.variables['tsurf'][:] # (time x lev x lat x lon)
tsi = file_tsi.variables['tsi'][:] # (time x lev x lat x lon)
siced = file_siced.variables['siced'][:] # (time x lev x lat x lon)
sni = file_sni.variables['sni'][:] # (time x lev x lat x lon)

# diagnose conductive flux
heff=siced+ki*rhow/(ks*rhos)*sni
fcon=ki*(tsi-T0)/heff

# save file as netCDF
file_cond = Dataset(path_cond, "w", format='NETCDF4_CLASSIC')

# copy attributes from tsurf file
file_cond.setncatts(file_tsurf.__dict__)

# copy dimensions from tsurf file
for name, dimension in file_tsurf.dimensions.items():
    file_cond.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_tsurf.variables.items():
    if any(name in s for s in ['tsurf']):
        continue
    
    x = file_cond.createVariable(name, variable.datatype, variable.dimensions)
    file_cond[name].setncatts(file_tsurf[name].__dict__)
    file_cond[name][:] = file_tsurf[name][:]
    
cond = file_cond.createVariable('cond', 'f4', ("time","lat","lon"))
cond.units = "W m**-2"
cond.long_name = "conductive heat flux through sea ice"
cond[:,:,:] = fcon

file_cond.close()
