import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;
Lf = 3.337e5 # J kg**-1 latent heat of fusion
rhoi = 917 # kg m**-3 density of ice

# paths to ps and siced files
path_tsurf = sys.argv[1]
path_tsi = sys.argv[2]
path_siced = sys.argv[3]
path_stend = sys.argv[4]

# open files
file_tsurf = Dataset(path_tsurf, 'r')
file_tsi = Dataset(path_tsi, 'r')
file_siced = Dataset(path_siced, 'r')

# read data
tsurf = file_tsurf.variables['tsurf'][:] # (time x lev x lat x lon)
tsi = file_tsi.variables['tsi'][:] # (time x lev x lat x lon)
siced = file_siced.variables['siced'][:] # (time x lev x lat x lon)

if tsurf.shape[0]*30 != siced.shape[0]:
    error('This script uses daily sea ice thickness data. Please check that sea ice contains daily data.')

# take time stendency
# dt = 360*86400/12
# dt = 86400
dsiceddt = np.empty(siced.shape) 
dsiceddt[0] = (siced[1,:,:]-siced[0,:,:])#/dt
dsiceddt[-1] = (siced[-1,:,:]-siced[-2,:,:])#/dt

# -- center difference
dsiceddt[1:-1,:,:] = (siced[2:,:,:]-siced[:-2,:,:])/2
# -- forward difference
# dsiceddt[1:-1,:,:] = (siced[2:,:,:]-siced[1:-1,:,:])#/(dt)

# keep only regions of melt (dhdt < 0)
dsiceddt = np.minimum(np.zeros_like(dsiceddt), dsiceddt)
# nomelt=np.where(dsiceddt > 0)
# dsiceddt[nomelt] = 0

# sum total melt each month
fmelt=np.empty_like(tsurf)
for m in tqdm(range(tsurf.shape[0])):
    db=30*m # index of first day of this month
    de=30*(m+1) # index of first day of next month
    fmelt[m,...]=-Lf*rhoi*np.nansum(dsiceddt[db:de,...],axis=0)/(30*86400)

# fmelt = np.maximum(np.zeros_like(fmelt), fmelt)

# save file as netCDF
file_stend = Dataset(path_stend, "w", format='NETCDF4_CLASSIC')

# copy attributes from tsurf file
file_stend.setncatts(file_tsurf.__dict__)

# copy dimensions from tsurf file
for name, dimension in file_tsurf.dimensions.items():
    file_stend.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_tsurf.variables.items():
    if any(name in s for s in ['tsurf']):
        continue
    
    x = file_stend.createVariable(name, variable.datatype, variable.dimensions)
    file_stend[name].setncatts(file_tsurf[name].__dict__)
    file_stend[name][:] = file_tsurf[name][:]
    
stend = file_stend.createVariable('smelt', 'f4', ("time","lat","lon"))
stend.units = "W m**-2"
stend.long_name = "Sea ice melt flux (diagnosed)"
stend[:,:,:] = fmelt

file_stend.close()
