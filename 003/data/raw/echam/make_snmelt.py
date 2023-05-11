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
rhoh2o = 1000 # kg m**-3 density of liquid water
rhos = 300 # kg m**-3 density of snow

# paths to ps and sni files
path_tsurf = sys.argv[1]
path_tsi = sys.argv[2]
path_sni = sys.argv[3]
path_stend = sys.argv[4]

# open files
file_tsurf = Dataset(path_tsurf, 'r')
file_tsi = Dataset(path_tsi, 'r')
file_sni = Dataset(path_sni, 'r')

# read data
tsurf = file_tsurf.variables['tsurf'][:] # (time x lev x lat x lon)
tsi = file_tsi.variables['tsi'][:] # (time x lev x lat x lon)
sni = file_sni.variables['sni'][:] # (time x lev x lat x lon)

if tsurf.shape[0]*30 != sni.shape[0]:
    error('This script uses daily sea ice thickness data. Please check that sea ice contains daily data.')

# take time stendency
# dt = 360*86400/12
# dt = 86400
dsnidt = np.empty(sni.shape) 
dsnidt[0] = (sni[1,:,:]-sni[0,:,:])#/dt
dsnidt[-1] = (sni[-1,:,:]-sni[-2,:,:])#/dt

# -- center difference
dsnidt[1:-1,:,:] = (sni[2:,:,:]-sni[:-2,:,:])/2
# -- forward difference
# dsnidt[1:-1,:,:] = (sni[2:,:,:]-sni[1:-1,:,:])#/(dt)

# keep only regions of melt (dhdt < 0)
# dsnidt = np.minimum(np.zeros_like(dsnidt), dsnidt)
# nomelt=np.where(dsnidt > 0)
# dsnidt[nomelt] = 0

# sum total melt each month
fmelt=np.empty_like(tsurf)
for m in tqdm(range(tsurf.shape[0])):
    db=30*m # index of first day of this month
    de=30*(m+1) # index of first day of next month
    # fmelt[m,...]=-Lf*rhoh2o*np.nansum(dsnidt[db:de,...],axis=0)/(30*86400)
    fmelt[m,...]=-Lf*rhoh2o*np.nansum(dsnidt[db:de,...],axis=0)/(30*86400)

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
    
stend = file_stend.createVariable('snmelt', 'f4', ("time","lat","lon"))
stend.units = "W m**-2"
stend.long_name = "Snow melt flux (diagnosed)"
stend[:,:,:] = fmelt

file_stend.close()
