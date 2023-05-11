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
rhow = 1000 # kg m**-3 density of water
rhos = 300 # kg m**-3 density of snow
ki = 2.1656 # ice thermal conductivity [W/K/m]
ks = 0.31 # snow thermal conductivity [W/K/m]

# paths to ps and siced files
path_seaice = sys.argv[1]
path_tsi = sys.argv[2]
path_siced = sys.argv[3]
path_sni = sys.argv[4]
path_cond = sys.argv[5]
path_aprs = sys.argv[6]
path_seaice = sys.argv[7]
path_tmelt = sys.argv[8]
path_smelt = sys.argv[9]
path_snmelt = sys.argv[10]

# open files
file_seaice = Dataset(path_seaice, 'r')
file_tsi = Dataset(path_tsi, 'r')
file_siced = Dataset(path_siced, 'r')
file_sni = Dataset(path_sni, 'r')
file_cond = Dataset(path_cond, 'r')
file_aprs = Dataset(path_aprs, 'r')
file_seaice = Dataset(path_seaice, 'r')

# read data
seaice = file_seaice.variables['seaice'][:] # (time x lev x lat x lon)
tsi = file_tsi.variables['tsi'][:] # (time x lev x lat x lon)
siced = file_siced.variables['siced'][:] # (time x lev x lat x lon)
sni = file_sni.variables['sni'][:] # (time x lev x lat x lon)
cond = file_cond.variables['cond'][:] # (time x lev x lat x lon)
aprs = file_aprs.variables['aprs'][:] # (time x lev x lat x lon)
seaice = file_seaice.variables['seaice'][:] # (time x lev x lat x lon)

if seaice.shape[0]*30 != siced.shape[0]:
    error('This script uses daily sea ice thickness data. Please check that sea ice contains daily data.')

# take time tmeltency
dt=30*86400
dsiceddt = np.empty(siced.shape) 
dsiceddt[0] = (siced[1,:,:]-siced[0,:,:])/dt
dsiceddt[-1] = (siced[-1,:,:]-siced[-2,:,:])/dt

dsnidt = np.empty(sni.shape) 
dsnidt[0] = (sni[1,:,:]-sni[0,:,:])/dt
dsnidt[-1] = (sni[-1,:,:]-sni[-2,:,:])/dt

# # -- forward difference
# dsiceddt[1:-1,:,:] = (siced[2:,:,:]-siced[1:-1,:,:])#/(dt)
# dsnidt[1:-1,:,:] = (sni[2:,:,:]-sni[1:-1,:,:])#/(dt)

# # -- backward difference
# dsiceddt[1:-1,:,:] = (-siced[:-2,:,:]+siced[1:-1,:,:])#/(dt)
# dsnidt[1:-1,:,:] = (-sni[:-2,:,:]+sni[1:-1,:,:])#/(dt)

# -- centered difference
dsiceddt[1:-1,:,:] = (siced[2:,:,:]-siced[:-2,:,:])/(2*dt)
dsnidt[1:-1,:,:] = (sni[2:,:,:]-sni[:-2,:,:])/(2*dt)

# sea ice melt matters for surface budget only when surface is bare ice (sni == 0)
sc = np.where(sni > 0) # where there is a snow layer
dsiceddt[sc] = 0

# keep only regions of sea ice melt (cond >= 0)
# frz = np.where(cond < 0) # where sea ice is growing
# dsiceddt[frz] = 0
# frz = np.where(np.logical_and(cond == 0, dsiceddt>0)) # where sea ice is growing
# dsiceddt[frz] = 0

dsiceddt = np.minimum(np.zeros_like(dsiceddt), dsiceddt)
dsnidt = np.minimum(np.zeros_like(dsnidt), dsnidt)

# sum total melt each month
fmelt=np.empty_like(seaice)
seaicemelt=np.empty_like(seaice)
snowmelt=np.empty_like(seaice)
for m in tqdm(range(seaice.shape[0])):
    db=30*m # index of first day of this month
    de=30*(m+1) # index of first day of next month
    seaicemelt[m,...]=-Lf*rhoi*np.nansum(dsiceddt[db:de,...],axis=0) # bare sea ice melt
    snowmelt[m,...]=-Lf*rhow*np.nansum(dsnidt[db:de,...],axis=0) # snow melt
    fmelt[m,...]=seaicemelt[m,...]+snowmelt[m,...] # total melt
    # fmelt[m,...]=fmelt[m,...]+Lf*np.nansum(aprs[db:de,...],axis=0)/30
    # fmelt[m,...]=fmelt[m,...]-np.nansum(cond[db:de,...],axis=0)/30

# save file as netCDF
file_tmelt = Dataset(path_tmelt, "w", format='NETCDF4_CLASSIC')
file_smelt = Dataset(path_smelt, "w", format='NETCDF4_CLASSIC')
file_snmelt = Dataset(path_snmelt, "w", format='NETCDF4_CLASSIC')

# copy attributes from seaice file
file_tmelt.setncatts(file_seaice.__dict__)
file_smelt.setncatts(file_seaice.__dict__)
file_snmelt.setncatts(file_seaice.__dict__)

# copy dimensions from seaice file
for name, dimension in file_seaice.dimensions.items():
    file_tmelt.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_smelt.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_snmelt.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_seaice.variables.items():
    if any(name in s for s in ['seaice']):
        continue
    
    x = file_tmelt.createVariable(name, variable.datatype, variable.dimensions)
    file_tmelt[name].setncatts(file_seaice[name].__dict__)
    file_tmelt[name][:] = file_seaice[name][:]
    
    x = file_smelt.createVariable(name, variable.datatype, variable.dimensions)
    file_smelt[name].setncatts(file_seaice[name].__dict__)
    file_smelt[name][:] = file_seaice[name][:]
    
    x = file_snmelt.createVariable(name, variable.datatype, variable.dimensions)
    file_snmelt[name].setncatts(file_seaice[name].__dict__)
    file_snmelt[name][:] = file_seaice[name][:]
    
tmelt = file_tmelt.createVariable('tmelt', 'f4', ("time","lat","lon"))
tmelt.units = "W m**-2"
tmelt.long_name = "Total (sea ice and snow) melt flux (diagnosed)"
tmelt[:,:,:] = fmelt
file_tmelt.close()

smelt = file_smelt.createVariable('smelt', 'f4', ("time","lat","lon"))
smelt.units = "W m**-2"
smelt.long_name = "Bare sea ice melt flux (diagnosed)"
smelt[:,:,:] = seaicemelt
file_smelt.close()

snmelt = file_snmelt.createVariable('snmelt', 'f4', ("time","lat","lon"))
snmelt.units = "W m**-2"
snmelt.long_name = "Snow melt flux (diagnosed)"
snmelt[:,:,:] = snowmelt
file_snmelt.close()
