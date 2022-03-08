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
path_evap = sys.argv[1]
path_prec = sys.argv[2]
path_qtend = sys.argv[3]
path_qaht = sys.argv[4]

# open files
file_evap = Dataset(path_evap, 'r')
file_prec = Dataset(path_prec, 'r')
file_qtend = Dataset(path_qtend, 'r')

# read data
evap = file_evap.variables['evspsbl'][:] # (mon x lat x lon)
prec = file_prec.variables['pr'][:] # (mon x lat x lon)
qtend = file_qtend.variables['qtend'][:] # (mon x lat x lon)

# infer humidity divergence plus storage (dps) as the residual
dps = evap - prec

# check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
lat2d = file_evap.variables['lat'][:] 
lat3d = file_qtend.variables['lat'][:] 
if not np.array_equal(lat2d,lat3d):
    filleddps = dps.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filleddps, axis=1)
    dps = f(filledlat3d)
    filleddps = None; f = None;

lat=lat3d
rlat = np.radians(lat)
clat = np.cos(rlat)

# isolate flux divergence from dps by subtracting the humidity tendency
fluxdiv = dps - qtend

# take zonal mean
fluxdiv_z = np.nanmean(fluxdiv, axis=2)

# subtract global mean
fluxdiv_g = np.trapz(clat*fluxdiv_z, rlat, axis=1) / np.trapz(clat, rlat)
fluxdiv_z = fluxdiv_z - np.transpose(np.tile(fluxdiv_g, [len(lat), 1]))

# meridionally integrate
qaht = 2*np.pi*a**2*L * integrate.cumtrapz(clat*fluxdiv_z, rlat, axis=1, initial=0)
if lat[1]-lat[0]<0:
    qaht = -qaht

# save file as netCDF
file_qaht = Dataset(path_qaht, "w", format='NETCDF4_CLASSIC')

# copy attributes from qtend file
file_qaht.setncatts(file_qtend.__dict__)

# copy dimensions from qtend file
for name, dimension in file_qtend.dimensions.items():
    file_qaht.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from qtend file
for name, variable in file_qtend.variables.items():
    if any(name in s for s in ['qtend']):
        continue
    
    x = file_qaht.createVariable(name, variable.datatype, variable.dimensions)
    file_qaht[name].setncatts(file_qtend[name].__dict__)
    file_qaht[name][:] = file_qtend[name][:]
    
vqE = file_qaht.createVariable('qaht', 'f4', ("time","lat"))
vqE.units = "W"
vqE.long_name = "vertically integrated moist energy flux transport"
vqE[:,:] = qaht

file_qaht.close()
