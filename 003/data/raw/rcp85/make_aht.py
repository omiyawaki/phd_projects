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
path_rlut = sys.argv[1]
path_rsut = sys.argv[2]
path_rsdt = sys.argv[3]
path_rlus = sys.argv[4]
path_rlds = sys.argv[5]
path_rsds = sys.argv[6]
path_rsus = sys.argv[7]
path_hfls = sys.argv[8]
path_hfss = sys.argv[9]
path_tend = sys.argv[10]
path_aht = sys.argv[11]

# open files
file_rlut = Dataset(path_rlut, 'r')
file_rsut = Dataset(path_rsut, 'r')
file_rsdt = Dataset(path_rsdt, 'r')
file_rlus = Dataset(path_rlus, 'r')
file_rlds = Dataset(path_rlds, 'r')
file_rsds = Dataset(path_rsds, 'r')
file_rsus = Dataset(path_rsus, 'r')
file_hfls = Dataset(path_hfls, 'r')
file_hfss = Dataset(path_hfss, 'r')
file_tend = Dataset(path_tend, 'r')

# read data
rlut = file_rlut.variables['rlut'][:] # (mon x lat x lon)
rsut = file_rsut.variables['rsut'][:] # (mon x lat x lon)
rsdt = file_rsdt.variables['rsdt'][:] # (mon x lat x lon)
rlus = file_rlus.variables['rlus'][:] # (mon x lat x lon)
rlds = file_rlds.variables['rlds'][:] # (mon x lat x lon)
rsds = file_rsds.variables['rsds'][:] # (mon x lat x lon)
rsus = file_rsus.variables['rsus'][:] # (mon x lat x lon)
hfls = file_hfls.variables['hfls'][:] # (mon x lat x lon)
hfss = file_hfss.variables['hfss'][:] # (mon x lat x lon)
tend = file_tend.variables['tend'][:] # (mon x lat x lon)
lat = file_tend.variables['lat'][:] # (lat)
rlat = np.radians(lat)
clat = np.cos(rlat)

# infer total energy flux divergence as the residual
fluxdiv = rsdt - rsut - rlut + rsus - rsds + rlus - rlds + hfls + hfss - tend

# zonal mean
fluxdiv_z = np.nanmean(fluxdiv, axis=2)

# subtract global mean
fluxdiv_g = np.trapz(clat*fluxdiv_z, rlat, axis=1) / np.trapz(clat, rlat)
fluxdiv_z = fluxdiv_z - np.transpose(np.tile(fluxdiv_g, [len(lat), 1]))

# meridionally integrate
aht = 2*np.pi*a**2 * integrate.cumtrapz(clat*fluxdiv_z, rlat, axis=1, initial=0)
if lat[1]-lat[0]<0:
    aht = -aht

# save file as netCDF
file_aht = Dataset(path_aht, "w", format='NETCDF4_CLASSIC')

# copy attributes from rlut file
file_aht.setncatts(file_rlut.__dict__)

# copy dimensions from rlut file
for name, dimension in file_rlut.dimensions.items():
    if any(name in s for s in ['plev']):
        continue
    file_aht.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_rlut.variables.items():
    if any(name in s for s in ['rlut' 'plev']):
        continue
    
    x = file_aht.createVariable(name, variable.datatype, variable.dimensions)
    file_aht[name].setncatts(file_rlut[name].__dict__)
    file_aht[name][:] = file_rlut[name][:]
    
vE = file_aht.createVariable('vE', 'f4', ("time","lat"))
vE.units = "W"
vE.long_name = "vertically integrated total energy flux transport"
vE[:,:] = aht

file_aht.close()
