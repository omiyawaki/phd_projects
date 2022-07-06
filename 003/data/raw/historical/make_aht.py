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

# infer energy flux divergence plus storage (dps) as the residual
dps = rsdt - rsut - rlut + rsus - rsds + rlus - rlds + hfls + hfss

# check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
lat2d = file_rlut.variables['lat'][:] 
lat3d = file_tend.variables['lat'][:] 
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

# isolate flux divergence from dps by subtracting the MSE tendency
fluxdiv = dps - tend

# take zonal mean
fluxdiv_z = np.transpose( np.tile( np.nanmean(fluxdiv, axis=2), [fluxdiv.shape[2],1,1] ), [1,2,0] )

# subtract global mean
fluxdiv_g = np.trapz(clat[None,:,None]*fluxdiv_z, rlat, axis=1) / np.trapz(clat, rlat)
fluxdiv_z = fluxdiv_z - np.transpose(np.tile(fluxdiv_g, [len(lat), 1, 1]), [1,0,2])

# meridionally integrate
aht = 2*np.pi*a**2 * integrate.cumtrapz(clat[None,:,None]*fluxdiv_z, rlat, axis=1, initial=0)
if lat[1]-lat[0]<0:
    aht = -aht

# save file as netCDF
file_aht = Dataset(path_aht, "w", format='NETCDF4_CLASSIC')

# copy attributes from tend file
file_aht.setncatts(file_tend.__dict__)

# copy dimensions from tend file
for name, dimension in file_tend.dimensions.items():
    file_aht.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from tend file
for name, variable in file_tend.variables.items():
    if any(name in s for s in ['tend']):
        continue
    
    x = file_aht.createVariable(name, variable.datatype, variable.dimensions)
    file_aht[name].setncatts(file_tend[name].__dict__)
    file_aht[name][:] = file_tend[name][:]
    
vE = file_aht.createVariable('aht', 'f4', ("time","lat","lon"))
vE.units = "W"
vE.long_name = "vertically integrated total energy flux transport"
vE[:] = aht

file_aht.close()
