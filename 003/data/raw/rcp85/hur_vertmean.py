import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;
p0 = 1100.e2 # bottom integral bound [hPa]
pt = 0. # top integral bound [hPa]

# Paths to ps and hur files
path_ps = sys.argv[1]
path_hur = sys.argv[2]
path_beta = sys.argv[3]
path_vhur = sys.argv[4]

# Open files
file_ps = Dataset(path_ps, 'r')
file_hur = Dataset(path_hur, 'r')
file_beta = Dataset(path_beta, 'r')

# Read data
ps = np.squeeze(file_ps.variables['ps'][:]) # (lat x lon)
hur = file_hur.variables['hur'][:] # (mon x lev x lat x lon)
beta = file_beta.variables['beta'][:] # (1 x lev x lat x lon)

# For datasets that fill data below surface as missing data, set as nans
hur = hur.filled(fill_value=np.nan)

# Grid data
lat2d = file_ps.variables['lat'][:] # (mon x lat x lon)
lat3d = file_hur.variables['lat'][:] # (mon x lev x lat x lon)
rlat = np.radians(lat3d)
clat = np.cos(rlat)

lon = file_hur.variables['lon'][:]
dlon = np.pi / 180 * (lon[1] - lon[0])

plev = file_hur.variables['plev'][:]
plev_half = 1/2 * (plev[1:] + plev[:-1])
plev_half = np.sort(np.append(plev_half, [pt, p0]))
plev_full = np.sort(np.concatenate((plev, plev_half)))
if plev[1]-plev[0]<0:
    plev_half = plev_half[::-1]
    plev_full = plev_full[::-1]
dplev = plev_half[1:] - plev_half[:-1]

# Check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
if not np.array_equal(lat2d,lat3d):
    print('\nInterpolating ps to 3d lat grid...\n')
    filledps = ps.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filledps, axis=0)
    ps = f(filledlat3d)
    filledps = None; f = None;
    print('\nDone.\n')

# Zonal surface pressure
ps_z = np.nanmean(ps, axis=1)
ps_tile = np.tile(ps, [hur.shape[0], hur.shape[1], 1, 1])

# Mask out data below surface using long-term annual mean surface pressure
pa = np.transpose(np.tile(plev, [hur.shape[0], hur.shape[2], hur.shape[3], 1]), [0,3,1,2])
subsurf = pa > ps_tile
hur[subsurf] = np.nan
beta[subsurf[[0],...]] = np.nan
pa = None

# compute vertical mean hur
dplev_ext = dplev[None,:,None,None] # extend dplev dimensions (mon x lev x lat x lon)
denom_v = np.nansum(beta*dplev_ext, axis=1, keepdims=True)
vhur = np.nansum(beta*hur*dplev_ext, axis=1, keepdims=True) / denom_v
print(vhur.shape)
vhur = np.squeeze(vhur)
print(vhur.shape)

# save file as netCDF
file_vhur = Dataset(path_vhur, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_vhur.setncatts(file_hur.__dict__)

# copy dimensions from hur file
for name, dimension in file_hur.dimensions.items():
    if any(name in s for s in ['plev']):
        continue

    file_vhur.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from hur file
for name, variable in file_hur.variables.items():
    if any(name in s for s in ['hur', 'plev', 'plev_bnds']):
        continue
    
    x = file_vhur.createVariable(name, variable.datatype, variable.dimensions)
    file_vhur[name].setncatts(file_hur[name].__dict__)
    file_vhur[name][:] = file_hur[name][:]

mean_hur = file_vhur.createVariable('vhur', 'f4', ("time","lat","lon"))
mean_hur.units = "%"
mean_hur.long_name = "Vertical-mean relative humidity"
mean_hur[:,:,:] = vhur

file_vhur.close()
