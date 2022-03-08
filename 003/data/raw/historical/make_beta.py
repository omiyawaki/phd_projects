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

# paths to ps and va files
path_ps = sys.argv[1]
path_va = sys.argv[2]
path_beta = sys.argv[3]

# open files
file_ps = Dataset(path_ps, 'r')
file_va = Dataset(path_va, 'r')

# read data
ps = np.squeeze(file_ps.variables['ps'][:]) # (lat x lon)
va = np.squeeze(file_va.variables['va'][:]) # (mon x lev x lat x lon)

plev = file_va.variables['plev'][:]
plev_half = 1/2 * (plev[1:] + plev[:-1])
plev_half = np.sort(np.append(plev_half, [pt, p0]))
plev_full = np.sort(np.concatenate((plev, plev_half)))
if plev[1]-plev[0]<0:
    plev_half = plev_half[::-1]
    plev_full = plev_full[::-1]
dplev = plev_half[1:] - plev_half[:-1]

# check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
lat2d = file_ps.variables['lat'][:] # (mon x lat x lon)
lat3d = file_va.variables['lat'][:] # (mon x lev x lat x lon)
if not np.array_equal(lat2d,lat3d):
    print('\nInterpolating ps to 3d lat grid...\n')
    filledps = ps.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filledps, axis=0)
    ps = f(filledlat3d)
    filledps = None; f = None;
    print('\nDone.\n')

# Create the beta parameter (for dealing with subsurface masking following Trenberth 1991)
b = np.empty( [1, va.shape[1], va.shape[2], va.shape[3]] )
for ilon in tqdm(range(ps.shape[1])):
    for ilat in range(ps.shape[0]):
        ps_local = ps[ilat, ilon]

        for ilev in range(len(plev)):
            ilevf = 2*ilev+1 # corresponding index in the plev_full array
            if plev_full[ilevf-1] < ps_local:
                b[0, ilev, ilat, ilon] = 1
            if plev_full[ilevf+1] > ps_local:
                b[0, ilev, ilat, ilon] = 0
            if (plev_full[ilevf-1] > ps_local) and (plev_full[ilevf+1] < ps_local):
                b[0, ilev, ilat, ilon] = (ps_local - plev_full[ilevf+1]) / (plev_full[ilevf-1] - plev_full[ilevf+1])

# save file as netCDF
file_beta = Dataset(path_beta, "w", format='NETCDF4_CLASSIC')

# copy attributes from va file
file_beta.setncatts(file_va.__dict__)

# copy dimensions from va file 
for name, dimension in file_va.dimensions.items():
    if any(name in s for s in ['time', 'time_bnds']): # copy time dimension from ps
        file_beta.createDimension(name, file_ps.dimensions[name].size)
    else:
        file_beta.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_va.variables.items():
    if any(name in s for s in ['va']): # don't copy va data
        continue
    elif any(name in s for s in ['time', 'time_bnds']): # copy time data from ps
        x = file_beta.createVariable(name, variable.datatype, variable.dimensions)
        file_beta[name].setncatts(file_ps[name].__dict__)
        file_beta[name][:] = file_ps[name][:]
    else: # copy all other axis data from va
        x = file_beta.createVariable(name, variable.datatype, variable.dimensions)
        file_beta[name].setncatts(file_va[name].__dict__)
        file_beta[name][:] = file_va[name][:]
    
beta = file_beta.createVariable('beta', 'f4', ("time","plev","lat","lon"))
beta.units = "unitless"
beta.long_name = "subsurface data weighting parameter for taking vertical integrals (see Trenberth 1991)"
beta[:,:] = b

file_beta.close()
