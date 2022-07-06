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
path_pp0 = sys.argv[1]
path_psic = sys.argv[2]
path_beta = sys.argv[3]
path_ubeta = sys.argv[4]

# open files
file_pp0 = Dataset(path_pp0, 'r')
file_psic = Dataset(path_psic, 'r')
file_beta = Dataset(path_beta, 'r')

# read data
pp0 = np.squeeze(file_pp0.variables['pp0'][:]) # (mon x lat)
psic = np.squeeze(file_psic.variables['psic'][:]) # (mon x lev x lat)
beta = np.squeeze(file_beta.variables['beta'][:]) # (lev x lat x lon)

plev = file_psic.variables['plev'][:]
plev_half = 1/2 * (plev[1:] + plev[:-1])
plev_half = np.sort(np.append(plev_half, [pt, p0]))
plev_full = np.sort(np.concatenate((plev, plev_half)))
if plev[1]-plev[0]<0:
    plev_half = plev_half[::-1]
    plev_full = plev_full[::-1]
dplev = plev_half[1:] - plev_half[:-1]

# Create the upper beta parameter 

b = np.empty_like(psic)
for itime in tqdm(range(psic.shape[0])):
    for ilat in range(psic.shape[2]):
        pp0_local = pp0[itime, ilat]

        for ilev in range(len(plev)):
            ilevf = 2*ilev+1 # corresponding index in the plev_full array
            if plev_full[ilevf-1] < pp0_local: # above pp0
                b[itime, ilev, ilat] = 1
            if plev_full[ilevf+1] > pp0_local: # below pp0
                b[itime, ilev, ilat] = 0
            if (plev_full[ilevf-1] > pp0_local) and (plev_full[ilevf+1] < pp0_local):
                b[itime, ilev, ilat] = (pp0_local - plev_full[ilevf+1]) / (plev_full[ilevf-1] - plev_full[ilevf+1])

ub = np.transpose( np.tile(b, [beta.shape[2],1,1,1]), [1,2,3,0] )

print(ub[0,:,0,0])

# save file as netCDF
file_ubeta = Dataset(path_ubeta, "w", format='NETCDF4_CLASSIC')

# copy attributes from psic file
file_ubeta.setncatts(file_psic.__dict__)

# copy dimensions from psic file 
for name, dimension in file_psic.dimensions.items():
    file_ubeta.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

# copy lon dimension from beta file 
for name, dimension in file_beta.dimensions.items():
    if name in ['lon', 'lon_bnds']:
        file_ubeta.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from psic file
for name, variable in file_psic.variables.items():
    if any(name in s for s in ['psic']): # don't copy psic data
        continue

    x = file_ubeta.createVariable(name, variable.datatype, variable.dimensions)
    file_ubeta[name].setncatts(file_psic[name].__dict__)
    file_ubeta[name][:] = file_psic[name][:]

# copy lon from beta
for name, variable in file_beta.variables.items():
    if name in ['lon', 'lon_bnds']:
        x = file_ubeta.createVariable(name, variable.datatype, variable.dimensions)
        file_ubeta[name].setncatts(file_beta[name].__dict__)
        file_ubeta[name][:] = file_beta[name][:]
    
ubeta = file_ubeta.createVariable('ubeta', 'f4', ("time","plev","lat","lon"))
ubeta.units = "unitless"
ubeta.long_name = "upper cell mask"
ubeta[:] = ub

file_ubeta.close()
