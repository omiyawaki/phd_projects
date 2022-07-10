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
path_full_r1 = sys.argv[1]
path_full_ra = sys.argv[2]
path_full_dyn = sys.argv[3]
path_ref_r1 = sys.argv[4]
path_ref_ra = sys.argv[5]
path_ref_dyn = sys.argv[6]
path_diff_r1 = sys.argv[7]
path_diff_ra = sys.argv[8]
path_diff_dyn = sys.argv[9]
path_diff_res = sys.argv[10]
path_full_stf = sys.argv[11]
path_ref_stf = sys.argv[12]
path_diff_stf = sys.argv[13]

# open files
file_full_r1 = Dataset(path_full_r1, 'r')
file_full_ra = Dataset(path_full_ra, 'r')
file_full_stf = Dataset(path_full_stf, 'r')
file_full_dyn = Dataset(path_full_dyn, 'r')
file_ref_r1 = Dataset(path_ref_r1, 'r')
file_ref_ra = Dataset(path_ref_ra, 'r')
file_ref_stf = Dataset(path_ref_stf, 'r')
file_ref_dyn = Dataset(path_ref_dyn, 'r')

# read data
r1 = file_full_r1.variables['r1'][:] # (mon x lat)
ra = file_full_ra.variables['ra'][:] # (mon x lat)
stf = file_full_stf.variables['stf'][:] # (mon x lat)
dyn = file_full_dyn.variables['stgadv'][:] # (mon x lat)
r1_r = file_ref_r1.variables['r1'][:] # (1 x lat)
ra_r = file_ref_ra.variables['ra'][:] # (1 x lat)
stf_r = file_ref_stf.variables['stf'][:] # (1 x lat)
dyn_r = file_ref_dyn.variables['stgadv'][:] # (1 x lat)

# take zonal means
ra = np.nanmean(ra, axis=2)
stf = np.nanmean(stf, axis=2)
dyn = np.nanmean(dyn, axis=2)
ra_r = np.nanmean(ra_r, axis=2)
stf_r = np.nanmean(stf_r, axis=2)
dyn_r = np.nanmean(dyn_r, axis=2)

# take the differences
dr1 = r1 - r1_r
dra = ra - ra_r
dstf = stf - stf_r
ddyn = dyn - dyn_r

# compute components
dcdyn = ddyn / ra_r
dcra = - dyn_r / (ra_r**2) * dra
dcres = dr1 - dcdyn - dcra
dcstf = (dstf) / ra_r

# save file as netCDF
file_diff_r1 = Dataset(path_diff_r1, "w", format='NETCDF4_CLASSIC')
file_diff_ra = Dataset(path_diff_ra, "w", format='NETCDF4_CLASSIC')
file_diff_stf = Dataset(path_diff_stf, "w", format='NETCDF4_CLASSIC')
file_diff_dyn = Dataset(path_diff_dyn, "w", format='NETCDF4_CLASSIC')
file_diff_res = Dataset(path_diff_res, "w", format='NETCDF4_CLASSIC')

# copy attributes from r1
file_diff_r1.setncatts(file_full_r1.__dict__)
file_diff_ra.setncatts(file_full_r1.__dict__)
file_diff_stf.setncatts(file_full_r1.__dict__)
file_diff_dyn.setncatts(file_full_r1.__dict__)
file_diff_res.setncatts(file_full_r1.__dict__)

# copy dimensions from ps file
for name, dimension in file_full_r1.dimensions.items():
    file_diff_r1.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_diff_ra.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_diff_stf.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_diff_dyn.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_diff_res.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_full_r1.variables.items():
    if any(name in s for s in ['r1']):
        continue
    
    x = file_diff_r1.createVariable(name, variable.datatype, variable.dimensions)
    file_diff_r1[name].setncatts(file_full_r1[name].__dict__)
    file_diff_r1[name][:] = file_full_r1[name][:]

    x = file_diff_ra.createVariable(name, variable.datatype, variable.dimensions)
    file_diff_ra[name].setncatts(file_full_r1[name].__dict__)
    file_diff_ra[name][:] = file_full_r1[name][:]

    x = file_diff_stf.createVariable(name, variable.datatype, variable.dimensions)
    file_diff_stf[name].setncatts(file_full_r1[name].__dict__)
    file_diff_stf[name][:] = file_full_r1[name][:]

    x = file_diff_dyn.createVariable(name, variable.datatype, variable.dimensions)
    file_diff_dyn[name].setncatts(file_full_r1[name].__dict__)
    file_diff_dyn[name][:] = file_full_r1[name][:]

    x = file_diff_res.createVariable(name, variable.datatype, variable.dimensions)
    file_diff_res[name].setncatts(file_full_r1[name].__dict__)
    file_diff_res[name][:] = file_full_r1[name][:]

diff_r1 = file_diff_r1.createVariable('dr1', 'f4', ("time","lat"))
diff_r1.units = "unitless"
diff_r1.long_name = "r1 deviation from reference state (historical)"
diff_r1[:,:] = dr1
file_diff_r1.close()

diff_ra = file_diff_ra.createVariable('dcra', 'f4', ("time","lat"))
diff_ra.units = "unitless"
diff_ra.long_name = "contribution of r1 deviation from reference state (historical) due to radiative cooling"
diff_ra[:,:] = dcra
file_diff_ra.close()

diff_stf = file_diff_stf.createVariable('dcstf', 'f4', ("time","lat"))
diff_stf.units = "unitless"
diff_stf.long_name = "contribution of r1 deviation from reference state (historical) due to surface turbulent fluxes"
diff_stf[:,:] = dcstf
file_diff_stf.close()

diff_dyn = file_diff_dyn.createVariable('dcdyn', 'f4', ("time","lat"))
diff_dyn.units = "unitless"
diff_dyn.long_name = "contribution of r1 deviation from reference state (historical) due to storage plus advection"
diff_dyn[:,:] = dcdyn
file_diff_dyn.close()

diff_res = file_diff_res.createVariable('dcres', 'f4', ("time","lat"))
diff_res.units = "unitless"
diff_res.long_name = "contribution of r1 deviation from reference state (historical) due to nonlinear interactions"
diff_res[:,:] = dcres
file_diff_res.close()

