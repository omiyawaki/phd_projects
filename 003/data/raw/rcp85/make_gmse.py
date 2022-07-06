import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts/misc')
from tvregdiff import TVRegDiff
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP
cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps and mse files
path_mse = sys.argv[1]
path_gmse = sys.argv[2]

# open files
file_mse = Dataset(path_mse, 'r')

# read data
mse = file_mse.variables['mse'][:] # (mon x lat x lon)
lat = file_mse.variables['lat'][:] # (lat)

# for datasets that fill data below surface as missing data, fill with nans
mse = mse.filled(fill_value=np.nan)
lat = lat.filled(fill_value=np.nan)

# take meridional derivative
rlat = np.radians(lat)

gmse = np.empty_like(mse)

######################################
### DERIVATIVE OF AN INTERPOLATION
######################################

#n = len(rlat)
#k = 3 # cubic spline
## s = n - np.sqrt(2*n) # smoothing factor
#s=n
#prefactor = 1/(2*np.pi*a**2*clat)
#for itime in tqdm(range(vm_mmc.shape[0])):
# spl_vE = interpolate.UnivariateSpline(rlat, vE[itime,:], k=k, s=s).derivative(1)
# spl_vm_mmc = interpolate.UnivariateSpline(rlat, vm_mmc[itime,:], k=k, s=s).derivative(1)
# spl_vm_se = interpolate.UnivariateSpline(rlat, vm_se[itime,:], k=k, s=s).derivative(1)
# spl_vm_te = interpolate.UnivariateSpline(rlat, vm_te[itime,:], k=k, s=s).derivative(1)

# div_vE[itime,:] = prefactor*spl_vE(rlat)
# div_vm_mmc[itime,:] = prefactor*spl_vm_mmc(rlat)
# div_vm_se[itime,:] = prefactor*spl_vm_se(rlat)
# div_vm_te[itime,:] = prefactor*spl_vm_te(rlat)

# # div_vm_mmc[itime,:] = prefactor*TVRegDiff(vm_mmc[itime,:], 1, 1e-1, dx=np.mean(rlat[1]-rlat[0]), plotflag=True)

#####################################
## CENTERED DIFFERENCE
#####################################
rlat_half = 1/2*(rlat[1:]+rlat[:-1])

gmse_half = (mse[:,1:,:] - mse[:,:-1,:]) / (rlat[None,1:,None] - rlat[None,:-1,None])
f_gmse = interpolate.interp1d(rlat_half, gmse_half, axis=1, bounds_error=False)
gmse = f_gmse(rlat)

#####################################
## FORWARD AND BACKWARD DIFFERENCE at ENDPOINTS
#####################################
gmse[:,0,:] = (mse[:,1,:] - mse[:,0,:]) / (rlat[1] - rlat[0])
gmse[:,1,:] = (mse[:,1,:] - mse[:,2,:]) / (rlat[-1] - rlat[-2])

# save file as netCDF
file_gmse = Dataset(path_gmse, "w", format='NETCDF4_CLASSIC')

# copy attributes
file_gmse.setncatts(file_mse.__dict__)

# copy dimensions
for name, dimension in file_mse.dimensions.items():
    file_gmse.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables
for name, variable in file_mse.variables.items():
    if any(name in s for s in ['mse']):
        continue
    
    x = file_gmse.createVariable(name, variable.datatype, variable.dimensions)
    file_gmse[name].setncatts(file_mse[name].__dict__)
    file_gmse[name][:] = file_mse[name][:]
    
var_gmse = file_gmse.createVariable('gmse', 'f4', ("time","lat","lon"))
# var_gmse.units = "J kg**-1 m**-1"
var_gmse.units = "J kg**-1"
var_gmse.long_name = "92500 Pa meridional MSE gradient"
var_gmse[:] = gmse
file_gmse.close()
