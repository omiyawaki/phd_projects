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
path_vmte = sys.argv[1]
path_divte = sys.argv[2]

# open files
file_vmte = Dataset(path_vmte, 'r')

# read data
vm_te = file_vmte.variables['vmte'][:] # (mon x lat)
lat = file_vmte.variables['lat'][:] # (mon x lat x lon)

# # for datasets that fill data below surface as missing data, fill with nans
vm_te = vm_te.filled(fill_value=np.nan)
lat = lat.filled(fill_value=np.nan)

# take meridional derivative
rlat = np.radians(lat)
clat = np.cos(rlat)

div_vm_te = np.empty_like(vm_te)

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

div_vm_te_half = 1/(2*np.pi*a**2 * 1/2*(clat[1:]+clat[:-1])**2 ) * (clat[1:]*vm_te[...,1:] - clat[:-1]*vm_te[...,:-1]) / (rlat[1:] - rlat[:-1])
f_div_vm_te = interpolate.interp1d(rlat_half, div_vm_te_half, bounds_error=False)
div_vm_te = f_div_vm_te(rlat)

# save file as netCDF
file_te = Dataset(path_divte, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_te.setncatts(file_vmte.__dict__)

# copy dimensions from vmte file
for name, dimension in file_vmte.dimensions.items():
    if any(name in s for s in ['plev']):
        continue
    file_te.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_vmte.variables.items():
    if any(name in s for s in ['vmte' 'plev']):
        continue
    
    x = file_te.createVariable(name, variable.datatype, variable.dimensions)
    file_te[name].setncatts(file_vmte[name].__dict__)
    file_te[name][:] = file_vmte[name][:]
    
te = file_te.createVariable('divte', 'f4', ("time","lat"))
te.units = "W m**-2"
te.long_name = "vertically integrated moist static energy flux divergence due to transient eddies"
te[:,:] = div_vm_te

file_te.close()
