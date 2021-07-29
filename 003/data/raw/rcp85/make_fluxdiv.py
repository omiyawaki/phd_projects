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
path_aht = sys.argv[1]
path_vmmmc = sys.argv[2]
path_vmse = sys.argv[3]
path_vmte = sys.argv[4]
path_divaht = sys.argv[5]
path_divmmc = sys.argv[6]
path_divse = sys.argv[7]
path_divte = sys.argv[8]

# open files
file_aht = Dataset(path_aht, 'r')
file_vmmmc = Dataset(path_vmmmc, 'r')
file_vmse = Dataset(path_vmse, 'r')
file_vmte = Dataset(path_vmte, 'r')

# read data
vE = file_aht.variables['vE'][:] # (mon x lat)
vm_mmc = file_vmmmc.variables['vmmmc'][:] # (mon x lat)
vm_se = file_vmse.variables['vmse'][:] # (mon x lat)
vm_te = file_vmte.variables['vmte'][:] # (mon x lat)
lat = file_vmmmc.variables['lat'][:] # (mon x lat x lon)

# # for datasets that fill data below surface as missing data, fill with nans
vE = vE.filled(fill_value=np.nan)
vm_mmc = vm_mmc.filled(fill_value=np.nan)
vm_se = vm_se.filled(fill_value=np.nan)
vm_te = vm_te.filled(fill_value=np.nan)
lat = lat.filled(fill_value=np.nan)

# take meridional derivative
rlat = np.radians(lat)
clat = np.cos(rlat)

div_vE = np.empty_like(vE)
div_vm_mmc = np.empty_like(vm_mmc)
div_vm_se = np.empty_like(vm_se)
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

div_vE_half = 1/(2*np.pi*a**2 * 1/2*(clat[1:]+clat[:-1]) ) * (vE[...,1:] - vE[...,:-1]) / (rlat[1:] - rlat[:-1])
f_div_vE = interpolate.interp1d(rlat_half, div_vE_half, bounds_error=False)
div_vE = f_div_vE(rlat)

div_vm_mmc_half = 1/(2*np.pi*a**2 * 1/2*(clat[1:]+clat[:-1])**2 ) * (clat[1:]*vm_mmc[...,1:] - clat[:-1]*vm_mmc[...,:-1]) / (rlat[1:] - rlat[:-1])
f_div_vm_mmc = interpolate.interp1d(rlat_half, div_vm_mmc_half, bounds_error=False)
div_vm_mmc = f_div_vm_mmc(rlat)

div_vm_se_half = 1/(2*np.pi*a**2 * 1/2*(clat[1:]+clat[:-1])**2 ) * (clat[1:]*vm_se[...,1:] - clat[:-1]*vm_se[...,:-1]) / (rlat[1:] - rlat[:-1])
f_div_vm_se = interpolate.interp1d(rlat_half, div_vm_se_half, bounds_error=False)
div_vm_se = f_div_vm_se(rlat)

div_vm_te_half = 1/(2*np.pi*a**2 * 1/2*(clat[1:]+clat[:-1])**2 ) * (clat[1:]*vm_te[...,1:] - clat[:-1]*vm_te[...,:-1]) / (rlat[1:] - rlat[:-1])
f_div_vm_te = interpolate.interp1d(rlat_half, div_vm_te_half, bounds_error=False)
div_vm_te = f_div_vm_te(rlat)

# save file as netCDF
file_tot = Dataset(path_divaht, "w", format='NETCDF4_CLASSIC')
file_mmc = Dataset(path_divmmc, "w", format='NETCDF4_CLASSIC')
file_se = Dataset(path_divse, "w", format='NETCDF4_CLASSIC')
file_te = Dataset(path_divte, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_tot.setncatts(file_aht.__dict__)
file_mmc.setncatts(file_vmmmc.__dict__)
file_se.setncatts(file_vmse.__dict__)
file_te.setncatts(file_vmte.__dict__)

# copy dimensions from vmmmc file
for name, dimension in file_vmmmc.dimensions.items():
    if any(name in s for s in ['plev']):
        continue
    file_tot.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_mmc.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_se.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_te.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_vmmmc.variables.items():
    if any(name in s for s in ['vmmmc' 'plev']):
        continue
    
    x = file_tot.createVariable(name, variable.datatype, variable.dimensions)
    file_tot[name].setncatts(file_vmmmc[name].__dict__)
    file_tot[name][:] = file_vmmmc[name][:]

    x = file_mmc.createVariable(name, variable.datatype, variable.dimensions)
    file_mmc[name].setncatts(file_vmmmc[name].__dict__)
    file_mmc[name][:] = file_vmmmc[name][:]

    x = file_se.createVariable(name, variable.datatype, variable.dimensions)
    file_se[name].setncatts(file_vmmmc[name].__dict__)
    file_se[name][:] = file_vmmmc[name][:]

    x = file_te.createVariable(name, variable.datatype, variable.dimensions)
    file_te[name].setncatts(file_vmmmc[name].__dict__)
    file_te[name][:] = file_vmmmc[name][:]
    
tot = file_tot.createVariable('divaht', 'f4', ("time","lat"))
tot.units = "W m**-2"
tot.long_name = "vertically integrated total energy flux divergence"
tot[:,:] = div_vE

mmc = file_mmc.createVariable('divmmc', 'f4', ("time","lat"))
mmc.units = "W m**-2"
mmc.long_name = "vertically integrated moist static energy flux divergence due to mean meridional circulation"
mmc[:,:] = div_vm_mmc

se = file_se.createVariable('divse', 'f4', ("time","lat"))
se.units = "W m**-2"
se.long_name = "vertically integrated moist static energy flux divergence due to stationary eddies"
se[:,:] = div_vm_se

te = file_te.createVariable('divte', 'f4', ("time","lat"))
te.units = "W m**-2"
te.long_name = "vertically integrated moist static energy flux divergence due to transient eddies"
te[:,:] = div_vm_te

file_tot.close()
file_mmc.close()
file_se.close()
file_te.close()
