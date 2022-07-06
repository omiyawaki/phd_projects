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
path_full_vmte = sys.argv[1]
path_full_dvp = sys.argv[2]
path_full_gmse = sys.argv[3]
path_ref_vmte = sys.argv[4]
path_ref_dvp = sys.argv[5]
path_ref_gmse = sys.argv[6]
path_diff_vmte = sys.argv[7]
path_diff_dvp = sys.argv[8]
path_diff_gmse = sys.argv[9]
path_diff_res = sys.argv[10]
clat = float(sys.argv[11])

# open files
file_full_vmte = Dataset(path_full_vmte, 'r')
file_full_dvp = Dataset(path_full_dvp, 'r')
file_full_gmse = Dataset(path_full_gmse, 'r')
file_ref_vmte = Dataset(path_ref_vmte, 'r')
file_ref_dvp = Dataset(path_ref_dvp, 'r')
file_ref_gmse = Dataset(path_ref_gmse, 'r')

# read data
vmte = file_full_vmte.variables['vmte_sm'][:] # (mon x lat)
dvp = file_full_dvp.variables['dvp'][:] # (mon x lat)
gmse = file_full_gmse.variables['gmse_sm'][:] # (mon x lat)
vmte_r = file_ref_vmte.variables['vmte'][:] # (1 x lat)
dvp_r = file_ref_dvp.variables['dvp'][:] # (1 x lat)
gmse_r = file_ref_gmse.variables['gmse'][:] # (1 x lat)

# take the differences
dvmte = vmte - vmte_r
ddvp = dvp - dvp_r
dgmse = gmse - gmse_r

# compute components
prefix = -2*np.pi*a**2*np.cos(np.radians(clat))
dcdvp = prefix * ddvp * gmse_r
dcgmse = prefix * dvp_r * dgmse
dcres = dvmte - dcdvp - dcgmse

# save file as netCDF
file_diff_vmte = Dataset(path_diff_vmte, "w", format='NETCDF4_CLASSIC')
file_diff_dvp = Dataset(path_diff_dvp, "w", format='NETCDF4_CLASSIC')
file_diff_gmse = Dataset(path_diff_gmse, "w", format='NETCDF4_CLASSIC')
file_diff_res = Dataset(path_diff_res, "w", format='NETCDF4_CLASSIC')

# copy attributes from r1
file_diff_vmte.setncatts(file_full_vmte.__dict__)
file_diff_dvp.setncatts(file_full_vmte.__dict__)
file_diff_gmse.setncatts(file_full_vmte.__dict__)
file_diff_res.setncatts(file_full_vmte.__dict__)

# copy dimensions from ps file
for name, dimension in file_full_vmte.dimensions.items():
    file_diff_vmte.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_diff_dvp.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_diff_gmse.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_diff_res.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_full_vmte.variables.items():
    if any(name in s for s in ['vmte']):
        continue
    
    x = file_diff_vmte.createVariable(name, variable.datatype, variable.dimensions)
    file_diff_vmte[name].setncatts(file_full_vmte[name].__dict__)
    file_diff_vmte[name][:] = file_full_vmte[name][:]

    x = file_diff_dvp.createVariable(name, variable.datatype, variable.dimensions)
    file_diff_dvp[name].setncatts(file_full_vmte[name].__dict__)
    file_diff_dvp[name][:] = file_full_vmte[name][:]

    x = file_diff_gmse.createVariable(name, variable.datatype, variable.dimensions)
    file_diff_gmse[name].setncatts(file_full_vmte[name].__dict__)
    file_diff_gmse[name][:] = file_full_vmte[name][:]

    x = file_diff_res.createVariable(name, variable.datatype, variable.dimensions)
    file_diff_res[name].setncatts(file_full_vmte[name].__dict__)
    file_diff_res[name][:] = file_full_vmte[name][:]

diff_vmte = file_diff_vmte.createVariable('dvmte', 'f4', ("time"))
diff_vmte.units = "W"
diff_vmte.long_name = "vmte deviation from reference state (historical)"
diff_vmte[:] = dvmte
file_diff_vmte.close()

diff_dvp = file_diff_dvp.createVariable('dcdvp', 'f4', ("time"))
diff_dvp.units = "W"
diff_dvp.long_name = "contribution of vmte deviation from reference state (historical) due to diffusivity change"
diff_dvp[:] = dcdvp
file_diff_dvp.close()

diff_gmse = file_diff_gmse.createVariable('dcgmse', 'f4', ("time"))
diff_gmse.units = "W"
diff_gmse.long_name = "contribution of vmte deviation from reference state (historical) due to MSE gradient change"
diff_gmse[:] = dcgmse
file_diff_gmse.close()

diff_res = file_diff_res.createVariable('dcres', 'f4', ("time"))
diff_res.units = "W"
diff_res.long_name = "contribution of vmte deviation from reference state (historical) due to nonlinear interactions"
diff_res[:] = dcres
file_diff_res.close()

