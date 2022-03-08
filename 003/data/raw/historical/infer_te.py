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
path_vmmmc = sys.argv[1]
path_vmse = sys.argv[2]
path_aht = sys.argv[3]
path_vmte = sys.argv[4]

# open files
print('Loading vmmmc file...')
file_vmmmc = Dataset(path_vmmmc, 'r')
print('Done.\n')

print('Loading vmse file...')
file_vmse = Dataset(path_vmse, 'r')
print('Done.\n')

print('Loading aht file...')
file_aht = Dataset(path_aht, 'r')
print('Done.\n')

# read data
vmmmc = file_vmmmc.variables['vmmmc'][:] # (mon x lat)
vmse = file_vmse.variables['vmse'][:] # (mon x lat)
aht = file_aht.variables['aht'][:] # (mon x lat)

# infer transient eddy transport as residual
vm_te_vint = aht - vmmmc - vmse

# save file as netCDF
file_vmte = Dataset(path_vmte, "w", format='NETCDF4_CLASSIC')

# copy attributes from vmmmcfile
file_vmte.setncatts(file_vmmmc.__dict__)

# copy dimensions from vmmmc file
for name, dimension in file_vmmmc.dimensions.items():
    # if any(name in s for s in ['lon']):
    #     continue
    file_vmte.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from vmmmc file
for name, variable in file_vmmmc.variables.items():
    if any(name in s for s in ['vmmmc']):
        continue
    
    x = file_vmte.createVariable(name, variable.datatype, variable.dimensions)
    file_vmte[name].setncatts(file_vmmmc[name].__dict__)
    file_vmte[name][:] = file_vmmmc[name][:]
    
vmte = file_vmte.createVariable('vmte', 'f4', ("time","lat"))
vmte.units = "W"
vmte.long_name = "vertically integrated moist static energy flux transport due to transient eddies"
vmte[:,:] = vm_te_vint

file_vmte.close()
