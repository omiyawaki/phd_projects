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
path_rlus = sys.argv[1]
path_rlds = sys.argv[2]
path_rsds = sys.argv[3]
path_rsus = sys.argv[4]
path_hfls = sys.argv[5]
path_hfss = sys.argv[6]
path_fsfc = sys.argv[7]

# open files
print('Loading rlus file...')
file_rlus = Dataset(path_rlus, 'r')
print('Done.\n')

print('Loading rlds file...')
file_rlds = Dataset(path_rlds, 'r')
print('Done.\n')

print('Loading rsds file...')
file_rsds = Dataset(path_rsds, 'r')
print('Done.\n')

print('Loading rsus file...')
file_rsus = Dataset(path_rsus, 'r')
print('Done.\n')

print('Loading hfls file...')
file_hfls = Dataset(path_hfls, 'r')
print('Done.\n')

print('Loading hfss file...')
file_hfss = Dataset(path_hfss, 'r')
print('Done.\n')

# read data
rlus = file_rlus.variables['rlus'][:] # (mon x lat x lon)
rlds = file_rlds.variables['rlds'][:] # (mon x lat x lon)
rsds = file_rsds.variables['rsds'][:] # (mon x lat x lon)
rsus = file_rsus.variables['rsus'][:] # (mon x lat x lon)
hfls = file_hfls.variables['hfls'][:] # (mon x lat x lon)
hfss = file_hfss.variables['hfss'][:] # (mon x lat x lon)

# compute net surface flux
sfc = - (rsus - rsds + rlus - rlds) - hfls - hfss # positive heats surface

# save file as netCDF
file_fsfc = Dataset(path_fsfc, "w", format='NETCDF4_CLASSIC')

# copy attributes from rsus file
file_fsfc.setncatts(file_rsus.__dict__)

# copy dimensions from rsus file
for name, dimension in file_rsus.dimensions.items():
    file_fsfc.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from rsus file
for name, variable in file_rsus.variables.items():
    if any(name in s for s in ['rsus']):
        continue
    
    x = file_fsfc.createVariable(name, variable.datatype, variable.dimensions)
    file_fsfc[name].setncatts(file_rsus[name].__dict__)
    file_fsfc[name][:] = file_rsus[name][:]
    
fsfc = file_fsfc.createVariable('fsfc', 'f4', ("time","lat","lon"))
fsfc.units = "W m**-2"
fsfc.long_name = "net surface flux"
fsfc[:] = sfc

file_fsfc.close()
