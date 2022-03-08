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
path_vqmmc = sys.argv[1]
path_vqse = sys.argv[2]
path_qaht = sys.argv[3]
path_vqte = sys.argv[4]

# open files
print('Loading vqmmc file...')
file_vqmmc = Dataset(path_vqmmc, 'r')
print('Done.\n')

print('Loading vqse file...')
file_vqse = Dataset(path_vqse, 'r')
print('Done.\n')

print('Loading qaht file...')
file_qaht = Dataset(path_qaht, 'r')
print('Done.\n')

# read data
vqmmc = file_vqmmc.variables['vqmmc'][:] # (mon x lat)
vqse = file_vqse.variables['vqse'][:] # (mon x lat)
qaht = file_qaht.variables['qaht'][:] # (mon x lat)

# infer transient eddy transport as residual
vq_te_vint = qaht - vqmmc - vqse

# save file as netCDF
file_vqte = Dataset(path_vqte, "w", format='NETCDF4_CLASSIC')

# copy attributes from vqmmcfile
file_vqte.setncatts(file_vqmmc.__dict__)

# copy dimensions from vqmmc file
for name, dimension in file_vqmmc.dimensions.items():
    # if any(name in s for s in ['lon']):
    #     continue
    file_vqte.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from vqmmc file
for name, variable in file_vqmmc.variables.items():
    if any(name in s for s in ['vqmmc']):
        continue
    
    x = file_vqte.createVariable(name, variable.datatype, variable.dimensions)
    file_vqte[name].setncatts(file_vqmmc[name].__dict__)
    file_vqte[name][:] = file_vqmmc[name][:]
    
vqte = file_vqte.createVariable('vqte', 'f4', ("time","lat"))
vqte.units = "W"
vqte.long_name = "vertically integrated latent energy flux transport due to transient eddies"
vqte[:,:] = vq_te_vint

file_vqte.close()
