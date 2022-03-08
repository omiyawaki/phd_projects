import sys
import numpy as np
from scipy import interpolate, integrate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to various transport files
path_aht = sys.argv[1]
# path_vmmmc = sys.argv[2]
# path_vmse = sys.argv[3]
path_vmte = sys.argv[4]

path_qaht = sys.argv[5]
# path_vqmmc = sys.argv[6]
# path_vqse = sys.argv[7]
path_vqte = sys.argv[8]

path_saht = sys.argv[9]
# path_vsmmc = sys.argv[10]
# path_vsse = sys.argv[11]
path_vste = sys.argv[12]

# open files
file_aht = Dataset(path_aht, 'r')
# file_vmmmc = Dataset(path_vmmmc, 'r')
# file_vmse = Dataset(path_vmse, 'r')
file_vmte = Dataset(path_vmte, 'r')

file_qaht = Dataset(path_qaht, 'r')
# file_vqmmc = Dataset(path_vqmmc, 'r')
# file_vqse = Dataset(path_vqse, 'r')
file_vqte = Dataset(path_vqte, 'r')

# read data
aht = file_aht.variables['aht'][:] # (mon x lat x lon)
# vmmmc = file_vmmmc.variables['vmmmc'][:] # (mon x lat x lon)
# vmse = file_vmse.variables['vmse'][:] # (mon x lat x lon)
vmte = file_vmte.variables['vmte'][:] # (mon x lat x lon)

qaht = file_qaht.variables['qaht'][:] # (mon x lat x lon)
# vqmmc = file_vqmmc.variables['vqmmc'][:] # (mon x lat x lon)
# vqse = file_vqse.variables['vqse'][:] # (mon x lat x lon)
vqte = file_vqte.variables['vqte'][:] # (mon x lat x lon)

# infer DSE energy flux divergence for all components
saht = aht - qaht
# vsmmc = vmmmc - vqmmc
# vsse = vmse - vqse
vste = vmte - vqte

# save file as netCDF
file_saht = Dataset(path_saht, "w", format='NETCDF4_CLASSIC')
file_vste = Dataset(path_vste, "w", format='NETCDF4_CLASSIC')

# copy attributes from aht file
file_saht.setncatts(file_aht.__dict__)
file_vste.setncatts(file_aht.__dict__)

# copy dimensions from aht file
for name, dimension in file_aht.dimensions.items():
    file_saht.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_vste.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from aht file
for name, variable in file_aht.variables.items():
    if any(name in s for s in ['aht']):
        continue
    
    x = file_saht.createVariable(name, variable.datatype, variable.dimensions)
    file_saht[name].setncatts(file_aht[name].__dict__)
    file_saht[name][:] = file_aht[name][:]

    x = file_vste.createVariable(name, variable.datatype, variable.dimensions)
    file_vste[name].setncatts(file_aht[name].__dict__)
    file_vste[name][:] = file_aht[name][:]

vsE = file_saht.createVariable('saht', 'f4', ("time","lat"))
vsE.units = "W"
vsE.long_name = "vertically integrated dry static energy flux transport"
vsE[:,:] = saht

te = file_vste.createVariable('vste', 'f4', ("time","lat"))
te.units = "W"
te.long_name = "vertically integrated dry transient eddy energy flux transport"
te[:,:] = vste

file_saht.close()
file_vste.close()
