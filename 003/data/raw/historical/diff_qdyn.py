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
path_qaht = sys.argv[1]
path_vqmmc = sys.argv[2]
path_vqse = sys.argv[3]
path_vqte = sys.argv[4]
path_rqaht = sys.argv[5]
path_rvqmmc = sys.argv[6]
path_rvqse = sys.argv[7]
path_rvqte = sys.argv[8]
path_dqaht = sys.argv[9]
path_dvqmmc = sys.argv[10]
path_dvqse = sys.argv[11]
path_dvqte = sys.argv[12]

# open files
print('Loading qaht...')
file_qaht = Dataset(path_qaht, 'r')
print('Done.\n')
print('Loading vqmmc...')
file_vqmmc = Dataset(path_vqmmc, 'r')
print('Done.\n')
print('Loading vqse...')
file_vqse = Dataset(path_vqse, 'r')
print('Done.\n')
print('Loading vqte...')
file_vqte = Dataset(path_vqte, 'r')
print('Done.\n')
print('Loading rqaht...')
file_rqaht = Dataset(path_rqaht, 'r')
print('Done.\n')
print('Loading rvqmmc...')
file_rvqmmc = Dataset(path_rvqmmc, 'r')
print('Done.\n')
print('Loading rvqse...')
file_rvqse = Dataset(path_rvqse, 'r')
print('Done.\n')
print('Loading rvqte...')
file_rvqte = Dataset(path_rvqte, 'r')
print('Done.\n')

# read data
qaht = file_qaht.variables['qaht'][:] # (mon x lat)
vqmmc = file_vqmmc.variables['vqmmc'][:] # (mon x lat)
vqse = file_vqse.variables['vqse'][:] # (mon x lat)
vqte = file_vqte.variables['vqte'][:] # (mon x lat)
rqaht = file_rqaht.variables['qaht'][:] # (mon x lat)
rvqmmc = file_rvqmmc.variables['vqmmc'][:] # (mon x lat)
rvqse = file_rvqse.variables['vqse'][:] # (mon x lat)
rvqte = file_rvqte.variables['vqte'][:] # (mon x lat)

# take the differences
diff_qaht = qaht - rqaht
diff_vqmmc = vqmmc - rvqmmc
diff_vqse = vqse - rvqse
diff_vqte = vqte - rvqte

# save file as netCDF
file_dqaht = Dataset(path_dqaht, "w", format='NETCDF4_CLASSIC')
file_dvqmmc = Dataset(path_dvqmmc, "w", format='NETCDF4_CLASSIC')
file_dvqse = Dataset(path_dvqse, "w", format='NETCDF4_CLASSIC')
file_dvqte = Dataset(path_dvqte, "w", format='NETCDF4_CLASSIC')

# copy attributes from qaht
file_dqaht.setncatts(file_qaht.__dict__)
file_dvqmmc.setncatts(file_vqmmc.__dict__)
file_dvqse.setncatts(file_vqse.__dict__)
file_dvqte.setncatts(file_vqte.__dict__)

# copy dimensions from ps file
for name, dimension in file_qaht.dimensions.items():
    # if any(name in s for s in ['lon']):
    #     continue
    file_dqaht.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_dvqmmc.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_dvqse.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_dvqte.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_qaht.variables.items():
    if any(name in s for s in ['qaht']):
        continue
    
    x = file_dqaht.createVariable(name, variable.datatype, variable.dimensions)
    file_dqaht[name].setncatts(file_qaht[name].__dict__)
    file_dqaht[name][:] = file_qaht[name][:]
    
    x = file_dvqmmc.createVariable(name, variable.datatype, variable.dimensions)
    file_dvqmmc[name].setncatts(file_qaht[name].__dict__)
    file_dvqmmc[name][:] = file_qaht[name][:]

    x = file_dvqse.createVariable(name, variable.datatype, variable.dimensions)
    file_dvqse[name].setncatts(file_qaht[name].__dict__)
    file_dvqse[name][:] = file_qaht[name][:]

    x = file_dvqte.createVariable(name, variable.datatype, variable.dimensions)
    file_dvqte[name].setncatts(file_qaht[name].__dict__)
    file_dvqte[name][:] = file_qaht[name][:]

dqaht = file_dqaht.createVariable('dqaht', 'f4', ("time","lat"))
dqaht.units = "W"
dqaht.long_name = "change in vertically integrated latent energy flux transport relative to historical"
dqaht[:,:] = diff_qaht
file_dqaht.close()

dvqmmc = file_dvqmmc.createVariable('dvqmmc', 'f4', ("time","lat"))
dvqmmc.units = "W"
dvqmmc.long_name = "change in vertically integrated latent energy flux transport due to mean meridional circulation relative to historical"
dvqmmc[:,:] = diff_vqmmc
file_dvqmmc.close()

dvqse = file_dvqse.createVariable('dvqse', 'f4', ("time","lat"))
dvqse.units = "W"
dvqse.long_name = "change in vertically integrated latent energy flux transport due to stationary eddies relative to historical"
dvqse[:,:] = diff_vqse
file_dvqse.close()

dvqte = file_dvqte.createVariable('dvqte', 'f4', ("time","lat"))
dvqte.units = "W"
dvqte.long_name = "change in vertically integrated latent energy flux transport due to transient eddies relative to historical"
dvqte[:,:] = diff_vqte
file_dvqte.close()
