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
path_aht = sys.argv[1]
path_vmmmc = sys.argv[2]
path_vmse = sys.argv[3]
path_vmte = sys.argv[4]
path_raht = sys.argv[5]
path_rvmmmc = sys.argv[6]
path_rvmse = sys.argv[7]
path_rvmte = sys.argv[8]
path_daht = sys.argv[9]
path_dvmmmc = sys.argv[10]
path_dvmse = sys.argv[11]
path_dvmte = sys.argv[12]

# open files
print('Loading aht...')
file_aht = Dataset(path_aht, 'r')
print('Done.\n')
print('Loading vmmmc...')
file_vmmmc = Dataset(path_vmmmc, 'r')
print('Done.\n')
print('Loading vmse...')
file_vmse = Dataset(path_vmse, 'r')
print('Done.\n')
print('Loading vmte...')
file_vmte = Dataset(path_vmte, 'r')
print('Done.\n')
print('Loading raht...')
file_raht = Dataset(path_raht, 'r')
print('Done.\n')
print('Loading rvmmmc...')
file_rvmmmc = Dataset(path_rvmmmc, 'r')
print('Done.\n')
print('Loading rvmse...')
file_rvmse = Dataset(path_rvmse, 'r')
print('Done.\n')
print('Loading rvmte...')
file_rvmte = Dataset(path_rvmte, 'r')
print('Done.\n')

# read data
aht = file_aht.variables['aht'][:] # (mon x lat)
vmmmc = file_vmmmc.variables['vmmmc'][:] # (mon x lat)
vmse = file_vmse.variables['vmse'][:] # (mon x lat)
vmte = file_vmte.variables['vmte'][:] # (mon x lat)
raht = file_raht.variables['aht'][:] # (mon x lat)
rvmmmc = file_rvmmmc.variables['vmmmc'][:] # (mon x lat)
rvmse = file_rvmse.variables['vmse'][:] # (mon x lat)
rvmte = file_rvmte.variables['vmte'][:] # (mon x lat)

# take the differences
diff_aht = aht - raht
diff_vmmmc = vmmmc - rvmmmc
diff_vmse = vmse - rvmse
diff_vmte = vmte - rvmte

# save file as netCDF
file_daht = Dataset(path_daht, "w", format='NETCDF4_CLASSIC')
file_dvmmmc = Dataset(path_dvmmmc, "w", format='NETCDF4_CLASSIC')
file_dvmse = Dataset(path_dvmse, "w", format='NETCDF4_CLASSIC')
file_dvmte = Dataset(path_dvmte, "w", format='NETCDF4_CLASSIC')

# copy attributes from aht
file_daht.setncatts(file_aht.__dict__)
file_dvmmmc.setncatts(file_vmmmc.__dict__)
file_dvmse.setncatts(file_vmse.__dict__)
file_dvmte.setncatts(file_vmte.__dict__)

# copy dimensions from ps file
for name, dimension in file_aht.dimensions.items():
    # if any(name in s for s in ['lon']):
    #     continue
    file_daht.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_dvmmmc.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_dvmse.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_dvmte.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_aht.variables.items():
    if any(name in s for s in ['aht']):
        continue
    
    x = file_daht.createVariable(name, variable.datatype, variable.dimensions)
    file_daht[name].setncatts(file_aht[name].__dict__)
    file_daht[name][:] = file_aht[name][:]
    
    x = file_dvmmmc.createVariable(name, variable.datatype, variable.dimensions)
    file_dvmmmc[name].setncatts(file_aht[name].__dict__)
    file_dvmmmc[name][:] = file_aht[name][:]

    x = file_dvmse.createVariable(name, variable.datatype, variable.dimensions)
    file_dvmse[name].setncatts(file_aht[name].__dict__)
    file_dvmse[name][:] = file_aht[name][:]

    x = file_dvmte.createVariable(name, variable.datatype, variable.dimensions)
    file_dvmte[name].setncatts(file_aht[name].__dict__)
    file_dvmte[name][:] = file_aht[name][:]

daht = file_daht.createVariable('daht', 'f4', ("time","lat"))
daht.units = "W"
daht.long_name = "change in vertically integrated moist static energy flux transport relative to historical"
daht[:,:] = diff_aht
file_daht.close()

dvmmmc = file_dvmmmc.createVariable('dvmmmc', 'f4', ("time","lat"))
dvmmmc.units = "W"
dvmmmc.long_name = "change in vertically integrated moist static energy flux transport due to mean meridional circulation relative to historical"
dvmmmc[:,:] = diff_vmmmc
file_dvmmmc.close()

dvmse = file_dvmse.createVariable('dvmse', 'f4', ("time","lat"))
dvmse.units = "W"
dvmse.long_name = "change in vertically integrated moist static energy flux transport due to stationary eddies relative to historical"
dvmse[:,:] = diff_vmse
file_dvmse.close()

dvmte = file_dvmte.createVariable('dvmte', 'f4', ("time","lat"))
dvmte.units = "W"
dvmte.long_name = "change in vertically integrated moist static energy flux transport due to transient eddies relative to historical"
dvmte[:,:] = diff_vmte
file_dvmte.close()
