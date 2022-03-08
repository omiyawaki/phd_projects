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
path_saht = sys.argv[1]
path_vsmmc = sys.argv[2]
path_vsse = sys.argv[3]
path_vste = sys.argv[4]
path_rsaht = sys.argv[5]
path_rvsmmc = sys.argv[6]
path_rvsse = sys.argv[7]
path_rvste = sys.argv[8]
path_dsaht = sys.argv[9]
path_dvsmmc = sys.argv[10]
path_dvsse = sys.argv[11]
path_dvste = sys.argv[12]

# open files
print('Loading saht...')
file_saht = Dataset(path_saht, 'r')
print('Done.\n')
print('Loading vsmmc...')
file_vsmmc = Dataset(path_vsmmc, 'r')
print('Done.\n')
print('Loading vsse...')
file_vsse = Dataset(path_vsse, 'r')
print('Done.\n')
print('Loading vste...')
file_vste = Dataset(path_vste, 'r')
print('Done.\n')
print('Loading rsaht...')
file_rsaht = Dataset(path_rsaht, 'r')
print('Done.\n')
print('Loading rvsmmc...')
file_rvsmmc = Dataset(path_rvsmmc, 'r')
print('Done.\n')
print('Loading rvsse...')
file_rvsse = Dataset(path_rvsse, 'r')
print('Done.\n')
print('Loading rvste...')
file_rvste = Dataset(path_rvste, 'r')
print('Done.\n')

# read data
saht = file_saht.variables['saht'][:] # (mon x lat)
vsmmc = file_vsmmc.variables['vsmmc'][:] # (mon x lat)
vsse = file_vsse.variables['vsse'][:] # (mon x lat)
vste = file_vste.variables['vste'][:] # (mon x lat)
rsaht = file_rsaht.variables['saht'][:] # (mon x lat)
rvsmmc = file_rvsmmc.variables['vsmmc'][:] # (mon x lat)
rvsse = file_rvsse.variables['vsse'][:] # (mon x lat)
rvste = file_rvste.variables['vste'][:] # (mon x lat)

# take the differences
diff_saht = saht - rsaht
diff_vsmmc = vsmmc - rvsmmc
diff_vsse = vsse - rvsse
diff_vste = vste - rvste

# save file as netCDF
file_dsaht = Dataset(path_dsaht, "w", format='NETCDF4_CLASSIC')
file_dvsmmc = Dataset(path_dvsmmc, "w", format='NETCDF4_CLASSIC')
file_dvsse = Dataset(path_dvsse, "w", format='NETCDF4_CLASSIC')
file_dvste = Dataset(path_dvste, "w", format='NETCDF4_CLASSIC')

# copy attributes from saht
file_dsaht.setncatts(file_saht.__dict__)
file_dvsmmc.setncatts(file_vsmmc.__dict__)
file_dvsse.setncatts(file_vsse.__dict__)
file_dvste.setncatts(file_vste.__dict__)

# copy dimensions from ps file
for name, dimension in file_saht.dimensions.items():
    # if any(name in s for s in ['lon']):
    #     continue
    file_dsaht.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_dvsmmc.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_dvsse.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
    file_dvste.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_saht.variables.items():
    if any(name in s for s in ['saht']):
        continue
    
    x = file_dsaht.createVariable(name, variable.datatype, variable.dimensions)
    file_dsaht[name].setncatts(file_saht[name].__dict__)
    file_dsaht[name][:] = file_saht[name][:]
    
    x = file_dvsmmc.createVariable(name, variable.datatype, variable.dimensions)
    file_dvsmmc[name].setncatts(file_saht[name].__dict__)
    file_dvsmmc[name][:] = file_saht[name][:]

    x = file_dvsse.createVariable(name, variable.datatype, variable.dimensions)
    file_dvsse[name].setncatts(file_saht[name].__dict__)
    file_dvsse[name][:] = file_saht[name][:]

    x = file_dvste.createVariable(name, variable.datatype, variable.dimensions)
    file_dvste[name].setncatts(file_saht[name].__dict__)
    file_dvste[name][:] = file_saht[name][:]

dsaht = file_dsaht.createVariable('dsaht', 'f4', ("time","lat"))
dsaht.units = "W"
dsaht.long_name = "change in vertically integrated dry static energy flux transport relative to historical"
dsaht[:,:] = diff_saht
file_dsaht.close()

dvsmmc = file_dvsmmc.createVariable('dvsmmc', 'f4', ("time","lat"))
dvsmmc.units = "W"
dvsmmc.long_name = "change in vertically integrated dry static energy flux transport due to mean meridional circulation relative to historical"
dvsmmc[:,:] = diff_vsmmc
file_dvsmmc.close()

dvsse = file_dvsse.createVariable('dvsse', 'f4', ("time","lat"))
dvsse.units = "W"
dvsse.long_name = "change in vertically integrated dry static energy flux transport due to stationary eddies relative to historical"
dvsse[:,:] = diff_vsse
file_dvsse.close()

dvste = file_dvste.createVariable('dvste', 'f4', ("time","lat"))
dvste.units = "W"
dvste.long_name = "change in vertically integrated dry static energy flux transport due to transient eddies relative to historical"
dvste[:,:] = diff_vste
file_dvste.close()
