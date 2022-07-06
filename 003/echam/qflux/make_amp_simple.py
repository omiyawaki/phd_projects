import sys
import numpy as np
from scipy import interpolate, integrate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;
Lf = 3.337e5 # latent heat of fusion J kg**-1
rhow = 997; rhoi = 917; # density, kg m**-3
cw = 4182; ci = 2108; # specific heat capacity, J kg**-1 K**-1
mldsic = 40 # bare ocean mixed layer depth

# paths to ps and mse files
path_sic = sys.argv[1]
path_yice = sys.argv[2]
path_nice = sys.argv[3]
path_qflx = sys.argv[4]

# open files
file_sic = Dataset(path_sic, 'r')
file_yice = Dataset(path_yice, 'r')
file_nice = Dataset(path_nice, 'r')

# read data
sic = file_sic.variables['seaice'][:] # sea ice fraction, unitless (mon x lat x lon)
fsfci = file_yice.variables['fsfc'][:] # net surface energy flux for ice run, W m**-2 (mon x lat x lon)
fsfcw = file_nice.variables['fsfc'][:] # net surface energy flux for no ice run, W m**-2 (mon x lat x lon)

# compute effective mixed layer depth (mldeff)
F = np.ones_like(sic)
F[sic==1] = 10 # q-flux amplification factor where sea ice exists
qflux = fsfci*(1-F)

##################################################
# save qflx as netCDF
##################################################
file_qflx = Dataset(path_qflx, "w", format='NETCDF4_CLASSIC')

# copy attributes from sic file
file_qflx.setncatts(file_sic.__dict__)

# copy dimensions from sic file
for name, dimension in file_sic.dimensions.items():
    file_qflx.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from sic file
for name, variable in file_sic.variables.items():
    if any(name in s for s in ['seaice']):
        continue
    
    x = file_qflx.createVariable(name, variable.datatype, variable.dimensions)
    file_qflx[name].setncatts(file_sic[name].__dict__)
    file_qflx[name][:] = file_sic[name][:]
    
qflx = file_qflx.createVariable('aflux', 'f4', ("time","lat","lon"))
qflx.units = "W m**-2"
qflx.long_name = "LW flux over water"
qflx[:] = qflux

file_qflx.close()

