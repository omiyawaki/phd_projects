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
path_sit = sys.argv[2]
path_melt = sys.argv[3]
path_yice = sys.argv[4]
path_nice = sys.argv[5]
path_qflx = sys.argv[6]
path_deff = sys.argv[7]
path_ampf = sys.argv[8]

# open files
file_sic = Dataset(path_sic, 'r')
file_sit = Dataset(path_sit, 'r')
file_melt = Dataset(path_melt, 'r')
file_yice = Dataset(path_yice, 'r')
file_nice = Dataset(path_nice, 'r')

# read data
sic = file_sic.variables['seaice'][:] # sea ice fraction, unitless (mon x lat x lon)
sit = file_sit.variables['siced'][:] # sea ice depth, m (mon x lat x lon)
melt = file_melt.variables['ahfres'][:] # energy flux associated with melting sea ice, W m**-2 (mon x lat x lon)
fsfci = file_yice.variables['fsfc'][:] # net surface energy flux for ice run, W m**-2 (mon x lat x lon)
fsfcw = file_nice.variables['fsfc'][:] # net surface energy flux for no ice run, W m**-2 (mon x lat x lon)

# compute effective mixed layer depth (mldeff)
Ci = ci*rhoi*sit # ice heat capacity J m**-2 K**-1
Cw = cw*rhow*mldsic # water heat capacity J m**-2 K**-1
mldeff = ( sic*Ci + (1-sic)*Cw ) / (cw*rhow) # effective MLD
F = mldsic / mldeff # q-flux amplification factor
F[melt>0] = 1 # don't amplify where melt is occuring

# enthalpy of sea ice
hi = Lf*rhoi*sit # J m**-2

# energy associated with melting sea ice
qmelt = melt*30*86400 # J m**-2

# residual energy after complete melt
qres = qmelt - hi

# Q flux to apply to no ice run to target ice Fsfc 
qflux = - (F*fsfci - fsfcw)

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

##################################################
# save DEFF as netCDF
##################################################
file_deff = Dataset(path_deff, "w", format='NETCDF4_CLASSIC')

# copy attributes from sic file
file_deff.setncatts(file_sic.__dict__)

# copy dimensions from sic file
for name, dimension in file_sic.dimensions.items():
    file_deff.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from sic file
for name, variable in file_sic.variables.items():
    if any(name in s for s in ['seaice']):
        continue
    
    x = file_deff.createVariable(name, variable.datatype, variable.dimensions)
    file_deff[name].setncatts(file_sic[name].__dict__)
    file_deff[name][:] = file_sic[name][:]
    
deff = file_deff.createVariable('deff', 'f4', ("time","lat","lon"))
deff.units = "m"
deff.long_name = "effective mixed layer depth"
deff[:] = mldeff

file_deff.close()

##################################################
# save ampf as netCDF
##################################################
file_ampf = Dataset(path_ampf, "w", format='NETCDF4_CLASSIC')

# copy attributes from sic file
file_ampf.setncatts(file_sic.__dict__)

# copy dimensions from sic file
for name, dimension in file_sic.dimensions.items():
    file_ampf.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from sic file
for name, variable in file_sic.variables.items():
    if any(name in s for s in ['seaice']):
        continue
    
    x = file_ampf.createVariable(name, variable.datatype, variable.dimensions)
    file_ampf[name].setncatts(file_sic[name].__dict__)
    file_ampf[name][:] = file_sic[name][:]
    
ampf = file_ampf.createVariable('ampf', 'f4', ("time","lat","lon"))
ampf.units = "unitless"
ampf.long_name = "q flux amplification factor"
ampf[:] = F

file_ampf.close()
