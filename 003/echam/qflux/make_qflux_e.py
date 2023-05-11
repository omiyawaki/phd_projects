import sys
import numpy as np
from scipy import interpolate, integrate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;
rhosea = 1025 # sea water density, kg m**-3 (from Giorgetta et al 2013)
csea = 3994 # sea water specific heat capacity, J kg**-1 K**-1 (from Giorgetta et al 2013)
h = 40 # slab ocean mixed layer depth
albedoo = 0.07 # ocean albedoedo

# paths to ps and mse files
path_temp = sys.argv[1]
path_fsfc = sys.argv[2]
path_swref = sys.argv[3]
path_albedoi = sys.argv[4]
path_sradsd = sys.argv[5]
path_qflx = sys.argv[6]
path_qsw = sys.argv[7]

# open files
file_temp = Dataset(path_temp, 'r')
file_fsfc = Dataset(path_fsfc, 'r')
file_swref = Dataset(path_swref, 'r')
file_albedoi = Dataset(path_albedoi, 'r')
file_sradsd = Dataset(path_sradsd, 'r')

# read data
ts = file_temp.variables['tsurf'][:] # target surface temperature, K (mon x lat x lon)
fsfc = file_fsfc.variables['fsfc'][:] # target net surface energy flux, W m**-2 (mon x lat x lon)
swref = file_swref.variables['swref'][:] # target net surface energy flux, W m**-2 (mon x lat x lon)
albedoi = file_albedoi.variables['albedo'][:] # target net surface energy flux, W m**-2 (mon x lat x lon)
sradsd = file_sradsd.variables['sradsd'][:] # target net surface energy flux, W m**-2 (mon x lat x lon)

ts=ts.filled(fill_value=np.nan)
fsfc=fsfc.filled(fill_value=np.nan)
swref=swref.filled(fill_value=np.nan)
albedoi=albedoi.filled(fill_value=np.nan)
sradsd=sradsd.filled(fill_value=np.nan)

# compute surf temp tendency
dt = 30*86400
ttend = (ts[2:,...] - ts[:-2,...])/(2*dt)
ttend = np.insert(ttend, 0, ( (ts[1,...] - ts[-1,...])/(2*dt) )[None,...], axis=0)
ttend = np.append(ttend, ( (ts[0,...] - ts[-2,...])/(2*dt) )[None,...], axis=0)

# compute qflux
qc=fsfc-rhosea*h*csea*ttend 
qs=sradsd/(1-albedoo*swref)*(albedoi-albedoo)*(1-swref)
qs[np.isnan(qs)]=0
qflux = qc + qs

##################################################
# save qflx as netCDF
##################################################
file_qflx = Dataset(path_qflx, "w", format='NETCDF4_CLASSIC')

# copy attributes from temp file
file_qflx.setncatts(file_temp.__dict__)

# copy dimensions from temp file
for name, dimension in file_temp.dimensions.items():
    file_qflx.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from temp file
for name, variable in file_temp.variables.items():
    if any(name in s for s in ['tsurf']):
        continue
    
    x = file_qflx.createVariable(name, variable.datatype, variable.dimensions)
    file_qflx[name].setncatts(file_temp[name].__dict__)
    file_qflx[name][:] = file_temp[name][:]
    
qflx = file_qflx.createVariable('aflux', 'f4', ("time","lat","lon"))
qflx.units = "W m**-2"
qflx.long_name = "LW flux over water"
qflx[:] = qflux

file_qflx.close()

file_qsw = Dataset(path_qsw, "w", format='NETCDF4_CLASSIC')

# copy attributes from temp file
file_qsw.setncatts(file_temp.__dict__)

# copy dimensions from temp file
for name, dimension in file_temp.dimensions.items():
    file_qsw.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from temp file
for name, variable in file_temp.variables.items():
    if any(name in s for s in ['tsurf']):
        continue
    
    x = file_qsw.createVariable(name, variable.datatype, variable.dimensions)
    file_qsw[name].setncatts(file_temp[name].__dict__)
    file_qsw[name][:] = file_temp[name][:]
    
qsw = file_qsw.createVariable('aflux', 'f4', ("time","lat","lon"))
qsw.units = "W m**-2"
qsw.long_name = "LW flux over water"
qsw[:] = qs

file_qsw.close()

