import sys
import numpy as np
from scipy import interpolate, integrate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;
ampf = 25 # unitless, amplification factor = Cw/Ci
Lf = 3.337e5 # J kg**-1, latent heat of fusion
rhoi = 905 # kg m**-3, density of ice

# paths to ps and mse files
path_ref = sys.argv[1]
path_sic = sys.argv[2]
path_sit = sys.argv[3]
path_fice = sys.argv[4]
path_fnoice = sys.argv[5]
path_qflux = sys.argv[6]

# open files
file_ref = Dataset(path_ref, 'r')
file_sic = Dataset(path_sic, 'r')
file_sit = Dataset(path_sit, 'r')
file_fice = Dataset(path_fice, 'r')
file_fnoice = Dataset(path_fnoice, 'r')

# read data
ref = file_ref.variables['aflux'][:] # (mon x lat x lon)
sic = file_sic.variables['seaice'][:] # (mon x lat x lon) sea ice concentration (1)
sit = file_sit.variables['siced'][:] # (mon x lat x lon) sea ice thickness (m)
fice = file_fice.variables['fsfc'][:] # (mon x lat x lon) net surface flux w/ ice
fnoice = file_fnoice.variables['fsfc'][:] # (mon x lat x lon) net surface flux w/o i

# compute ice thickness tendency
dt = 30*86400
stend = (sit[2:,...] - sit[:-2,...])/(2*dt)
stend = np.insert(stend, 0, ( (sit[1,...] - sit[-1,...])/(2*dt) )[None,...], axis=0)
stend = np.append(stend, ( (sit[0,...] - sit[-2,...])/(2*dt) )[None,...], axis=0)

# compute qflux to apply
f = np.ones_like(sic)
f[sic==1] = ampf # q-flux amplification factor where sea ice exists
f[stend < 0] = 1 # where sea ice is melting set the heat capacity equal to water
# f[stend < 0] = 0 # where sea ice is melting set the heat capacity to infinity
aflux = fice - f*(fice + Lf*rhoi*stend)

# save file as netCDF
file_qflux = Dataset(path_qflux, "w", format='NETCDF4_CLASSIC')

# copy attributes from ref file
file_qflux.setncatts(file_ref.__dict__)

# copy dimensions from ref file
for name, dimension in file_ref.dimensions.items():
    file_qflux.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables from ref file
for name, variable in file_ref.variables.items():
    if any(name in s for s in ['aflux']):
        continue
    
    x = file_qflux.createVariable(name, variable.datatype, variable.dimensions)
    file_qflux[name].setncatts(file_ref[name].__dict__)
    file_qflux[name][:] = file_ref[name][:]
    
qflux = file_qflux.createVariable('aflux', 'f4', ("time","lat","lon"))
qflux.units = "W m**-2"
qflux.long_name = "LW flux over water"
qflux[:] = aflux

file_qflux.close()
