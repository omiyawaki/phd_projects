import sys
import numpy as np
from scipy import interpolate, integrate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;
d = 40 # mixed layer depth
# rhow = 1000 # kg m**-3, density of water
# cw = 4218 # J K**-1 kg**-1 specific heat of liquid water
rhow = 1025 # kg m**-3, density of sea water
cw = 3994 # J K**-1 kg**-1 specific heat of sea water
alo = 0.07 # surface albedo of ocean

# paths to ps and mse files
path_ref = sys.argv[1]
path_dtt = sys.argv[2]
path_swn = sys.argv[3]
path_alb = sys.argv[4]
path_lw = sys.argv[5]
path_lh = sys.argv[6]
path_sh = sys.argv[7]
path_qflux = sys.argv[8]

# open files
file_ref = Dataset(path_ref, 'r')
file_dtt = Dataset(path_dtt, 'r')
file_swn = Dataset(path_swn, 'r')
file_alb = Dataset(path_alb, 'r')
file_lw = Dataset(path_lw, 'r')
file_lh = Dataset(path_lh, 'r')
file_sh = Dataset(path_sh, 'r')

# read data
ref = file_ref.variables['aflux'][:] # (mon x lat x lon)
dtt = file_dtt.variables['tsurf'][:] # (mon x lat x lon) surface temperature (K)
swn = file_swn.variables['srads'][:] # (mon x lat x lon) net surface sw (wm-2)
alb = file_alb.variables['albedo'][:] # (mon x lat x lon) surface albedo (1)
lw = file_lw.variables['trads'][:] # (mon x lat x lon) net surface lw flux (wm-2)
lh = file_lh.variables['ahfl'][:] # (mon x lat x lon) surf latent heat flux (wm-2)
sh = file_sh.variables['ahfs'][:] # (mon x lat x lon) surf sensible heat flux (wm-2)

# compute surf temp tendency
dt = 30*86400
ttend = (dtt[2:,...] - dtt[:-2,...])/(2*dt)
ttend = np.insert(ttend, 0, ( (dtt[1,...] - dtt[-1,...])/(2*dt) )[None,...], axis=0)
ttend = np.append(ttend, ( (dtt[0,...] - dtt[-2,...])/(2*dt) )[None,...], axis=0)

# q flux
aflux = swn*(1-alo)/(1-alb) + lw + lh + sh  - rhow*cw*d*ttend

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
