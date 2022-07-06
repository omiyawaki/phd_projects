import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps and va files
path_psi = sys.argv[1]
path_pp0 = sys.argv[2]

# open files
file_psi = Dataset(path_psi, 'r')

# read data
psi = np.squeeze(file_psi.variables['psic'][:]) # (time x lev x lat)
p = np.squeeze(file_psi.variables['plev'][:])

# compute pressure where streamfunction is 0
pp0 = np.zeros([psi.shape[0], psi.shape[2]]) # (time x lat)
for t in tqdm(range(psi.shape[0])):
    for lat in range(psi.shape[2]):
        psi_l = psi[t,:,lat]
        f = np.polynomial.polynomial.Polynomial.fit(p, psi_l, 3) # cubic polynomial fit
        pcand = f.roots() # solve for the roots
        pcand = pcand[pcand > 0] # discard negative pressures (not physical)
        if len(pcand) == 0: # if there are no solutions, move on
            continue
        pcand = np.min(pcand) # if there are more than one, choose lower pressure
        pp0[t, lat] = pcand

print(pp0[0,:])

# save file as netCDF
file_pp0 = Dataset(path_pp0, "w", format='NETCDF4_CLASSIC')

# copy attributes from psi file
file_pp0.setncatts(file_psi.__dict__)

# copy dimensions from psi file
for name, dimension in file_psi.dimensions.items():
    if any(name in s for s in ['plev']):
        continue
    file_pp0.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from psi file
for name, variable in file_psi.variables.items():
    if any(name in s for s in ['psic', 'plev', 'plev_bnds']):
        continue
    
    x = file_pp0.createVariable(name, variable.datatype, variable.dimensions)
    file_pp0[name].setncatts(file_psi[name].__dict__)
    file_pp0[name][:] = file_psi[name][:]
    
pa = file_pp0.createVariable('pp0', 'f4', ("time","lat"))
pa.units = "Pa"
pa.long_name = "pressure level where streamfunction changes sign (where psi=0)"
pa[:,:] = pp0

file_pp0.close()
