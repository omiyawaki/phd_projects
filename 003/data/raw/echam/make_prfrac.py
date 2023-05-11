import sys
import numpy as np
from scipy import interpolate, integrate
from scipy.ndimage import uniform_filter1d
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps and mse files
path_pr = sys.argv[1]
path_prc = sys.argv[2]
path_prl = sys.argv[3]
path_prfrac = sys.argv[4]

# open files
print('Loading pr...')
file_pr = Dataset(path_pr, 'r')
print('Done.\n')

print('Loading prc...')
file_prc = Dataset(path_prc, 'r')
print('Done.\n')

# read data
pr = file_pr.variables['precip'][:]
prc = file_prc.variables['aprc'][:]

# compute large scale precip and convective precipitation fraction
prfrac = prc/pr
prfrac[pr == 0] = np.nan

# save prl and prfrac as netCDF
file_prfrac = Dataset(path_prfrac, "w", format='NETCDF4_CLASSIC')

# copy attributes from original file
file_prfrac.setncatts(file_pr.__dict__)

# copy dimensions from original file
for name, dimension in file_pr.dimensions.items():
    # if any(name in s for s in ['lon']):
    #     continue
    file_prfrac.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

# copy all variables except time from ps file
for name, variable in file_pr.variables.items():
    if any(name in s for s in ['precip']):
        continue
    
    x = file_prfrac.createVariable(name, variable.datatype, variable.dimensions)
    file_prfrac[name].setncatts(file_pr[name].__dict__)
    file_prfrac[name][:] = file_pr[name][:]

precip = file_pr.variables['precip']
pr_frac = file_prfrac.createVariable('prfrac', 'f4', precip.dimensions)
pr_frac.units = 'unitless'
pr_frac.long_name = 'convective precipitation fraction'
pr_frac[:] = prfrac
file_prfrac.close()

