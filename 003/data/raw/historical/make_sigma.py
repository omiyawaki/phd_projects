import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

# for CMIP

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps and mse files
path_ps = sys.argv[1]
path_tas = sys.argv[2]
path_temp = sys.argv[3]
path_tempsi = sys.argv[4]

# open files
file_ps = Dataset(path_ps, 'r')
file_tas = Dataset(path_tas, 'r')
file_temp = Dataset(path_temp, 'r')

# read data
ps = file_ps.variables['ps'][:] # (mon x lat x lon)
tas = file_tas.variables['tas'][:] # (mon x lat x lon)
temp = file_temp.variables['ta'][:] # (mon x lev x lat x lon)
plev = file_temp.variables['plev'][:] # (lev)

# # for datasets that fill data below surface as missing data, fill with nans
plev = plev.filled(fill_value=np.nan)
si = 1e-5*plev
ps = ps.filled(fill_value=np.nan)
tas = tas.filled(fill_value=np.nan)
temp = temp.filled(fill_value=np.nan)

# mask out data below surface
pa = np.transpose(np.tile(plev, [ps.shape[0], ps.shape[1], ps.shape[2], 1]), [0,3,1,2])
subsurf = pa > np.transpose(np.tile(ps, [len(plev),1,1,1]),[1,0,2,3])
temp[subsurf] = np.nan
pa[subsurf] = np.nan

# interpolate to sigma coord
tempsi = np.empty_like(temp)
for itime in tqdm(range(temp.shape[0])):
    for ilat in range(temp.shape[2]):
        for ilon in range(temp.shape[3]):
            ps_local = ps[itime, ilat, ilon]
            tas_local = tas[itime, ilat, ilon]

            plev_local = pa[itime, :, ilat, ilon]
            temp_local = temp[itime, :, ilat, ilon]

            plev_local = plev_local[~np.isnan(plev_local)]
            temp_local = temp_local[~np.isnan(temp_local)]

            if plev[1]-plev[0]>0: # if pressure increases with index
                if not plev_local[-1] == ps_local:
                    plev_local = np.append(plev_local, ps_local) 
                    temp_local = np.append(temp_local, tas_local)
            else:
                if not plev_local[0] == ps_local:
                    plev_local = np.insert(plev_local, 0, ps_local) 
                    temp_local = np.insert(temp_local, 0, tas_local)

            si_local = plev_local / ps_local

            f_temp = interpolate.interp1d(si_local, temp_local, kind='linear', fill_value="extrapolate")
            tempsi[itime, :, ilat, ilon] = f_temp(si)

# save file as netCDF
file_tempsi = Dataset(path_tempsi, "w", format='NETCDF4_CLASSIC')

# copy dimensions from temp file
for name, dimension in file_temp.dimensions.items():
    file_tempsi.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy variables from temp file
for name, variable in file_temp.variables.items():
    if any(name in s for s in ['ta']):
        continue

    x = file_tempsi.createVariable(name, variable.datatype, variable.dimensions)
    file_tempsi[name].setncatts(file_temp[name].__dict__)
    file_tempsi[name][:] = file_temp[name][:]

temp_si = file_tempsi.createVariable('tempsi', 'f4', ("time","plev","lat","lon"))
temp_si.units = "K"
temp_si.long_name = "temperature"
temp_si[:,:,:,:] = tempsi

file_tempsi.close()
