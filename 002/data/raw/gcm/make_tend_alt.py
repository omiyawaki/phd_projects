import sys
import numpy as np
from scipy import interpolate
import datetime
import time
from tqdm import tqdm
from netCDF4 import Dataset,num2date

cpd = 1005.7; Rd = 287; Rv = 461; L = 2.501e6; g = 9.81; a = 6357e3; eps = Rd/Rv;

# paths to ps, ta, zg, and hus files
path_ps = sys.argv[1]
path_ta = sys.argv[2]
path_zg = sys.argv[3]
path_hus = sys.argv[4]
path_tend = sys.argv[5]
path_tas = sys.argv[6]
path_orog = sys.argv[7]
path_huss = sys.argv[8]

# open files
file_ps = Dataset(path_ps, 'r')
file_ta = Dataset(path_ta, 'r')
file_zg = Dataset(path_zg, 'r')
file_hus = Dataset(path_hus, 'r')
file_tas = Dataset(path_tas, 'r')
file_orog = Dataset(path_orog, 'r')

# read data
orog = file_orog.variables['orog'][:] # (lat x lon)
ps = file_ps.variables['ps'][:] # (mon x lat x lon)
tas = file_tas.variables['tas'][:] # (mon x lev x lat x lon)
ta = file_ta.variables['ta'][:] # (day x lev x lat x lon)
zg = file_zg.variables['zg'][:] # (day x lev x lat x lon)
hus = file_hus.variables['hus'][:] # (day x lev x lat x lon)
plev = file_ta.variables['plev'][:] 

# check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
lat2d = file_ps.variables['lat'][:] # (mon x lat x lon)
lat3d = file_ta.variables['lat'][:] # (day x lev x lat x lon)
if not np.array_equal(lat2d,lat3d):
    filledps = ps.filled(fill_value=np.nan)
    filledtas = tas.filled(fill_value=np.nan)
    filledorog = orog.filled(fill_value=np.nan)
    filledlat2d = lat2d.filled(fill_value=np.nan)
    filledlat3d = lat3d.filled(fill_value=np.nan)
    f = interpolate.interp1d(filledlat2d, filledps, axis=1)
    ps = f(filledlat3d)
    f = interpolate.interp1d(filledlat2d, filledtas, axis=1)
    tas = f(filledlat3d)
    f = interpolate.interp1d(filledlat2d, filledorog, axis=0)
    orog = f(filledlat3d)
    filledps = None; filledtas = None; filledorog = None; f = None;

# resample monthly data at daily frequency
monthlytime = file_ps.variables['time'] # (day x lev x lat x lon)
monthlydate = num2date(monthlytime[:], units=monthlytime.units, calendar=monthlytime.calendar)
dailytime = file_ta.variables['time'] # (day x lev x lat x lon)
dailydate = num2date(dailytime[:], units=dailytime.units, calendar=dailytime.calendar)

monthly = np.empty([ps.shape[0],2])
print('Resampling monthly data at daily frequency...');
for i in tqdm(range(0,ps.shape[0])):
    monthly[i,0] = monthlydate[i].year
    monthly[i,1] = monthlydate[i].month

daily = np.empty([1,2])
psd = np.empty([ta.shape[0],ps.shape[1],ps.shape[2]])
for i in tqdm(range(0,ta.shape[0])):
    daily[0,0] = dailydate[i].year
    daily[0,1] = dailydate[i].month
    idx_ps = np.where((monthly==daily).all(axis=1))[0][0]
    psd[i,:,:] = ps[idx_ps,:,:]

# for datasets that fill data below surface as missing data, fill with nans
ta = ta.filled(fill_value=np.nan)
zg = zg.filled(fill_value=np.nan)
hus = hus.filled(fill_value=np.nan)
# replace subsurface data in general with nans (important for computing dhdt near the surface)
ps3d = np.tile(psd, [plev.size, 1, 1, 1])
ps3d = np.transpose( ps3d, [1, 0, 2, 3] )
pa3d = np.tile(plev, [psd.shape[0], psd.shape[1], psd.shape[2], 1])
pa3d = np.transpose( pa3d, [0, 3, 1, 2] )
idx_below = pa3d > ps3d
ps3d = None;

pa3d[idx_below]=np.nan
ta[idx_below]=np.nan
zg[idx_below]=np.nan
hus[idx_below]=np.nan

# calculate MSE
mse = cpd*ta + g*zg + L*hus
ta = None; zg = None; hus = None; # collect garbage

# take time tendency
# 3d mse tendency
dmsedt = np.empty(mse.shape) 
dmsedt[1:-1,:,:,:] = (mse[2:,:,:,:]-mse[0:-2,:,:,:])/(2*86400)
dmsedt[0] = (mse[1,:,:,:]-mse[0,:,:,:])/86400
dmsedt[-1] = (mse[-1,:,:,:]-mse[-2,:,:,:])/86400

# compute vertical integral
if plev[1]-plev[0]>0: # if pressure increases with index
    dvmsedt = 1/g*np.trapz(dmsedt, pa3d, axis=1)
else:
    dvmsedt = -1/g*np.trapz(dmsedt, pa3d, axis=1)

# save file as netCDF
file_tend = Dataset(path_tend, "w", format='NETCDF4_CLASSIC')

# copy attributes from ps file
file_tend.setncatts(file_ps.__dict__)

# copy dimensions from ta file
for name, dimension in file_ta.dimensions.items():
    if any(name in s for s in ['plev']):
        continue
    file_tend.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
 
# copy all variables except time from ps file
for name, variable in file_ta.variables.items():
    if any(name in s for s in ['ta' 'plev' 'plev_bnds']):
        continue
    
    x = file_tend.createVariable(name, variable.datatype, variable.dimensions)
    file_tend[name].setncatts(file_ta[name].__dict__)
    file_tend[name][:] = file_ta[name][:]
    
tend = file_tend.createVariable('tend', 'f4', ("time","lat","lon"))
tend.units = "W/m^2"
tend.long_name = "time tendency of vertically integrated moist static energy"
tend[:,:,:] = dvmsedt

file_tend.close()
