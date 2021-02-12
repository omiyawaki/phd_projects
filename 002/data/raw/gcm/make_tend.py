import sys
import numpy as np
from scipy import interpolate
import datetime
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

# try to open huss (not all models have huss)
try:
    file_huss = Dataset(path_huss, 'r')
except: # if huss is not available, compute using monthly hus
    path_husmon = path_huss.replace('huss', 'hus')
    file_husmon = Dataset(path_husmon, 'r')
    husmon = file_husmon.variables['husmon'][:] # (mon x lev x lat x lon)
    filledhusmon = husmon.filled(fill_value=np.nan)
    huss = np.empty(ps.shape)
    for itime in range(0,filledhusmon.shape[0]):
        for ilat in range(0,filledhusmon.shape[2]):
            for ilon in range(0,filledhusmon.shape[3]):
                f = interpolate.interp1d(plev, np.squeeze(filledhusmon[itime,:,ilat,ilon]))
                huss[itime,ilat,ilon] = f(ps[itime,ilat,ilon])
    husmon = None; filledhusmon = None;
else:
    huss = file_huss.variables['huss'][:] # (mon x lev x lat x lon)
    if not np.array_equal(lat2d,lat3d):
        filledhuss = huss.filled(fill_value=np.nan)
        f = interpolate.interp1d(filledlat2d, filledhuss, axis=1)
        huss = f(filledlat3d)
    filledhuss = None; filledlat2d = None; filledlat3d = None; 
    
monthlytime = file_ps.variables['time'] # (day x lev x lat x lon)
monthlydate = num2date(monthlytime[:], units=monthlytime.units, calendar=monthlytime.calendar)
dailytime = file_ta.variables['time'] # (day x lev x lat x lon)
dailydate = num2date(dailytime[:], units=dailytime.units, calendar=dailytime.calendar)

# calculate MSE
mse = cpd*ta + g*zg + L*hus
ta = None; zg = None; hus = None; # collect garbage

# resample monthly data at daily frequency
monthly = np.empty([ps.shape[0],2])
for i in range(0,ps.shape[0]):
    monthly[i,0] = monthlydate[i].year
    monthly[i,1] = monthlydate[i].month

daily = np.empty([1,2])
psd = np.empty([mse.shape[0],ps.shape[1],ps.shape[2]])
tasd = np.empty([mse.shape[0],ps.shape[1],ps.shape[2]])
hussd = np.empty([mse.shape[0],ps.shape[1],ps.shape[2]])
for i in range(0,mse.shape[0]):
    daily[0,0] = dailydate[i].year
    daily[0,1] = dailydate[i].month
    idx_ps = np.where((monthly==daily).all(axis=1))[0][0]
    psd[i,:,:] = ps[idx_ps,:,:]
    tasd[i,:,:] = tas[idx_ps,:,:]
    hussd[i,:,:] = huss[idx_ps,:,:]

ps = None; tas = None; huss = None;

# mask out data below surface
#plev = file_ta.variables['plev'][:] 
#plev4d = np.tile(plev[np.newaxis,:,np.newaxis,np.newaxis], (mse.shape[0],1,mse.shape[2],mse.shape[3]))
#psd4d = np.tile(psd[:,np.newaxis,:,:], (1,plev.size,1,1))
#mse[plev4d > psd4d] = 0
#psd4d = None;

# compute vertical integral
vmse = np.empty(psd.shape)
# track progress of loop using a progress bar
toolbar_width = 40
sys.stdout.write("[%s]" % (" " * toolbar_width))
sys.stdout.flush()
sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['

for itime in range(0,mse.shape[0]):
    for ilat in range(0,mse.shape[2]):
        for ilon in range(0,mse.shape[3]):
            # levels below surface
            psdcol = psd[itime,ilat,ilon]
            idx_below = np.where(plev > psdcol)
            
            msecol = np.squeeze(mse[itime,:,ilat,ilon])
            
            tascol = np.squeeze(tasd[itime,ilat,ilon])
            husscol = np.squeeze(hussd[itime,ilat,ilon])
            orogcol = np.squeeze(orog[ilat,ilon])
            mse0 = cpd*tascol + g*orogcol + L*husscol
            
            plevmsk = np.delete(plev, idx_below)
            msemsk = np.delete(msecol, idx_below)
            
            plevmsk = np.insert(plevmsk, 0, psdcol)
            msemsk = np.insert(msemsk, 0, mse0)
            
            if plev[1]-plev[0]>0: # if pressure increases with index
                vmse[itime,ilat,ilon] = 1/g*np.trapz(msemsk,plevmsk)
            else:
                vmse[itime,ilat,ilon] = -1/g*np.trapz(msemsk,plevmsk)
                
    # update the bar
    sys.stdout.write("-")
    sys.stdout.flush()
    
sys.stdout.write("]\n") # this ends the progress bar
    
#if plev[1]-plev[0]>0: # if pressure increases with index
#    vmse = 1/g*np.trapz(mse,plev4d,axis=1)
#else:
#    vmse = -1/g*np.trapz(mse,plev4d,axis=1)
    
# take time tendency
dvmsedt = np.empty(vmse.shape) 
dvmsedt[1:-1,:,:] = (vmse[2:,:,:]-vmse[0:-2,:,:])/(2*86400)
dvmsedt[0] = (vmse[1,:,:]-vmse[0,:,:])/86400
dvmsedt[-1] = (vmse[-1,:,:]-vmse[-2,:,:])/86400

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
