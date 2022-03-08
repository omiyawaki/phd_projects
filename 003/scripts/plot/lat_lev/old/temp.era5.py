import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from netCDF4 import Dataset

##############################################
# ANNUAL MEAN TRENDS
##############################################

file = Dataset('/project2/tas1/miyawaki/projects/002/data/raw/era5c/temp/era5c_tempsi_1980_2005.zonmean.yearmean.nc', 'r')
temp = np.squeeze(file.variables['tempsi'][:])
lev = 1e-3*np.squeeze(file.variables['level'][:])
lat = np.squeeze(file.variables['lat'][:])

time = range(temp.shape[0])

# compute trends
m = np.empty([temp.shape[1], temp.shape[2]])
for ilev in range(temp.shape[1]):
    for ilat in range(temp.shape[2]):
        A = np.vstack([time, np.ones(len(time))]).T
        m[ilev, ilat], c = np.linalg.lstsq(A, temp[:,ilev,ilat], rcond=None)[0]
        
[mesh_lev, mesh_lat] = np.meshgrid(lev, lat) # create mesh
plotname = 'ann.era5'
fig, ax = plt.subplots()
csf = ax.contourf(mesh_lat, mesh_lev, 10*np.transpose(m), np.arange(-0.5,0.5,0.02), cmap='RdBu_r', extend='both')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('ERA5, ANN')
ax.set_xlabel('Latitude (deg)')
ax.set_xticks(np.arange(-90,91,30))
ax.set_ylabel('$\sigma$ (unitless)')
ax.set_yticks(np.arange(0,1.1,0.1))
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(MultipleLocator(10))
cbar = plt.colorbar(csf)
cbar.set_label('$T$ trend (K decade$^{-1}$)')
# plt.savefig('%s.png' % (plotname), dpi=300)
plt.gca().invert_yaxis()
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
plt.show()
plt.close()

##############################################
# JANUARY TRENDS
##############################################

file = Dataset('/project2/tas1/miyawaki/projects/002/data/raw/era5c/temp/era5c_tempsi_1980_2005.zonmean.1.nc', 'r')
temp = np.squeeze(file.variables['tempsi'][:])
lev = 1e-3*np.squeeze(file.variables['level'][:])
lat = np.squeeze(file.variables['lat'][:])

time = range(temp.shape[0])

# compute trends
m = np.empty([temp.shape[1], temp.shape[2]])
for ilev in range(temp.shape[1]):
    for ilat in range(temp.shape[2]):
        A = np.vstack([time, np.ones(len(time))]).T
        m[ilev, ilat], c = np.linalg.lstsq(A, temp[:,ilev,ilat], rcond=None)[0]
        
[mesh_lev, mesh_lat] = np.meshgrid(lev, lat) # create mesh
plotname = 'jan.era5'
fig, ax = plt.subplots()
csf = ax.contourf(mesh_lat, mesh_lev, 10*np.transpose(m), np.arange(-0.5,0.5,0.02), cmap='RdBu_r', extend='both')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('ERA5, January')
ax.set_xlabel('Latitude (deg)')
ax.set_xticks(np.arange(-90,91,30))
ax.set_ylabel('$\sigma$ (unitless)')
ax.set_yticks(np.arange(0,1.1,0.1))
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(MultipleLocator(10))
cbar = plt.colorbar(csf)
cbar.set_label('$T$ trend (K decade$^{-1}$)')
# plt.savefig('%s.png' % (plotname), dpi=300)
plt.gca().invert_yaxis()
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
plt.show()
plt.close()

##############################################
# JUNE TRENDS
##############################################

file = Dataset('/project2/tas1/miyawaki/projects/002/data/raw/era5c/temp/era5c_tempsi_1980_2005.zonmean.6.nc', 'r')
temp = np.squeeze(file.variables['tempsi'][:])
lev = 1e-3*np.squeeze(file.variables['level'][:])
lat = np.squeeze(file.variables['lat'][:])

time = range(temp.shape[0])

# compute trends
m = np.empty([temp.shape[1], temp.shape[2]])
for ilev in range(temp.shape[1]):
    for ilat in range(temp.shape[2]):
        A = np.vstack([time, np.ones(len(time))]).T
        m[ilev, ilat], c = np.linalg.lstsq(A, temp[:,ilev,ilat], rcond=None)[0]
        
[mesh_lev, mesh_lat] = np.meshgrid(lev, lat) # create mesh
plotname = 'jun.era5'
fig, ax = plt.subplots()
csf = ax.contourf(mesh_lat, mesh_lev, 10*np.transpose(m), np.arange(-0.5,0.5,0.02), cmap='RdBu_r', extend='both')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('ERA5, June')
ax.set_xlabel('Latitude (deg)')
ax.set_xticks(np.arange(-90,91,30))
ax.set_ylabel('$\sigma$ (unitless)')
ax.set_yticks(np.arange(0,1.1,0.1))
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(MultipleLocator(10))
cbar = plt.colorbar(csf)
cbar.set_label('$T$ trend (K decade$^{-1}$)')
# plt.savefig('%s.png' % (plotname), dpi=300)
plt.gca().invert_yaxis()
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
plt.show()
plt.close()