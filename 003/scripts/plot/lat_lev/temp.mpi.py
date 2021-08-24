import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from netCDF4 import Dataset

file_rcp85 = Dataset('/project2/tas1/miyawaki/projects/003/data/raw/rcp85/MPI-ESM-LR/tasi_Amon_MPI-ESM-LR_rcp85_r1i1p1_200601-210012.ymonmean-30.nc', 'r')
temp_rcp85 = np.squeeze(np.nanmean(file_rcp85.variables['tempsi'][:], axis=3))
lev = 1e-5*np.squeeze(file_rcp85.variables['plev'][:])
lat = np.squeeze(file_rcp85.variables['lat'][:])

file_hist = Dataset('/project2/tas1/miyawaki/projects/003/data/raw/historical/MPI-ESM-LR/tasi_Amon_MPI-ESM-LR_historical_r1i1p1_186001-200512.ymonmean-30.nc', 'r')
temp_hist = np.squeeze(np.nanmean(file_hist.variables['tempsi'][:], axis=3))

##############################################
# ANNUAL MEAN TRENDS
##############################################

# compute annual mean difference
dtemp_ann = np.nanmean(temp_rcp85, axis=0) - np.nanmean(temp_hist, axis=0)
        
[mesh_lev, mesh_lat] = np.meshgrid(lev, lat) # create mesh
plotname = 'ann.mpi'
fig, ax = plt.subplots()
csf = ax.contourf(mesh_lat, mesh_lev, np.transpose(dtemp_ann), np.arange(-10,10,1), cmap='RdBu_r', extend='both')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('MPI-ESM-LR, RCP8.5$-$historical, ANN')
ax.set_xlabel('Latitude (deg)')
ax.set_xticks(np.arange(-90,91,30))
ax.set_ylabel('$\sigma$ (unitless)')
ax.set_yticks(np.arange(0,1.1,0.1))
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(MultipleLocator(10))
cbar = plt.colorbar(csf)
cbar.set_label('$\Delta T$ (K)')
# plt.savefig('%s.png' % (plotname), dpi=300)
plt.gca().invert_yaxis()
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
# plt.show()
plt.close()

##############################################
# JANUARY TRENDS
##############################################

# compute annual mean difference
dtemp_jan = temp_rcp85[0,...] - temp_hist[0,...]
        
[mesh_lev, mesh_lat] = np.meshgrid(lev, lat) # create mesh
plotname = 'jan.mpi'
fig, ax = plt.subplots()
csf = ax.contourf(mesh_lat, mesh_lev, np.transpose(dtemp_jan), np.arange(-10,10,0.5), cmap='RdBu_r', extend='both')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('MPI-ESM-LR, RCP8.5$-$historical, January')
ax.set_xlabel('Latitude (deg)')
ax.set_xticks(np.arange(-90,91,30))
ax.set_ylabel('$\sigma$ (unitless)')
ax.set_yticks(np.arange(0,1.1,0.1))
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(MultipleLocator(10))
cbar = plt.colorbar(csf)
cbar.set_label('$\Delta T$ (K)')
# plt.savefig('%s.png' % (plotname), dpi=300)
plt.gca().invert_yaxis()
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
# plt.show()
plt.close()

##############################################
# JUNE TRENDS
##############################################

# compute annual mean difference
dtemp_jun = temp_rcp85[5,...] - temp_hist[5,...]
        
[mesh_lev, mesh_lat] = np.meshgrid(lev, lat) # create mesh
plotname = 'jun.mpi'
fig, ax = plt.subplots()
csf = ax.contourf(mesh_lat, mesh_lev, np.transpose(dtemp_jun), np.arange(-10,10,0.5), cmap='RdBu_r', extend='both')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('MPI-ESM-LR, RCP8.5$-$historical, June')
ax.set_xlabel('Latitude (deg)')
ax.set_xticks(np.arange(-90,91,30))
ax.set_ylabel('$\sigma$ (unitless)')
ax.set_yticks(np.arange(0,1.1,0.1))
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(MultipleLocator(10))
cbar = plt.colorbar(csf)
cbar.set_label('$\Delta T$ (K)')
# plt.savefig('%s.png' % (plotname), dpi=300)
plt.gca().invert_yaxis()
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
# plt.show()
plt.close()

##############################################
# JULY TRENDS
##############################################

# compute annual mean difference
dtemp_jul = temp_rcp85[6,...] - temp_hist[6,...]
        
[mesh_lev, mesh_lat] = np.meshgrid(lev, lat) # create mesh
plotname = 'jul.mpi'
fig, ax = plt.subplots()
csf = ax.contourf(mesh_lat, mesh_lev, np.transpose(dtemp_jul), np.arange(-10,10,0.5), cmap='RdBu_r', extend='both')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('MPI-ESM-LR, RCP8.5$-$historical, July')
ax.set_xlabel('Latitude (deg)')
ax.set_xticks(np.arange(-90,91,30))
ax.set_ylabel('$\sigma$ (unitless)')
ax.set_yticks(np.arange(0,1.1,0.1))
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(MultipleLocator(10))
cbar = plt.colorbar(csf)
cbar.set_label('$\Delta T$ (K)')
# plt.savefig('%s.png' % (plotname), dpi=300)
plt.gca().invert_yaxis()
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
# plt.show()
plt.close()

##############################################
# DECEMBER TRENDS
##############################################

# compute annual mean difference
dtemp_dec = temp_rcp85[-1,...] - temp_hist[-1,...]
        
[mesh_lev, mesh_lat] = np.meshgrid(lev, lat) # create mesh
plotname = 'dec.mpi'
fig, ax = plt.subplots()
csf = ax.contourf(mesh_lat, mesh_lev, np.transpose(dtemp_dec), np.arange(-10,10,0.5), cmap='RdBu_r', extend='both')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('MPI-ESM-LR, RCP8.5$-$historical, December')
ax.set_xlabel('Latitude (deg)')
ax.set_xticks(np.arange(-90,91,30))
ax.set_ylabel('$\sigma$ (unitless)')
ax.set_yticks(np.arange(0,1.1,0.1))
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(MultipleLocator(10))
cbar = plt.colorbar(csf)
cbar.set_label('$\Delta T$ (K)')
# plt.savefig('%s.png' % (plotname), dpi=300)
plt.gca().invert_yaxis()
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
plt.show()
plt.close()
