import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from netCDF4 import Dataset

model = "bcc-csm1-1"

file_rcp85 = Dataset('/project2/tas1/miyawaki/projects/003/data/raw/rcp85/%s/wap_Amon_%s_rcp85_r1i1p1_200601-229912.ymonmean-30.nc' % (model, model), 'r')
wap_rcp85 = np.squeeze(np.nanmean(file_rcp85.variables['wap'][:], axis=3))
lev = np.squeeze(file_rcp85.variables['plev'][:])
lat = np.squeeze(file_rcp85.variables['lat'][:])

file_hist = Dataset('/project2/tas1/miyawaki/projects/003/data/raw/historical/%s/wap_Amon_%s_historical_r1i1p1_186001-200512.ymonmean-30.nc' % (model, model), 'r')
wap_hist = np.squeeze(np.nanmean(file_hist.variables['wap'][:], axis=3))

# convert to hPa / day
wap_rcp85 = 1e-2* 86400 * wap_rcp85
wap_hist = 1e-2* 86400 * wap_hist

# contours
cntrs_clim = np.arange(-50,50,1)
cntrs_delt = np.arange(-5,5,0.1)

##############################################
# ANNUAL MEAN CLIMATOLOGY
##############################################
[mesh_lev, mesh_lat] = np.meshgrid(lev, lat) # create mesh
plotname = 'ann.%s' % (model)
fig, ax = plt.subplots()
csf = ax.contourf(mesh_lat, mesh_lev, np.transpose(np.nanmean(wap_hist, axis=0)), cntrs_clim, cmap='RdBu_r', extend='both')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('%s, historical, ANN' % (model))
ax.set_xlabel('Latitude (deg)')
ax.set_xticks(np.arange(-90,91,30))
ax.set_ylabel('$\sigma$ (unitless)')
# ax.set_yticks(1e5*np.arange(0,1.1,0.1))
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
cbar = plt.colorbar(csf)
cbar.set_label('$\Delta\omega$ (Pa s$^{-1}$)')
# plt.savefig('%s.png' % (plotname), dpi=300)
plt.gca().invert_yaxis()
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
# plt.show()
plt.close()

##############################################
# ANNUAL MEAN TRENDS
##############################################

# compute annual mean difference
dwap_ann = np.nanmean(wap_rcp85, axis=0) - np.nanmean(wap_hist, axis=0)
        
[mesh_lev, mesh_lat] = np.meshgrid(lev, lat) # create mesh
plotname = 'dann.%s' % (model)
fig, ax = plt.subplots()
csf = ax.contourf(mesh_lat, mesh_lev, np.transpose(dwap_ann), cntrs_delt, cmap='RdBu_r', extend='both')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('%s, RCP8.5$-$historical, ANN' % (model))
ax.set_xlabel('Latitude (deg)')
ax.set_xticks(np.arange(-90,91,30))
ax.set_ylabel('$\sigma$ (unitless)')
# ax.set_yticks(1e5*np.arange(0,1.1,0.1))
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
cbar = plt.colorbar(csf)
cbar.set_label('$\Delta\omega$ (Pa s$^{-1}$)')
# plt.savefig('%s.png' % (plotname), dpi=300)
plt.gca().invert_yaxis()
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
# plt.show()
plt.close()

##############################################
# JANUARY TRENDS
##############################################

# compute annual mean difference
dwap_jan = wap_rcp85[0,...] - wap_hist[0,...]
        
[mesh_lev, mesh_lat] = np.meshgrid(lev, lat) # create mesh
plotname = 'djan.%s' % (model)
fig, ax = plt.subplots()
csf = ax.contourf(mesh_lat, mesh_lev, np.transpose(dwap_jan), cntrs_delt, cmap='RdBu_r', extend='both')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('%s, RCP8.5$-$historical, January' % (model))
ax.set_xlabel('Latitude (deg)')
ax.set_xticks(np.arange(-90,91,30))
ax.set_ylabel('$\sigma$ (unitless)')
# ax.set_yticks(1e5*np.arange(0,1.1,0.1))
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
cbar = plt.colorbar(csf)
cbar.set_label('$\Delta\omega$ (Pa s$^{-1}$)')
# plt.savefig('%s.png' % (plotname), dpi=300)
plt.gca().invert_yaxis()
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
# plt.show()
plt.close()

