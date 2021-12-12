import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from netCDF4 import Dataset

##############################
# TRANSIENT EDDIES (6H ML and DM PL)
##############################

divte_ml_file = Dataset('./divte_6h_ml_echr0001_103001.nc', 'r')
divte_ml = np.squeeze(np.nanmean(divte_ml_file.variables['divte'][:],axis=0))
lat = divte_ml_file.variables['lat'][:]

divte_pl_file = Dataset('./divte_dm_pl_echr0001_103001.nc', 'r')
divte_pl = np.squeeze(np.nanmean(divte_pl_file.variables['divte'][:],axis=0))

plotname = './divte_ml6h_pldm'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
lp_divte_ml = ax.plot(lat, divte_ml, color='r', label='ML (6h)')
lp_divte_pl = ax.plot(lat, divte_pl, ':', color='r', label='PL (daily)')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('$\partial_y(\overline{[v^{*\prime} m^{*\prime}]})$')
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('MSE flux divergence due to transient eddies (W m$^{-2}$)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlim([-90,90])
# ax.set_ylim(vmin_dev,vmax_dev)
ax.legend()
plt.tight_layout()
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
plt.show()
plt.close()

##############################
# TRANSIENT EDDIES (6H ML and 6H PL)
##############################

divte_ml_file = Dataset('./divte_6h_ml_echr0001_103001.nc', 'r')
divte_ml = np.squeeze(np.nanmean(divte_ml_file.variables['divte'][:],axis=0))
lat = divte_ml_file.variables['lat'][:]

divte_pl_file = Dataset('./divte_6h_pl_echr0001_103001.nc', 'r')
divte_pl = np.squeeze(np.nanmean(divte_pl_file.variables['divte'][:],axis=0))

plotname = './divte_ml6h_pl6h'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
lp_divte_ml = ax.plot(lat, divte_ml, color='r', label='ML (6h)')
lp_divte_pl = ax.plot(lat, divte_pl, ':', color='r', label='PL (6h)')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('$\partial_y(\overline{[v^{*\prime} m^{*\prime}]})$')
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('MSE flux divergence due to transient eddies (W m$^{-2}$)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlim([-90,90])
# ax.set_ylim(vmin_dev,vmax_dev)
ax.legend()
plt.tight_layout()
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
plt.show()
plt.close()
