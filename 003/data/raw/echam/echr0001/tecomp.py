import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from netCDF4 import Dataset

##############################
# TRANSIENT EDDIES (6H ML and DM PL)
##############################

vmte_ml_file = Dataset('./vmte_6h_ml_echr0001_103001.nc', 'r')
vmte_ml = np.squeeze(np.mean(vmte_ml_file.variables['vmte'][:],axis=0))
lat = vmte_ml_file.variables['lat'][:]

vmte_pl_file = Dataset('./vmte_dm_pl_echr0001_103001.nc', 'r')
vmte_pl = np.squeeze(np.mean(vmte_pl_file.variables['vmte'][:],axis=0))

plotname = './vmte_ml6h_pldm'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
lp_vmte_ml = ax.plot(lat, 1e-15*vmte_ml, color='r', label='ML (6h)')
lp_vmte_pl = ax.plot(lat, 1e-15*vmte_pl, ':', color='r', label='PL (daily)')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('$\overline{[v^{*\prime} m^{*\prime}]}$')
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('Energy transport (PW)')
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

vmte_ml_file = Dataset('./vmte_6h_ml_echr0001_103001.nc', 'r')
vmte_ml = np.squeeze(np.mean(vmte_ml_file.variables['vmte'][:],axis=0))
lat = vmte_ml_file.variables['lat'][:]

vmte_pl_file = Dataset('./vmte_6h_pl_echr0001_103001.nc', 'r')
vmte_pl = np.squeeze(np.mean(vmte_pl_file.variables['vmte'][:],axis=0))

plotname = './vmte_ml6h_pl6h'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
lp_vmte_ml = ax.plot(lat, 1e-15*vmte_ml, color='r', label='ML (6h)')
lp_vmte_pl = ax.plot(lat, 1e-15*vmte_pl, ':', color='r', label='PL (6h)')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('$\overline{[v^{*\prime} m^{*\prime}]}$')
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('Energy transport (PW)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlim([-90,90])
# ax.set_ylim(vmin_dev,vmax_dev)
ax.legend()
plt.tight_layout()
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
plt.show()
plt.close()
