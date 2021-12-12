import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from netCDF4 import Dataset

##############################
# MMC
##############################

vmmmc_ml_file = Dataset('./vmmmc_6h_ml_echr0001_103001.nc', 'r')
vmmmc_ml = np.squeeze(np.mean(vmmmc_ml_file.variables['vmmmc'][:],axis=0))
lat = vmmmc_ml_file.variables['lat'][:]

vmmmc_pl_file = Dataset('./vmmmc_6h_pl_echr0001_103001.nc', 'r')
vmmmc_pl = np.squeeze(np.mean(vmmmc_pl_file.variables['vmmmc'][:],axis=0))

plotname = './vmmmc'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
lp_vmmmc_ml = ax.plot(lat, 1e-15*vmmmc_ml, color='g', label='ML')
lp_vmmmc_pl = ax.plot(lat, 1e-15*vmmmc_pl, ':', color='g', label='PL')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('MMC')
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
# STATIONARY EDDIES
##############################

vmse_ml_file = Dataset('./vmse_6h_ml_echr0001_103001.nc', 'r')
vmse_ml = np.squeeze(np.mean(vmse_ml_file.variables['vmse'][:],axis=0))

vmse_pl_file = Dataset('./vmse_6h_pl_echr0001_103001.nc', 'r')
vmse_pl = np.squeeze(np.mean(vmse_pl_file.variables['vmse'][:],axis=0))

plotname = './vmse'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
lp_vmse_ml = ax.plot(lat, 1e-15*vmse_ml, color='b', label='ML')
lp_vmse_pl = ax.plot(lat, 1e-15*vmse_pl, ':', color='b', label='PL')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('SE')
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

###############################
## TRANSIENT EDDIES
###############################

#vmte_ml_file = Dataset('./vmte_6h_ml_echr0001_103001.nc', 'r')
#vmte_ml = np.squeeze(np.mean(vmte_ml_file.variables['vmte'][:],axis=0))

#vmte_pl_file = Dataset('./vmte_6h_pl_echr0001_103001.nc', 'r')
#vmte_pl = np.squeeze(np.mean(vmte_pl_file.variables['vmte'][:],axis=0))

#plotname = './vmte'
#fig, ax = plt.subplots()
#ax.axhline(0, color='k', linewidth=0.5)
#lp_vmte_ml = ax.plot(lat, 1e-15*vmte_ml, color='r', label='ML')
#lp_vmte_pl = ax.plot(lat, 1e-15*vmte_pl, ':', color='r', label='PL')
#ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
#ax.set_title('TE')
#ax.set_xlabel('Latitude (deg)')
#ax.set_ylabel('Energy transport (PW)')
#ax.xaxis.set_minor_locator(MultipleLocator(10))
#ax.yaxis.set_minor_locator(AutoMinorLocator())
#ax.set_xlim([-90,90])
## ax.set_ylim(vmin_dev,vmax_dev)
#ax.legend()
#plt.tight_layout()
#plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
#plt.show()
#plt.close()
