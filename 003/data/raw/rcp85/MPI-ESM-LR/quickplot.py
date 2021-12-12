import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from netCDF4 import Dataset

##############################
# plot all components
##############################


vmmmc_file = Dataset('./vmmmc_Amon_MPI-ESM-LR_rcp85_r1i1p1_200601-200612.nc', 'r')
vmmmc = vmmmc_file.variables['vmmmc'][:][0,...]
# vmmmc = np.squeeze(np.mean(vmmmc_file.variables['vmmmc'][:],axis=0))

vmse_file = Dataset('./vmse_Amon_MPI-ESM-LR_rcp85_r1i1p1_200601-200612.nc', 'r')
vmse = vmse_file.variables['vmse'][:][0,...]
# vmse = np.squeeze(np.mean(vmse_file.variables['vmse'][:],axis=0))

lat = vmmmc_file.variables['lat'][:]

plotname = './trans_test'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
lp_vmmmc = ax.plot(lat, 1e-15*vmmmc, color='b', label='MMC')
lp_vmse = ax.plot(lat, 1e-15*vmse, color='g', label='SE')
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
