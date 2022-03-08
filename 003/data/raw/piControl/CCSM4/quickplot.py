import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from netCDF4 import Dataset

##############################
# plot all components
##############################

aht_file = Dataset('./aht_Amon_CCSM4_piControl_r1i1p1_125001-129912.timmean.nc', 'r')
aht = aht_file.variables['aht'][:][0,...]
# aht = np.squeeze(np.mean(aht_file.variables['aht'][:],axis=0))

vmmmc_file = Dataset('./vmmmc_Amon_CCSM4_piControl_r1i1p1_125001-129912.timmean.nc', 'r')
vmmmc = vmmmc_file.variables['vmmmc'][:][0,...]
# vmmmc = np.squeeze(np.mean(vmmmc_file.variables['vmmmc'][:],axis=0))

vmse_file = Dataset('./vmse_Amon_CCSM4_piControl_r1i1p1_125001-129912.timmean.nc', 'r')
vmse = vmse_file.variables['vmse'][:][0,...]
# vmse = np.squeeze(np.mean(vmse_file.variables['vmse'][:],axis=0))

vmte_file = Dataset('./vmte_Amon_CCSM4_piControl_r1i1p1_125001-129912.timmean.nc', 'r')
vmte = vmte_file.variables['vmte'][:][0,...]
# vmte = np.squeeze(np.mean(vmte_file.variables['vmte'][:],axis=0))

lat = vmmmc_file.variables['lat'][:]

plotname = './aaron_comp'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
lp_aht = ax.plot(lat, 1e-15*aht, color='k', label='AHT')
lp_vmmmc = ax.plot(lat, 1e-15*vmmmc, color='g', label='MMC')
lp_vmse = ax.plot(lat, 1e-15*vmse, color='b', label='SE')
lp_vmte = ax.plot(lat, 1e-15*vmte, color='r', label='TE')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('CCSM4 piControl')
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
