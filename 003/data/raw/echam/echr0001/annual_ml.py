import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from netCDF4 import Dataset

##############################
# ANNUAL MEAN, ALL TRANSPORTS
##############################

vmmmc_file = Dataset('./vmmmc_6h_ml_echr0001_1030.nc', 'r')
vmmmc = np.squeeze(np.mean(vmmmc_file.variables['vmmmc'][:],axis=0))
lat = vmmmc_file.variables['lat'][:]

vmse_file = Dataset('./vmse_6h_ml_echr0001_1030.nc', 'r')
vmse = np.squeeze(np.mean(vmse_file.variables['vmse'][:],axis=0))

vmte_file = Dataset('./vmte_6h_ml_echr0001_1030.nc', 'r')
vmte = np.squeeze(np.mean(vmte_file.variables['vmte'][:],axis=0))

aht_file = Dataset('./aht_6h_ml_echr0001_1030.monmean.nc', 'r')
aht = np.squeeze(np.mean(aht_file.variables['aht'][:],axis=0))
lat_aht = aht_file.variables['lat'][:]

vmtot = vmmmc + vmse + vmte

plotname = './annual_ml'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
lp_vmmmc = ax.plot(lat, 1e-15*vmmmc, color='g', label='MMC')
lp_vmse = ax.plot(lat, 1e-15*vmse, color='b', label='SE')
lp_vmte = ax.plot(lat, 1e-15*vmte, color='r', label='TE')
lp_vmtot = ax.plot(lat, 1e-15*vmtot, color='k', label='Total AHT')
lp_aht = ax.plot(lat_aht, 1e-15*aht, ':', color='k', label='Total AHT (inferred)')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_title('ECHAM6, ANN')
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
