import sys
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

co2 = np.arange(300,1201,100,dtype='float')

# load files
ra = np.empty_like(co2)
qra = np.empty_like(co2)

i = 0
for c in co2:
    path_ra = '../../../data/proc/ra/co2vmr_%g.ra.djfmean.nc' % (c)
    path_qra = '../../../data/proc/ra/qco2vmr_%g.ra.djfmean.nc' % (c)

    file_ra = Dataset(path_ra, 'r')
    file_qra = Dataset(path_qra, 'r')

    ra[i] = file_ra.variables['ra'][:] # (time)
    qra[i] = file_qra.variables['ra'][:] # (time)

    i = i + 1

# plot

# ra
fig, ax = plt.subplots()
ax.plot(co2, ra, color='tab:gray', label='prog q')
ax.plot(co2, qra, color='tab:blue', label='fix q')
ax.set_xlabel('pCO2 (ppmv)')
ax.set_ylabel('$R_a$ (W m$^{-2}$)')
ax.legend()
fig.set_size_inches(4,3)
plt.tight_layout()
plt.savefig('p_ra.pdf', format='pdf', dpi=300)
