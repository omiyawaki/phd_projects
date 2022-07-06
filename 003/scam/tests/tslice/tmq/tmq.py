import sys
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

co2 = np.arange(300,1201,100, dtype='float')
var = "TMQ"

# load files
tmq = np.empty_like(co2)
qtmq = np.empty_like(co2)

i = 0
for c in co2:
    path_tmq = '../../../data/proc/%s/co2vmr_%g.%s.djfmean.nc' % (var, c, var)
    path_qtmq = '../../../data/proc/%s/qco2vmr_%g.%s.djfmean.nc' % (var, c, var)

    file_tmq = Dataset(path_tmq, 'r')
    file_qtmq = Dataset(path_qtmq, 'r')

    tmq[i] = file_tmq.variables[var][:] # (time)
    qtmq[i] = file_qtmq.variables[var][:] # (time)

    i = i + 1

# plot

# tmq
fig, ax = plt.subplots()
ax.plot(co2, tmq, color='tab:gray', label='prog q')
ax.plot(co2, qtmq, color='tab:blue', label='fix q')
ax.set_xlabel('pCO2 (ppmv)')
ax.set_ylabel(r'$\langle q \rangle$ (kg m$^{-2}$)')
ax.legend()
fig.set_size_inches(4,3)
plt.tight_layout()
plt.savefig('p_tmq.pdf', format='pdf', dpi=300)
