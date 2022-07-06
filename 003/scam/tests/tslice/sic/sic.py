import sys
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

co2 = np.arange(300,1201,100, dtype='float')
var = "ICEFRAC"

# load files
sic = np.empty_like(co2)
qsic = np.empty_like(co2)

i = 0
for c in co2:
    path_sic = '../../../data/proc/%s/co2vmr_%g.%s.djfmean.nc' % (var, c, var)
    path_qsic = '../../../data/proc/%s/qco2vmr_%g.%s.djfmean.nc' % (var, c, var)

    file_sic = Dataset(path_sic, 'r')
    file_qsic = Dataset(path_qsic, 'r')

    sic[i] = np.squeeze(file_sic.variables[var][:]) # (time)
    qsic[i] = np.squeeze(file_qsic.variables[var][:]) # (time)

    i = i + 1

# plot

# sic
fig, ax = plt.subplots()
ax.plot(co2, sic, color='tab:gray', label='prog q')
ax.plot(co2, qsic, color='tab:blue', label='fix q')
ax.set_xlabel('pCO2 (ppmv)')
ax.set_ylabel('$f_{ice}$ (unitless)')
ax.legend()
fig.set_size_inches(4,3)
plt.tight_layout()
plt.savefig('p_sic.pdf', format='pdf', dpi=300)
