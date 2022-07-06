import sys
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

co2 = np.arange(300,1201,100, dtype='float')
var = "TS"

# load files
ts = np.empty_like(co2)
qts = np.empty_like(co2)

i = 0
for c in co2:
    path_ts = '../../../data/proc/%s/co2vmr_%g.%s.djfmean.nc' % (var, c, var)
    path_qts = '../../../data/proc/%s/qco2vmr_%g.%s.djfmean.nc' % (var, c, var)

    file_ts = Dataset(path_ts, 'r')
    file_qts = Dataset(path_qts, 'r')

    ts[i] = file_ts.variables[var][:] # (time)
    qts[i] = file_qts.variables[var][:] # (time)

    i = i + 1

# plot

# ts
fig, ax = plt.subplots()
ax.plot(co2, ts, color='tab:gray', label='prog q')
ax.plot(co2, qts, color='tab:blue', label='fix q')
ax.set_xlabel('pCO2 (ppmv)')
ax.set_ylabel('$T_s$ (K)')
ax.legend()
fig.set_size_inches(4,3)
plt.tight_layout()
plt.savefig('p_ts.pdf', format='pdf', dpi=300)
