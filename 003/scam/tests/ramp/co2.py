import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# load files
path_r8 = '../data/tmp_rcp85.17.djfmean.nc'
path_fq = '../data/tmp_fixq.17.djfmean.nc'

file_r8 = Dataset(path_r8, 'r')
file_fq = Dataset(path_fq, 'r')

# load variables
v_r8 = {}
v_r8['co2'] = file_r8.variables['co2vmr'][:] # (time)

# plot

# co2
print(v_r8['co2'])
fig, ax = plt.subplots()
ax.plot(1e6*v_r8['co2'], '-k')
ax.set_xlabel('time (years)')
ax.set_ylabel('$pCO_2$ (ppmv)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.savefig('p_co2.pdf', format='pdf', dpi=300)

## REPEAT FOR FULL TIMESERIES
# load files
path_r8 = '../data/tmp_rcp85.17.nc'
path_fq = '../data/tmp_fixq.17.nc'

file_r8 = Dataset(path_r8, 'r')
file_fq = Dataset(path_fq, 'r')

# load variables
v_r8 = {}
v_r8['co2'] = file_r8.variables['co2vmr'][:] # (time)

# plot

# co2
print(v_r8['co2'])
fig, ax = plt.subplots()
ax.plot(1e6*v_r8['co2'], '-k')
ax.set_xlabel('time (months)')
ax.set_ylabel('$pCO_2$ (ppmv)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.savefig('p_co2_mon.pdf', format='pdf', dpi=300)
