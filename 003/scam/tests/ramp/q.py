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
v_r8['TMQ'] = np.squeeze(file_r8.variables['TMQ'][:]) # (time)

v_fq = {}
v_fq['TMQ'] = np.squeeze(file_fq.variables['TMQ'][:]) # (time)

# plot

# tmq
fig, ax = plt.subplots()
ax.plot(v_r8['TMQ'], color='tab:gray', label='all')
ax.plot(v_fq['TMQ'], color='tab:blue', label='fix q')
ax.set_xlabel('time (years)')
ax.set_ylabel(r'$\langle q \rangle$ (kg m$^{-2}$)')
ax.legend()
fig.set_size_inches(4,3)
plt.tight_layout()
plt.savefig('p_tmq.pdf', format='pdf', dpi=300)
