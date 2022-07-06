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
v_r8['hfls'] = np.squeeze(file_r8.variables['LHFLX'][:]) # (time)
v_r8['hfss'] = np.squeeze(file_r8.variables['SHFLX'][:]) # (time)

v_fq = {}
v_fq['hfls'] = np.squeeze(file_fq.variables['LHFLX'][:]) # (time)
v_fq['hfss'] = np.squeeze(file_fq.variables['SHFLX'][:]) # (time)

# plot

# ra
fig, ax = plt.subplots()
ax.plot(v_r8['hfls'], color='tab:blue', label='LH')
ax.plot(v_fq['hfss'], color='tab:orange', label='SH')
ax.set_xlabel('time (years)')
ax.set_ylabel('energy flux (W m$^{-2}$)')
ax.legend()
fig.set_size_inches(4,3)
plt.tight_layout()
plt.savefig('p_stf.pdf', format='pdf', dpi=300)
