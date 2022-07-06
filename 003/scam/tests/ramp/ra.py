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
v_r8['FSNT'] = file_r8.variables['FSNT'][:] # (time)
v_r8['FLNT'] = file_r8.variables['FLNT'][:] # (time)
v_r8['FSNS'] = file_r8.variables['FSNS'][:] # (time)
v_r8['FLNS'] = file_r8.variables['FLNS'][:] # (time)

v_fq = {}
v_fq['FSNT'] = file_fq.variables['FSNT'][:] # (time)
v_fq['FLNT'] = file_fq.variables['FLNT'][:] # (time)
v_fq['FSNS'] = file_fq.variables['FSNS'][:] # (time)
v_fq['FLNS'] = file_fq.variables['FLNS'][:] # (time)

# post processing
p_r8 = {}
p_r8['ra'] = np.squeeze(v_r8['FSNT'] - v_r8['FLNT'] - v_r8['FSNS'] + v_r8['FLNS'])

p_fq = {}
p_fq['ra'] = np.squeeze(v_fq['FSNT'] - v_fq['FLNT'] - v_fq['FSNS'] + v_fq['FLNS'])

# plot

# ra
fig, ax = plt.subplots()
ax.plot(p_r8['ra'], color='tab:gray', label='all')
ax.plot(p_fq['ra'], color='tab:blue', label='fix q')
ax.set_xlabel('time (years)')
ax.set_ylabel('$R_a$ (W m$^{-2}$)')
ax.legend()
fig.set_size_inches(4,3)
plt.tight_layout()
plt.savefig('p_ra.pdf', format='pdf', dpi=300)
