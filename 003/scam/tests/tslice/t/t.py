import sys
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

co2 = np.arange(300,1201,100, dtype='float')
var = "T"

# load files

# lev
path_t = '../../../data/proc/%s/co2vmr_%g.%s.djfmean.nc' % (var, co2[0], var)
file_t = Dataset(path_t, 'r')
lev = file_t.variables['lev'][:] # (lev)

t = np.empty([len(co2), len(lev)])
qt = np.empty([len(co2), len(lev)])

i = 0
for c in co2:
    path_t = '../../../data/proc/%s/co2vmr_%g.%s.djfmean.nc' % (var, c, var)
    path_qt = '../../../data/proc/%s/qco2vmr_%g.%s.djfmean.nc' % (var, c, var)

    file_t = Dataset(path_t, 'r')
    file_qt = Dataset(path_qt, 'r')

    t[i,:] = np.squeeze(file_t.variables[var][:]) # (time x lev)
    qt[i,:] = np.squeeze(file_qt.variables[var][:]) # (time x lev)

    i = i + 1


# plot
[mesh_lev, mesh_co2] = np.meshgrid(lev, co2) # create mesh
    ## csf = ax.contourf(mesh_lat, mesh_lev, circ['wap'], np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')

# t
vmin=200
vmax=360
fig, ax = plt.subplots()
csf=ax.contourf(mesh_co2, mesh_lev, t, np.arange(vmin,vmax,10))
ax.set_xlabel('pCO2 (ppmv)')
ax.set_ylabel('p (hPa))')
ax.set_ylim(ax.get_ylim()[::-1])
cbar = plt.colorbar(csf, ax=ax)
cbar.set_label('T (K)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.savefig('p_t_prog_q.pdf', format='pdf', dpi=300)

# t
vmin=-20
vmax=20
fig, ax = plt.subplots()
csf=ax.contourf(mesh_co2, mesh_lev, t-t[0,:], np.arange(vmin,vmax,1), vmin=vmin, vmax=vmax, cmap='RdBu_r')
ax.set_xlabel('pCO2 (ppmv)')
ax.set_ylabel('p (hPa))')
ax.set_ylim(ax.get_ylim()[::-1])
cbar = plt.colorbar(csf, ax=ax)
cbar.set_label('T (K)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.savefig('p_dt_prog_q.pdf', format='pdf', dpi=300)
