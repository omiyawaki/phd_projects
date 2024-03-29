import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

ds = xr.open_dataset('./qflux_ice40m_e.nc')
q = np.nanmean(ds.aflux.data, axis=2) # zonal mean
lat = ds.lat.data
mon = np.arange(12)
q_ann = np.nanmean(q, axis=0)

[mesh_lat, mesh_mon] = np.meshgrid(lat, mon)

fig, ax = plt.subplots()
vmax = np.max(q)
vmin = -vmax
csf=ax.contourf(mesh_mon, mesh_lat, q, levels=np.arange(-700,701,100), cmap='RdBu', vmin=vmin, vmax=vmax)
# ax.contour(mesh_mon, mesh_lat, q, levels=[0], colors='k')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_xticks(np.arange(12))
ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.set_yticks(np.arange(-90,91,30))
ax.set_ylabel('Latitude (deg)')
ax.yaxis.set_minor_locator(MultipleLocator(10))
cbar = plt.colorbar(csf)
cbar.set_label('$Q_{\mathrm{flux}}$ (W m$^{-2}$)')
fig.set_size_inches(5,3)
plt.tight_layout()
plt.savefig('./qflux_e.pdf', format='pdf', dpi=300)

fig, ax = plt.subplots()
vmax = np.max(q)
vmin = -vmax
ax.axhline(0,color='k',linewidth=0.5)
ax.plot(lat, q_ann, '-k')
hl = np.where(lat > 80)[0]
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
fig.set_size_inches(4,3)
plt.tight_layout()
plt.savefig('./qflux_ann_e.pdf', format='pdf', dpi=300)

# orig

ds = xr.open_dataset('./qflux_ice40m_b.nc')
qorig = np.nanmean(ds.aflux.data, axis=2) # zonal mean
qo_ann = np.nanmean(qorig, axis=0)


# DIFF
dq = q - qorig
fig, ax = plt.subplots()
vmax = np.max(dq)
vmin = -vmax
csf=ax.contourf(mesh_mon, mesh_lat, dq, cmap='RdBu', vmin=vmin, vmax=vmax)
ax.contour(mesh_mon, mesh_lat, dq, levels=[0], colors='k')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_xticks(np.arange(12))
ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.set_yticks(np.arange(-90,91,30))
ax.set_ylabel('Latitude (deg)')
ax.yaxis.set_minor_locator(MultipleLocator(10))
cbar = plt.colorbar(csf)
cbar.set_label('$\Delta Q_{\mathrm{flux}}$ (W m$^{-2}$)')
fig.set_size_inches(5,3)
plt.tight_layout()
plt.savefig('./qflux_e_diff.pdf', format='pdf', dpi=300)

fig, ax = plt.subplots()
ax.axhline(0,color='k',linewidth=0.5)
ax.plot(lat, q_ann-qo_ann, '-k')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
fig.set_size_inches(4,3)
plt.tight_layout()
plt.savefig('./qflux_ann_e_diff.pdf', format='pdf', dpi=300)

