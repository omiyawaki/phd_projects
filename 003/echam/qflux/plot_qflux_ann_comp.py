import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# ds = xr.open_dataset('./qflux_ice40m.nc')
# q = np.nanmean(ds.aflux.data, axis=2) # zonal mean
# q_ann = np.nanmean(q, axis=0)

ds = xr.open_dataset('./qflux_ice40m_b.nc')
lat = ds.lat.data
q = np.nanmean(ds.aflux.data, axis=2) # zonal mean
q_ann_b = np.nanmean(q, axis=0)

ds = xr.open_dataset('./qflux_ice40m_c.nc')
q = np.nanmean(ds.aflux.data, axis=2) # zonal mean
q_ann_c = np.nanmean(q, axis=0)

fig, ax = plt.subplots()
vmax = np.max(q)
vmin = -vmax
ax.plot(lat, q_ann_b, label='$SW_{SFC,\Delta Q}-SW_{SFC,i}$')
ax.plot(lat, q_ann_c, label=r'$SW_{SFC,i}\frac{\beta_o}{\beta_i}$')
# ax.contour(mesh_mon, mesh_lat, q, levels=[0], colors='k')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
# ax.set_xticks(np.arange(12))
# ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.set_ylabel(r'$\overline{Q}$ (W m$^{-2}$)')
plt.legend()
fig.set_size_inches(4,3)
plt.tight_layout()
plt.savefig('./qflux_ann_comp.pdf', format='pdf', dpi=300)
