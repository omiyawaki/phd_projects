import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# q0
ds = xr.open_dataset('./qflux_ice40m_b.nc')
lat = ds.lat.data
q0 = np.nanmean(ds.aflux.data, axis=2) # zonal mean
q0 = np.roll(q0,1,axis=0)[0:3,:]
q0a = np.nanmean(q0, axis=0)
# q1
ds = xr.open_dataset('./qflux_ice40m_d.nc')
q1 = np.nanmean(ds.aflux.data, axis=2) # zonal mean
q1 = np.roll(q1,1,axis=0)[0:3,:]
q1a = np.nanmean(q1, axis=0)
# q2
ds = xr.open_dataset('./qflux_ice40m_e.nc')
q2 = np.nanmean(ds.aflux.data, axis=2) # zonal mean
q2 = np.roll(q2,1,axis=0)[0:3,:]
q2a = np.nanmean(q2, axis=0)
# q3
ds = xr.open_dataset('./qflux_ice40m_f.nc')
q3 = np.nanmean(ds.aflux.data, axis=2) # zonal mean
q3 = np.roll(q3,1,axis=0)[0:3,:]
q3a = np.nanmean(q3, axis=0)
# q4
ds = xr.open_dataset('./qflux_ice40m_g.nc')
q4 = np.nanmean(ds.aflux.data, axis=2) # zonal mean
q4 = np.roll(q4,1,axis=0)[0:3,:]
q4a = np.nanmean(q4, axis=0)

fig, ax = plt.subplots()
ax.axhline(0,color='k',linewidth=0.5)
ax.plot(lat, q0a, label=r'$Q^0$')
ax.plot(lat, q1a, label=r'$Q^1$')
ax.plot(lat, q2a, label=r'$Q^2$')
ax.plot(lat, q3a, label=r'$Q^3$')
ax.plot(lat, q4a, label=r'$Q^4$')
ax.set_xticks(np.arange(-90,90+30,30))
ax.set_xlim([-90,90])
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('Q flux (W m$^{-2}$)')
ax.set_title('DJF mean')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
fig.set_size_inches(4,3)
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig('./comp_qflux_djf.pdf', format='pdf', dpi=300)

fig, ax = plt.subplots()
ax.axhline(0,color='k',linewidth=0.5)
ax.plot(lat, q1a-q0a,color='tab:orange', label=r'$Q^1-Q^0$')
ax.plot(lat, q2a-q0a,color='tab:green', label=r'$Q^2-Q^0$')
ax.plot(lat, q3a-q0a,color='tab:red', label=r'$Q^3-Q^0$')
ax.plot(lat, q4a-q0a,color='tab:red', label=r'$Q^4-Q^0$')
ax.set_xticks(np.arange(-90,90+30,30))
ax.set_xlim([-90,90])
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('$\Delta$ Q flux (W m$^{-2}$)')
ax.set_title('DJF mean')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
fig.set_size_inches(4,3)
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig('./comp_qflux_djf_diff.pdf', format='pdf', dpi=300)

# q0
ds = xr.open_dataset('./qsw_ice40m_b.nc')
lat = ds.lat.data
q0 = np.nanmean(ds.aflux.data, axis=2) # zonal mean
q0 = np.roll(q0,1,axis=0)[0:3,:]
q0a = np.nanmean(q0, axis=0)
# q1
ds = xr.open_dataset('./qsw_ice40m_d.nc')
q1 = np.nanmean(ds.aflux.data, axis=2) # zonal mean
q1 = np.roll(q1,1,axis=0)[0:3,:]
q1a = np.nanmean(q1, axis=0)
# q2
ds = xr.open_dataset('./qsw_ice40m_e.nc')
q2 = np.nanmean(ds.aflux.data, axis=2) # zonal mean
q2 = np.roll(q2,1,axis=0)[0:3,:]
q2a = np.nanmean(q2, axis=0)
# q3
ds = xr.open_dataset('./qsw_ice40m_f.nc')
q3 = np.nanmean(ds.aflux.data, axis=2) # zonal mean
q3 = np.roll(q3,1,axis=0)[0:3,:]
q3a = np.nanmean(q3, axis=0)
# q4
ds = xr.open_dataset('./qsw_ice40m_g.nc')
q4 = np.nanmean(ds.aflux.data, axis=2) # zonal mean
q4 = np.roll(q4,1,axis=0)[0:3,:]
q4a = np.nanmean(q4, axis=0)

fig, ax = plt.subplots()
ax.axhline(0,color='k',linewidth=0.5)
ax.plot(lat, q0a, label=r'$Q^0_{SW}$')
ax.plot(lat, q1a, label=r'$Q^1_{SW}$')
ax.plot(lat, q2a, label=r'$Q^2_{SW}$')
ax.plot(lat, q3a, label=r'$Q^3_{SW}$')
ax.plot(lat, q4a, label=r'$Q^4_{SW}$')
ax.set_xticks(np.arange(-90,90+30,30))
ax.set_xlim([-90,90])
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('Shortwave component of Q flux (W m$^{-2}$)')
ax.set_title('DJF mean')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
fig.set_size_inches(4,3)
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig('./comp_qsw_djf.pdf', format='pdf', dpi=300)

fig, ax = plt.subplots()
ax.axhline(0,color='k',linewidth=0.5)
ax.plot(lat, q1a-q0a,color='tab:orange', label=r'$Q^1_{SW}-Q^0_{SW}$')
ax.plot(lat, q2a-q0a,color='tab:green',  label=r'$Q^2_{SW}-Q^0_{SW}$')
ax.plot(lat, q3a-q0a,color='tab:red',    label=r'$Q^3_{SW}-Q^0_{SW}$')
ax.plot(lat, q4a-q0a,color='tab:red',    label=r'$Q^4_{SW}-Q^0_{SW}$')
ax.set_xticks(np.arange(-90,90+30,30))
ax.set_xlim([-90,90])
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('$\Delta$ shortwave component of Q flux (W m$^{-2}$)')
ax.set_title('DJF mean')
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
fig.set_size_inches(4,3)
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig('./comp_qsw_djf_diff.pdf', format='pdf', dpi=300)
