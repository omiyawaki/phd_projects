import os
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

var='t'
fam='echam'
mean='.djfmean'
meanlabel='DJF'
la=80

ice='rp000134'
tice='0020_0039'

noice='rp000135'
tnoice='0020_0039'

qc='rp000190a'
tqc='0020_0039'

q='rp000190b'
tq='0020_0039'

prefix='/project2/tas1/miyawaki/projects/003'
dd='%s/data/raw/%s/%s' % (prefix, fam, ice)
dd0='%s/data/raw/%s/%s' % (prefix, fam, noice)
ddqc='%s/data/raw/%s/%s' % (prefix, fam, qc)
ddq='%s/data/raw/%s/%s' % (prefix, fam, q)

ft='%s/%s_%s_%s%s.nc' % (dd, var, ice, tice, mean)
ft0='%s/%s_%s_%s%s.nc' % (dd0, var, noice, tnoice, mean)
ftqc='%s/%s_%s_%s%s.nc' % (ddqc, var, qc, tqc, mean)
ftq='%s/%s_%s_%s%s.nc' % (ddq, var, q, tq, mean)

fps='%s/%s_%s_%s%s.nc' % (dd, 'aps', ice, tice, mean)
fps0='%s/%s_%s_%s%s.nc' % (dd0, 'aps', noice, tnoice, mean)
fpsqc='%s/%s_%s_%s%s.nc' % (ddqc, 'aps', qc, tqc, mean)
fpsq='%s/%s_%s_%s%s.nc' % (ddq, 'aps', q, tq, mean)

pd='%s/plot/%s/%s/%s/%s' % (prefix, fam, ice, tice, var)
if not os.path.exists(pd):
    os.makedirs(pd)

ds=xr.open_dataset(ft)
ds0=xr.open_dataset(ft0)
dsqc=xr.open_dataset(ftqc)
dsq=xr.open_dataset(ftq)
dsp=xr.open_dataset(fps)
dsp0=xr.open_dataset(fps0)
dspqc=xr.open_dataset(fpsqc)
dspq=xr.open_dataset(fpsq)

plev=1e-2*ds['lev']
lat=ds['lat']

ds=ds.mean(['time', 'lon'],skipna=True)
ds0=ds0.mean(['time', 'lon'],skipna=True)
dsqc=dsqc.mean(['time', 'lon'],skipna=True)
dsq=dsq.mean(['time', 'lon'],skipna=True)

# mesh
mlev, mlat = np.meshgrid(plev, lat)

# ice
fig,ax=plt.subplots(figsize=(4,3))
clf=ax.contourf(mlat, mlev, np.transpose(ds['t'].data), np.arange(195,306,5), vmin=195, vmax=305, extend='both', cmap='Spectral_r')
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel(r'latitude (deg)')
ax.set_ylim([100,1000])
ax.set_yticks(np.arange(100,1001,100))
ax.set_ylabel(r'p (hPa)')
ax.set_title('AQUAice, %s' % meanlabel)
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator())
cb=plt.colorbar(clf)
cb.set_label('$T$ (K)')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig('%s/taz_ann%s.pdf' % (pd, mean), format='pdf', dpi=300)

# no ice
fig,ax=plt.subplots(figsize=(4,3))
clf=ax.contourf(mlat, mlev, np.transpose(ds0['t'].data), np.arange(195,306,5), vmin=195, vmax=305, extend='both', cmap='Spectral_r')
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel(r'latitude (deg)')
ax.set_ylim([100,1000])
ax.set_yticks(np.arange(100,1001,100))
ax.set_ylabel(r'p (hPa)')
ax.set_title('AQUAnoice, %s' % meanlabel)
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator())
cb=plt.colorbar(clf)
cb.set_label('$T$ (K)')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig('%s/taz_noice_ann%s.pdf' % (pd, mean), format='pdf', dpi=300)

# qc
fig,ax=plt.subplots(figsize=(4,3))
clf=ax.contourf(mlat, mlev, np.transpose(dsqc['t'].data), np.arange(195,306,5), vmin=195, vmax=305, extend='both', cmap='Spectral_r')
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel(r'latitude (deg)')
ax.set_ylim([100,1000])
ax.set_yticks(np.arange(100,1001,100))
ax.set_ylabel(r'p (hPa)')
ax.set_title(r'AQUAqflux ($Q_C$), %s' % meanlabel)
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator())
cb=plt.colorbar(clf)
cb.set_label('$T$ (K)')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig('%s/taz_qc_ann%s.pdf' % (pd, mean), format='pdf', dpi=300)

# q
fig,ax=plt.subplots(figsize=(4,3))
clf=ax.contourf(mlat, mlev, np.transpose(dsq['t'].data), np.arange(195,306,5), vmin=195, vmax=305, extend='both', cmap='Spectral_r')
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel(r'latitude (deg)')
ax.set_ylim([100,1000])
ax.set_yticks(np.arange(100,1001,100))
ax.set_ylabel(r'p (hPa)')
ax.set_title(r'AQUAqflux ($Q_C + Q_\alpha$), %s' % meanlabel)
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator())
cb=plt.colorbar(clf)
cb.set_label('$T$ (K)')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig('%s/taz_q_ann%s.pdf' % (pd, mean), format='pdf', dpi=300)

# DIFF

# ice minus no ice
fig,ax=plt.subplots(figsize=(4,3))
clf=ax.contourf(mlat, mlev, np.transpose(ds['t'].data-ds0['t'].data), np.arange(-20,20.1,1), vmin=-20, vmax=20, extend='both', cmap='RdBu_r')
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel(r'latitude (deg)')
ax.set_ylim([100,1000])
ax.set_yticks(np.arange(100,1001,100))
ax.set_ylabel(r'p (hPa)')
ax.set_title('AQUAice$-$AQUAnoice, %s' % meanlabel)
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator())
cb=plt.colorbar(clf)
cb.set_label('$\Delta T$ (K)')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig('%s/taz_noice_minus_ice_ann%s.pdf' % (pd, mean), format='pdf', dpi=300)

# qc minus no ice
fig,ax=plt.subplots(figsize=(4,3))
clf=ax.contourf(mlat, mlev, np.transpose(dsqc['t'].data-ds0['t'].data), np.arange(-20,20.1,1), vmin=-20, vmax=20, extend='both', cmap='RdBu_r')
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel(r'latitude (deg)')
ax.set_ylim([100,1000])
ax.set_yticks(np.arange(100,1001,100))
ax.set_ylabel(r'p (hPa)')
ax.set_title('AQUAqflux ($Q_C$)$-$AQUAnoice, %s' % meanlabel)
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator())
cb=plt.colorbar(clf)
cb.set_label('$\Delta T$ (K)')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig('%s/taz_noice_minus_qc_ann%s.pdf' % (pd, mean), format='pdf', dpi=300)

# q minus no ice
fig,ax=plt.subplots(figsize=(4,3))
clf=ax.contourf(mlat, mlev, np.transpose(dsq['t'].data-ds0['t'].data), np.arange(-20,20.1,1), vmin=-20, vmax=20, extend='both', cmap='RdBu_r')
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel(r'latitude (deg)')
ax.set_ylim([100,1000])
ax.set_yticks(np.arange(100,1001,100))
ax.set_ylabel(r'p (hPa)')
ax.set_title(r'AQUAqflux ($Q_C+Q_\alpha$)$-$AQUAnoice, %s' % meanlabel)
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator())
cb=plt.colorbar(clf)
cb.set_label('$\Delta T$ (K)')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig('%s/taz_noice_minus_q_ann%s.pdf' % (pd, mean), format='pdf', dpi=300)

# q minus ice
fig,ax=plt.subplots(figsize=(4,3))
clf=ax.contourf(mlat, mlev, np.transpose(dsq['t'].data-ds['t'].data), np.arange(-20,20.1,1), vmin=-20, vmax=20, extend='both', cmap='RdBu_r')
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel(r'latitude (deg)')
ax.set_ylim([100,1000])
ax.set_yticks(np.arange(100,1001,100))
ax.set_ylabel(r'p (hPa)')
ax.set_title(r'AQUAqflux ($Q_C+Q_\alpha$)$-$AQUAice, %s' % meanlabel)
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(AutoMinorLocator())
cb=plt.colorbar(clf)
cb.set_label('$\Delta T$ (K)')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig('%s/taz_ice_minus_q_ann%s.pdf' % (pd, mean), format='pdf', dpi=300)

