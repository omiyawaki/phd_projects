import os
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

var='t'
fam='echam'
mean='.djfmean'
la=80

run='rp000191b'
trun='0040_0239'

ref='rp000190b'
tref='0020_0039'
yref=1987

prefix='/project2/tas1/miyawaki/projects/003'
dd='%s/data/raw/%s/%s' % (prefix, fam, run)
dd0='%s/data/raw/%s/%s' % (prefix, fam, ref)

ft='%s/%s_%s_%s%s.nc' % (dd, var, run, trun, mean)
ft0='%s/%s_%s_%s%s.nc' % (dd0, var, ref, tref, mean)

fps='%s/%s_%s_%s%s.nc' % (dd, 'aps', run, trun, mean)
fps0='%s/%s_%s_%s%s.nc' % (dd0, 'aps', ref, tref, mean)

pd='%s/plot/%s/%s/%s/%s' % (prefix, fam, run, trun, var)
if not os.path.exists(pd):
    os.makedirs(pd)

ds=xr.open_dataset(ft)
ds0=xr.open_dataset(ft0)
dsp=xr.open_dataset(fps)
dsp0=xr.open_dataset(fps0)

plev=ds['plev']

ds=ds.mean(['lon'],skipna=True)
ds0=ds0.mean(['lon'],skipna=True)
# ds=ds.where(ds.plev < dsp.aps).mean(['lon'],skipna=True)
# ds0=ds0.where(ds0.lev < dsp0.aps).mean(['lon'],skipna=True)

# high latitude mean
ds=ds.where(ds.lat>=la, drop=True)
t=ds.t.weighted(np.cos(np.deg2rad(ds.lat))).mean('lat')

ds0=ds0.where(ds0.lat>=la, drop=True)
t0=ds0.t.weighted(np.cos(np.deg2rad(ds0.lat))).mean('lat')

# historical 20 year mean
t0=np.squeeze(t0.mean(['time'],keepdims=True,skipna=True))

# last 20 years of second to last century
t20n2=np.squeeze(np.nanmean(t[-121:-101,:], axis=0))

# last 20 years of last century
t20n1=np.squeeze(np.nanmean(t[-21:-1,:], axis=0))

fig,ax=plt.subplots(figsize=(4,3))
ax.plot(t0, 1e-2*plev, label='Year $20--40$')
ax.plot(t20n2, 1e-2*plev, label='Year $120--140$')
ax.plot(t20n1, 1e-2*plev, label='Year $220--240$')
# ax.set_xticks(np.arange(-90,91,30))
ax.set_ylim([200,1000])
ax.set_yticks(np.arange(200,1001,200))
ax.set_xlabel('$T$ (K)')
ax.set_ylabel(r'p (hPa)')
ax.set_title('AQUAqflux, DJF')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig('%s/taslice%s.pdf' % (pd, mean), format='pdf', dpi=300)

fig,ax=plt.subplots(figsize=(4,3))
ax.axvline(0, linewidth=0.5, color='k')
ax.plot(t20n2-t0, 1e-2*plev, label='Year $120--140$')
ax.plot(t20n1-t0, 1e-2*plev, label='Year $220--240$')
# ax.set_xticks(np.arange(-90,91,30))
ax.set_xlim([-np.max(t20n1-t0),np.max(t20n1-t0)])
ax.set_ylim([200,1000])
ax.set_yticks(np.arange(200,1001,200))
ax.set_xlabel('$\Delta T$ (K)')
ax.set_ylabel(r'p (hPa)')
ax.set_title('AQUAqflux, DJF')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig('%s/dtaslice%s.pdf' % (pd, mean), format='pdf', dpi=300)
