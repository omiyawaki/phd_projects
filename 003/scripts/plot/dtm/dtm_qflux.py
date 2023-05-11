import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

var='tmax'
fam='echam'
mean='.djfmean'
la=80

run='rp000191b'
trun='0040_0339'

ref='rp000190b'
tref='0020_0039'
yref=1987

prefix='/project2/tas1/miyawaki/projects/003'
dd='%s/data/raw/%s/%s' % (prefix, fam, run)
dd0='%s/data/raw/%s/%s' % (prefix, fam, ref)

fts='%s/%s_%s_%s%s.nc' % (dd, var, run, trun, mean)
fts0='%s/%s_%s_%s%s.nc' % (dd0, var, ref, tref, mean)

pd='%s/plot/%s/%s/%s/%s' % (prefix, fam, run, trun, var)
if not os.path.exists(pd):
    os.makedirs(pd)

ds=xr.open_dataset(fts)
ds0=xr.open_dataset(fts0)

ds=ds.mean(['lon'])
ds0=ds0.mean(['lon'])

ds=ds.where(ds.lat>=la, drop=True)
ts=ds.tmax.weighted(np.cos(np.deg2rad(ds.lat))).mean('lat')
ds0=ds0.where(ds0.lat>=la, drop=True)
ts0=ds0.tmax.weighted(np.cos(np.deg2rad(ds0.lat))).mean(['time', 'lat'])

dts=ts-ts0
time=yref+np.arange(dts.size)

fix,ax=plt.subplots(figsize=(4,3))
ax.axhline(0, linewidth=0.5, color='k')
ax.plot(time[1:-1], dts[1:-1], '-k')
ax.set_xlabel('Time (yr)')
ax.set_ylabel(r'$\Delta T_{max}$ (K)')
ax.set_title('RCP8.5 AQUAqflux, DJF, $80-90^\circ$N')
plt.tight_layout()
plt.savefig('%s/dtm%s.pdf' % (pd, mean), format='pdf', dpi=300)
