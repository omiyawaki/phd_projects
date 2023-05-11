import os
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

var='t'
fam='echam'
mean='.djfmean'
la=80

run='rp000188'
trun='0040_0139'

ref='rp000134'
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

ds=ds.mean(['lon'],skipna=True)
ds0=ds0.mean(['lon'],skipna=True)
# ds=ds.where(ds.plev < dsp.aps).mean(['lon'],skipna=True)
# ds0=ds0.where(ds0.lev < dsp0.aps).mean(['lon'],skipna=True)

print(ds.t.shape)
print(ds0.t.shape)

t=ds.t
t0=ds0.t.mean(['time'],keepdims=True,skipna=True)

print(t.shape)
print(t0.shape)

dt=t.data-t0.data

# last 20 years
print(dt.shape)
dt20=np.squeeze(np.nanmean(dt[-21:-1,:,:], axis=0))
print(dt20.shape)

# mesh
[plev, lat] = np.meshgrid(1e-2*ds.plev, ds.lat)

fig,ax=plt.subplots(figsize=(4,3))
clf=ax.contourf(lat, plev, np.transpose(dt20), np.arange(-20,21,1), cmap='RdBu_r', vmin=-20, vmax=20, extend='both')
cb=plt.colorbar(clf)
cb.set_ticks(np.arange(-20,21,10))
cb.set_label(r'$\Delta T$ (K)')
ax.set_xticks(np.arange(-90,91,30))
ax.set_yticks(np.arange(0,1001,200))
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel(r'p (hPa)')
ax.set_title('RCP8.5 AQUAice, DJF')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig('%s/dta%s.pdf' % (pd, mean), format='pdf', dpi=300)
