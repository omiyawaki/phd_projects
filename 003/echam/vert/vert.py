import os
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

vn='t'; 

sn='ann'; sni=np.arange(12);
ts='ANN Temperature'

# sn='djf'; sni=[0,1,11];
# ts='DJF Temperature'

# sn='jja'; sni=[5,6,7];
# ts='JJA Temperature'

la=80
pd='./plot/%s' % (vn)
if not os.path.exists(pd):
    os.makedirs(pd)
pn='%s/%s.%s.pdf' % (pd, vn, sn)

ri='rp000134';  ti='0020_0039'; mi='.ymonmean-20';
ra='rp000190f'; ta='0020_0039'; ma='.ymonmean-20';
# rb='rp000190b'; tb='0001_0009'; mb='.ymonmean-5'

dd='/project2/tas1/miyawaki/projects/003/data/raw/echam'
ddi='%s/%s' % (dd, ri)
dda='%s/%s' % (dd, ra)
# ddb='%s/%s' % (dd, rb)

fni='%s/%s_%s_%s%s.nc' % (ddi, vn, ri, ti, mi)
fna='%s/%s_%s_%s%s.nc' % (dda, vn, ra, ta, ma)
# fnb='%s/%s_%s_%s%s.nc' % (ddb, vn, rb, tb, mb)

dsi=xr.open_dataset(fni)
dsa=xr.open_dataset(fna)
# dsb=xr.open_dataset(fnb)

# take djf mean and zonal mean
dsi=dsi.isel(time=sni).mean(['time','lon'])
dsa=dsa.isel(time=sni).mean(['time','lon'])
# dsb=dsb.isel(time=sni).mean(['time','lon'])

# take Arctic mean
dsi=dsi.where(dsi.lat >= la, drop=True)
tai=dsi.t.weighted(np.cos(np.deg2rad(dsi.lat))).mean('lat')
dsa=dsa.where(dsa.lat >= la, drop=True)
taa=dsa.t.weighted(np.cos(np.deg2rad(dsa.lat))).mean('lat')
# dsb=dsb.where(dsb.lat >= la, drop=True)
# tab=dsb.t.weighted(np.cos(np.deg2rad(dsb.lat))).mean('lat')

fig,ax = plt.subplots(figsize=(3.5,4))
ax.plot(tai, 1e-2*dsi.lev, color='tab:purple', label='AQUAice')
ax.plot(taa, 1e-2*dsi.lev, color='tab:blue', label=r'AQUAnoice')
# ax.plot(tab, 1e-2*dsi.lev, color='tab:purple', label=r'AQUAnoice with $Q_C + Q_{SW}$')
ax.set_ylim([100,1000])
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Pressure (hPa)')
ax.set_title(ts)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
if sn == 'ann':
    plt.legend(frameon=False, prop={'size':12})
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig(pn, format='pdf', dpi=300)
