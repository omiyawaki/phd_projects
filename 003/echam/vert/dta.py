import os
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

vn='t'; 

# sn='ann'; sni=np.arange(12);
# ts='ANN'
# tsi='AQUAice ANN'
# tsn='AQUAqflux ANN'

sn='djf'; sni=[0,1,11];
ts='DJF'
tsi='AQUAice DJF'
tsn='AQUAnoice DJF'

# sn='jja'; sni=[5,6,7];
# ts='AQUAice JJA'

la=70
pd='./plot/%s' % (vn)
if not os.path.exists(pd):
    os.makedirs(pd)

ri='rp000134';  ti='0020_0039'; mi='.ymonmean-20';
ra='rp000188'; ta='0040_0252'; ma='.ymonmean-20';
rn='rp000190f';  tn='0020_0039'; mn='.ymonmean-20';
rf='rp000191f'; tf='0040_0252'; mf='.ymonmean-20';

pi='%s/%s.%s.%s.pdf' % (pd, ri,vn, sn)
pn='%s/%s.%s.%s.pdf' % (pd, rn,vn, sn)
dp='%s/d%s.%s.pdf' % (pd, vn, sn)
dpf='%s/df%s.%s.pdf' % (pd, vn, sn)

dd='/project2/tas1/miyawaki/projects/003/data/raw/echam'
ddi='%s/%s' % (dd, ri)
dda='%s/%s' % (dd, ra)
ddn='%s/%s' % (dd, rn)
ddf='%s/%s' % (dd, rf)

fni='%s/%s_%s_%s%s.nc' % (ddi, vn, ri, ti, mi)
fna='%s/%s_%s_%s%s.nc' % (dda, vn, ra, ta, ma)
fnn='%s/%s_%s_%s%s.nc' % (ddn, vn, rn, tn, mn)
fnf='%s/%s_%s_%s%s.nc' % (ddf, vn, rf, tf, mf)

dsi=xr.open_dataset(fni)
dsa=xr.open_dataset(fna)
dsn=xr.open_dataset(fnn)
dsf=xr.open_dataset(fnf)

# take djf mean and zonal mean
dsi=dsi.isel(time=sni).mean(['time','lon'])
dsa=dsa.isel(time=sni).mean(['time','lon'])
dsn=dsn.isel(time=sni).mean(['time','lon'])
dsf=dsf.isel(time=sni).mean(['time','lon'])

# take Arctic mean
dsi=dsi.where(dsi.lat >= la, drop=True)
tai=dsi.t.weighted(np.cos(np.deg2rad(dsi.lat))).mean('lat')
dsa=dsa.where(dsa.lat >= la, drop=True)
taa=dsa.t.weighted(np.cos(np.deg2rad(dsa.lat))).mean('lat')
dsn=dsn.where(dsn.lat >= la, drop=True)
tan=dsn.t.weighted(np.cos(np.deg2rad(dsn.lat))).mean('lat')
dsf=dsf.where(dsf.lat >= la, drop=True)
taf=dsf.t.weighted(np.cos(np.deg2rad(dsf.lat))).mean('lat')

fig,ax = plt.subplots(figsize=(3.5,4))
ax.plot(tai, 1e-2*dsi.lev, color='tab:purple', label='Control')
ax.plot(taa, 1e-2*dsi.lev, color='tab:red', label=r'RCP8.5')
# ax.plot(tab, 1e-2*dsi.lev, color='tab:purple', label=r'AQUAnoice with $Q_C + Q_{SW}$')
ax.set_ylim([100,1000])
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Pressure (hPa)')
ax.set_title(tsi)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
if sn == 'ann':
    plt.legend(frameon=False, fontsize='small')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig(pi, format='pdf', dpi=300)

fig,ax = plt.subplots(figsize=(3.5,4))
ax.plot(tan, 1e-2*dsi.lev, color='tab:purple', label='Control')
ax.plot(taf, 1e-2*dsi.lev, color='tab:red', label=r'RCP8.5')
# ax.plot(tab, 1e-2*dsi.lev, color='tab:purple', label=r'AQUAnoice with $Q_C + Q_{SW}$')
ax.set_ylim([100,1000])
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Pressure (hPa)')
ax.set_title(tsn)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
if sn == 'ann':
    plt.legend(frameon=False, fontsize='small')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig(pn, format='pdf', dpi=300)

fig,ax = plt.subplots(figsize=(3.5,4))
ax.axvline(0,linewidth=0.5,color='k')
ax.plot(taa.data-tai.data, 1e-2*dsi.lev, color='tab:purple',label='AQUAice')
ax.plot(taf.data-tan.data, 1e-2*dsi.lev, color='tab:blue',label='AQUAqflux')
# ax.plot(tab, 1e-2*dsi.lev, color='tab:purple', label=r'AQUAnoice with $Q_C + Q_{SW}$')
ax.set_ylim([100,1000])
ax.set_xlabel('$\Delta$ Temperature (K)')
ax.set_ylabel('Pressure (hPa)')
ax.set_title(ts)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
# if sn == 'ann':
plt.legend(frameon=False, fontsize='small')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig(dp, format='pdf', dpi=300)

fig,ax = plt.subplots(figsize=(3.5,4))
ax.axvline(0,linewidth=0.5,color='k')
ax.axvline(1e2,color='tab:purple')
ax.plot(1e2*(taa.data-tai.data-taf.data+tan.data)/(taa.data-tai.data), 1e-2*dsi.lev,'--', color='tab:purple',label='AQUAice$-$AQUAqflux')
# ax.plot(taf.data-tan.data, 1e-2*dsi.lev, color='tab:blue',label='AQUAqflux')
# ax.plot(tab, 1e-2*dsi.lev, color='tab:purple', label=r'AQUAnoice with $Q_C + Q_{SW}$')
ax.set_xlim([-10,110])
ax.set_ylim([200,1000])
ax.set_xlabel('$\Delta T/\Delta T_{AQUAice}$ (%)')
ax.set_ylabel('Pressure (hPa)')
ax.set_title(ts)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
# if sn == 'ann':
plt.legend(frameon=False, fontsize='small')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig(dpf, format='pdf', dpi=300)
