import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

latarc = 80 # equatorward latitude of Arctic domain
Lf = 3.337e5 # J kg**-1 latent heat of fusion
rhoi = 917 # kg m**-3 density of ice

prefix = '/project2/tas1/miyawaki/projects/003/data/raw/echam'
varnames = ['albedo', 'siced', 'temp2', 'tsurf', 'fsfc', 'trads', 'srads', 'sradsu', 'ahfl', 'ahfs']

ice = {}
run = 'rp000134'
for varname in varnames:
    ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
    dsvar = getattr(ds, varname)
    dsvar = dsvar.mean(dim='lon') # take zonal mean
    # take area mean over high latitudes
    dsvar = dsvar.where(dsvar.lat >= latarc, drop=True)
    w = np.cos(np.deg2rad(dsvar.lat))
    dsvarw = dsvar.weighted(w)
    dsvar = dsvarw.mean(dim='lat')
    ice[varname] = dsvar.data

noice = {}
run = 'rp000135'
for varname in varnames:
    ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
    dsvar = getattr(ds, varname)
    dsvar = dsvar.mean(dim='lon') # take zonal mean
    # take area mean over high latitudes
    dsvar = dsvar.where(dsvar.lat >= latarc, drop=True)
    w = np.cos(np.deg2rad(dsvar.lat))
    dsvarw = dsvar.weighted(w)
    dsvar = dsvarw.mean(dim='lat')
    noice[varname] = dsvar.data

qflux = {}
run = 'rp000190'
for varname in varnames:
    ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
    dsvar = getattr(ds, varname)
    dsvar = dsvar.mean(dim='lon') # take zonal mean
    # take area mean over high latitudes
    dsvar = dsvar.where(dsvar.lat >= latarc, drop=True)
    w = np.cos(np.deg2rad(dsvar.lat))
    dsvarw = dsvar.weighted(w)
    dsvar = dsvarw.mean(dim='lat')
    qflux[varname] = dsvar.data

##################################################
# PLOTS
##################################################
mons = np.arange(12)

# plot surface temperature
plotname = 'tsurf'
fig, ax = plt.subplots()
ax.plot(mons, ice['tsurf'], '-', color='tab:blue', label='40 m ice')
ax.plot(mons, noice['tsurf'], '-k', label='40 m no ice')
ax.plot(mons, qflux['tsurf'], '-', color='tab:red', label='40 m qflux')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$T_{s}$ (K)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend()
plt.savefig('./plot/%s.pdf' % (plotname) )
plt.close()

# plot surface temperature
plotname = 'temp2'
fig, ax = plt.subplots()
ax.plot(mons, ice['temp2'], '-', color='tab:blue', label='40 m ice')
ax.plot(mons, noice['temp2'], '-k', label='40 m no ice')
ax.plot(mons, qflux['temp2'], '-', color='tab:red', label='40 m qflux')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$T_{2\,m}$ (K)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend()
plt.savefig('./plot/%s.pdf' % (plotname) )
plt.close()

# plot temp tendency
dt = 30*86400

ice['ttend'] = (ice['tsurf'][2:] - ice['tsurf'][:-2])/(2*dt)
ice['ttend'] = np.insert(ice['ttend'], 0, (ice['tsurf'][1] - ice['tsurf'][-1])/(2*dt) )
ice['ttend'] = np.append(ice['ttend'], (ice['tsurf'][0] - ice['tsurf'][-2])/(2*dt))

noice['ttend'] = (noice['tsurf'][2:] - noice['tsurf'][:-2])/(2*dt)
noice['ttend'] = np.insert(noice['ttend'], 0, (noice['tsurf'][1] - noice['tsurf'][-1])/(2*dt) )
noice['ttend'] = np.append(noice['ttend'], (noice['tsurf'][0] - noice['tsurf'][-2])/(2*dt))

qflux['ttend'] = (qflux['tsurf'][2:] - qflux['tsurf'][:-2])/(2*dt)
qflux['ttend'] = np.insert(qflux['ttend'], 0, (qflux['tsurf'][1] - qflux['tsurf'][-1])/(2*dt) )
qflux['ttend'] = np.append(qflux['ttend'], (qflux['tsurf'][0] - qflux['tsurf'][-2])/(2*dt))

plotname = 'ttend'
fig, ax = plt.subplots()
ax.plot(mons, ice['ttend'], '-', color='tab:blue', label='40 m ice')
ax.plot(mons, noice['ttend'], '-k', label='40 m no ice')
ax.plot(mons, qflux['ttend'], '-', color='tab:red', label='40 m qflux')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$\partial_t T_{s}$ (K s$^{-1}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend()
plt.savefig('./plot/%s.pdf' % (plotname) )
plt.close()

# plot ice thickness tendency
dt = 30*86400

ice['stend'] = (ice['siced'][2:] - ice['siced'][:-2])/(2*dt)
ice['stend'] = np.insert(ice['stend'], 0, (ice['siced'][1] - ice['siced'][-1])/(2*dt) )
ice['stend'] = np.append(ice['stend'], (ice['siced'][0] - ice['siced'][-2])/(2*dt))

noice['stend'] = (noice['siced'][2:] - noice['siced'][:-2])/(2*dt)
noice['stend'] = np.insert(noice['stend'], 0, (noice['siced'][1] - noice['siced'][-1])/(2*dt) )
noice['stend'] = np.append(noice['stend'], (noice['siced'][0] - noice['siced'][-2])/(2*dt))

qflux['stend'] = (qflux['siced'][2:] - qflux['siced'][:-2])/(2*dt)
qflux['stend'] = np.insert(qflux['stend'], 0, (qflux['siced'][1] - qflux['siced'][-1])/(2*dt) )
qflux['stend'] = np.append(qflux['stend'], (qflux['siced'][0] - qflux['siced'][-2])/(2*dt))

plotname = 'stend'
fig, ax = plt.subplots()
ax.plot(mons, ice['stend'], '-', color='tab:blue', label='40 m ice')
ax.plot(mons, noice['stend'], '-k', label='40 m no ice')
ax.plot(mons, qflux['stend'], '-', color='tab:red', label='40 m qflux')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$\partial_t h_{i}$ (m s$^{-1}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend()
plt.savefig('./plot/%s.pdf' % (plotname) )
plt.close()

plotname = 'rlstend'
fig, ax = plt.subplots()
ax.plot(mons, rhoi*Lf*ice['stend'], '-', color='tab:blue', label='40 m ice')
ax.plot(mons, rhoi*Lf*noice['stend'], '-k', label='40 m no ice')
ax.plot(mons, rhoi*Lf*qflux['stend'], '-', color='tab:red', label='40 m qflux')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$\partial_t h_{i}$ (m s$^{-1}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend()
plt.savefig('./plot/%s.pdf' % (plotname) )
plt.close()

# plot inferred mixed layer depth
cw = 4182; rhow = 997;
ice['d'] = (ice['fsfc'] + Lf*rhoi*ice['stend']) / (cw*rhow*ice['ttend'])
ice['d'][ice['stend']<0] = np.nan # depth is infinite where ice is melting
noice['d'] = (noice['fsfc'] + Lf*rhoi*noice['stend']) / (cw*rhow*noice['ttend'])
qflux['d'] = (qflux['fsfc'] + Lf*rhoi*qflux['stend']) / (cw*rhow*qflux['ttend'])

plotname = 'd'
fig, ax = plt.subplots()
ax.plot(mons, ice['d'], '-', color='tab:blue', label='40 m ice')
ax.plot(mons, noice['d'], '-k', label='40 m no ice')
# ax.plot(mons, qflux['d'], '-', color='tab:red', label='40 m qflux')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$d_{\mathrm{inferred}}$ (m)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend()
plt.savefig('./plot/%s.pdf' % (plotname) )
plt.close()

# plot inferred net FSFCice
dinf=0.25
plotname = 'fsfci_infer'
fig, ax = plt.subplots()
ax.plot(mons, ice['fsfc'], '-', color='tab:blue', label='actual')
ax.plot(mons, cw*rhow*dinf*ice['ttend'] - Lf*rhoi*ice['stend'], '--', color='tab:blue', label='inferred')
# ax.plot(mons, qflux['d'], '-', color='tab:red', label='40 m qflux')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$F_{SFC,i}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend()
plt.savefig('./plot/%s.pdf' % (plotname) )
plt.close()

# plot inferred net FSFCwater
plotname = 'fsfcw_infer'
fig, ax = plt.subplots()
ax.plot(mons, noice['fsfc'], '-', color='k', label='actual')
ax.plot(mons, cw*rhow*40*noice['ttend'], '--', color='k', label='inferred')
# ax.plot(mons, qflux['d'], '-', color='tab:red', label='40 m qflux')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$F_{SFC,w}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend()
plt.savefig('./plot/%s.pdf' % (plotname) )
plt.close()

# plot surface sw radiation
plotname = 'srads'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['srads'], '-', color='tab:blue', label='40 m ice')
ax.plot(mons, noice['srads'], '-k', label='40 m no ice')
ax.plot(mons, qflux['srads'], '-', color='tab:red', label='40 m qflux')
print('%s wm**-2' % np.mean(qflux['srads']-ice['srads']))
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend()
plt.savefig('./plot/%s.pdf' % (plotname) )
plt.close()

# plot net surface fluxes
plotname = 'fsfc'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['fsfc'], '-', color='tab:blue', label='40 m ice')
ax.plot(mons, noice['fsfc'], '-k', label='40 m no ice')
ax.plot(mons, qflux['fsfc'], '-', color='tab:red', label='40 m qflux')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$F_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend()
plt.savefig('./plot/%s.pdf' % (plotname) )
plt.close()

# plot surface latent heat flux 
plotname = 'ahfl'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['ahfl'], '-', color='tab:blue', label='40 m ice')
ax.plot(mons, noice['ahfl'], '-k', label='40 m no ice')
ax.plot(mons, qflux['ahfl'], '-', color='tab:red', label='40 m qflux')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$LH$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend()
plt.savefig('./plot/%s.pdf' % (plotname) )
plt.close()

# plot surface sensible heat flux
plotname = 'ahfs'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['ahfs'], '-', color='tab:blue', label='40 m ice')
ax.plot(mons, noice['ahfs'], '-k', label='40 m no ice')
ax.plot(mons, qflux['ahfs'], '-', color='tab:red', label='40 m qflux')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SH$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend()
plt.savefig('./plot/%s.pdf' % (plotname) )
plt.close()

# plot surface LW radiation
plotname = 'trads'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['trads'], '-', color='tab:blue', label='40 m ice')
ax.plot(mons, noice['trads'], '-k', label='40 m no ice')
ax.plot(mons, qflux['trads'], '-', color='tab:red', label='40 m qflux')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$LW_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend()
plt.savefig('./plot/%s.pdf' % (plotname) )
plt.close()

# plot surface sw radiation
plotname = 'srads'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['srads'], '-', color='tab:blue', label='40 m ice')
ax.plot(mons, noice['srads'], '-k', label='40 m no ice')
ax.plot(mons, qflux['srads'], '-', color='tab:red', label='40 m qflux')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend()
plt.savefig('./plot/%s.pdf' % (plotname) )
plt.close()

# plot downward surface sw
ice['sradsd'] = ice['srads'] - ice['sradsu']
noice['sradsd'] = noice['srads'] - noice['sradsu']
qflux['sradsd'] = qflux['srads'] - qflux['sradsu']
plotname = 'sradsd'
fig, ax = plt.subplots()
ax.plot(mons, ice['sradsd'], '-', color='tab:blue', label='40 m ice')
ax.plot(mons, noice['sradsd'], '-k', label='40 m no ice')
ax.plot(mons, qflux['sradsd'], '-', color='tab:red', label='40 m qflux')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW^{\downarrow}_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend()
plt.savefig('./plot/%s.pdf' % (plotname) )
plt.close()

# plot approx net shortwave
ice['asrads'] =(1-ice['albedo'])*ice['sradsd']
noice['asrads'] =(1-noice['albedo'])*noice['sradsd']
qflux['asrads'] =(1-qflux['albedo'])*qflux['sradsd']
plotname = 'srads_approx'
fig, ax = plt.subplots()
ax.plot(mons, ice['asrads'], '--', color='tab:blue', label='40 m ice')
ax.plot(mons, noice['asrads'], '--k', label='40 m no ice')
ax.plot(mons, qflux['asrads'], '--', color='tab:red', label='40 m qflux')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$(1-\alpha_s)SW^{\downarrow}_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend()
plt.savefig('./plot/%s.pdf' % (plotname) )
plt.close()

# plot approx net shortwave
plotname = 'srads_comp'
fig, ax = plt.subplots()
ax.plot(mons, ice['asrads'], '--', color='tab:blue', label='40 m ice')
ax.plot(mons, noice['asrads'], '--k', label='40 m no ice')
ax.plot(mons, qflux['asrads'], '--', color='tab:red', label='40 m qflux')
ax.plot(mons, ice['srads'], '-', color='tab:blue', label='40 m ice')
ax.plot(mons, noice['srads'], '-k', label='40 m no ice')
ax.plot(mons, qflux['srads'], '-', color='tab:red', label='40 m qflux')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
# plt.legend()
plt.savefig('./plot/%s.pdf' % (plotname) )
plt.close()

print('\nice d_eff = %g m' % np.nanmean(ice['d']))
print('no ice d_eff = %g m' % np.nanmean(noice['d']))
print('amplification factor = %g' % (np.nanmean(noice['d'])/np.nanmean(ice['d'])) )
