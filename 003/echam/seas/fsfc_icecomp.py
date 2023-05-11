import os
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

latarc = 80 # equatorward latitude of Arctic domain
Lf = 3.337e5 # J kg**-1 latent heat of fusion
rhoi = 917 # kg m**-3 density of ice

prefix = '/project2/tas1/miyawaki/projects/003/data/raw/echam'
# varnames = ['smelt', 'snmelt', 'aprs', 'sntend', 'sni', 'xlvi', 'xivi', 'aclcov', 'ttend', 'stend', 'titend', 'twtend', 'sn', 'tsi', 'tsw', 'ahfcon', 'ahfres', 'albedo', 'siced', 'temp2', 'tsurf', 'fsfc', 'trads', 'srads', 'sradsu', 'ahfl', 'ahfs']
varnames = ['sni', 'xlvi', 'xivi', 'aclcov', 'ttend', 'stend', 'titend', 'twtend', 'sn', 'tsi', 'tsw', 'ahfcon', 'ahfres', 'albedo', 'siced', 'temp2', 'tsurf', 'fsfc', 'trads', 'srads', 'sradsu', 'ahfl', 'ahfs']

qfsim='rp000190f'
qlabel=r'AQUAnoice with $Q^3$'

# qfsim='rp000190a'
# qlabel=r'AQUAnoice, $Q=Q_C$'

plotdir='./plot/%s' % (qfsim) 
if not os.path.exists(plotdir):
    os.mkdir(plotdir)

ice = {}
run = 'rp000134'
ivarnames = varnames.copy()
ivarnames.insert(0, 'cond')
ivarnames.insert(0, 'tmelt')
ivarnames.insert(0, 'seaice')
for varname in ivarnames:
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
run = qfsim
for varname in varnames:
    # ds = xr.open_dataset('%s/%s/%s_%s_0001_0009.ymonmean-5.nc' % (prefix, run, varname, run))
    ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
    dsvar = getattr(ds, varname)
    dsvar = dsvar.mean(dim='lon') # take zonal mean
    # take area mean over high latitudes
    dsvar = dsvar.where(dsvar.lat >= latarc, drop=True)
    w = np.cos(np.deg2rad(dsvar.lat))
    dsvarw = dsvar.weighted(w)
    dsvar = dsvarw.mean(dim='lat')
    qflux[varname] = dsvar.data

varname='qflux'
# ds = xr.open_dataset('%s/%s/%s_%s_0001_0009.ymonmean-5.nc' % (prefix, run, varname, run))
ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
dsvar = getattr(ds, 'aflux')
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

# plot column integrated cloud liquid
plotname = 'xlvi'
fig, ax = plt.subplots()
ax.plot(mons, 1e3*ice['xlvi'], '-', color='tab:blue', label='AQUAice')
ax.plot(mons, 1e3*noice['xlvi'], '-k', label='AQUAnoice with $Q^0$')
ax.plot(mons, 1e3*qflux['xlvi'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'Column integrated cloud liquid (g m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot column integrated cloud ice
plotname = 'xivi'
fig, ax = plt.subplots()
ax.plot(mons, 1e3*ice['xivi'], '-', color='tab:blue', label='AQUAice')
ax.plot(mons, 1e3*noice['xivi'], '-k', label='AQUAnoice with $Q^0$')
ax.plot(mons, 1e3*qflux['xivi'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'Column integrated cloud ice (g m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot cloud cover
plotname = 'aclcov'
fig, ax = plt.subplots()
ax.plot(mons, 1e2*ice['aclcov'], '-', color='tab:blue', label='AQUAice')
ax.plot(mons, 1e2*noice['aclcov'], '-k', label='AQUAnoice with $Q^0$')
ax.plot(mons, 1e2*qflux['aclcov'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'Cloud cover (%)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface temperature
plotname = 'tsurf'
fig, ax = plt.subplots()
ax.plot(mons, ice['tsurf'], '-', color='tab:blue', label='AQUAice')
# ax.plot(mons, noice['tsurf'], '-k', label='AQUAnoice')
ax.plot(mons, qflux['tsurf'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$T_{s}$ (K)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface temperature
plotname = 'tsurf_cl'
fig, ax = plt.subplots()
ax.plot(mons, ice['tsurf'], '-', color='tab:blue', label='AQUAice')
ax.plot(mons, noice['tsurf'], '-k', label='AQUAnoice')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$T_{s}$ (K)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface temperature
plotname = 'temp2'
fig, ax = plt.subplots()
ax.plot(mons, ice['temp2'], '-', color='tab:blue', label='AQUAice')
ax.plot(mons, noice['temp2'], '-k', label='AQUAnoice')
ax.plot(mons, qflux['temp2'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$T_{2\,m}$ (K)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot temp tendency
dt = 30*86400

# centepurple
# ice['ttend'] = (ice['tsurf'][2:] - ice['tsurf'][:-2])/(2*dt)
# ice['ttend'] = np.insert(ice['ttend'], 0, (ice['tsurf'][1] - ice['tsurf'][-1])/(2*dt) )
# ice['ttend'] = np.append(ice['ttend'], (ice['tsurf'][0] - ice['tsurf'][-2])/(2*dt))

# # forward
# ice['ttend'] = (ice['tsurf'][1:] - ice['tsurf'][:-1])/(dt)
# ice['ttend'] = np.append(ice['ttend'], (ice['tsurf'][0] - ice['tsurf'][-1])/(dt))

# noice['ttend'] = (noice['tsurf'][2:] - noice['tsurf'][:-2])/(2*dt)
# noice['ttend'] = np.insert(noice['ttend'], 0, (noice['tsurf'][1] - noice['tsurf'][-1])/(2*dt) )
# noice['ttend'] = np.append(noice['ttend'], (noice['tsurf'][0] - noice['tsurf'][-2])/(2*dt))

# qflux['ttend'] = (qflux['tsurf'][2:] - qflux['tsurf'][:-2])/(2*dt)
# qflux['ttend'] = np.insert(qflux['ttend'], 0, (qflux['tsurf'][1] - qflux['tsurf'][-1])/(2*dt) )
# qflux['ttend'] = np.append(qflux['ttend'], (qflux['tsurf'][0] - qflux['tsurf'][-2])/(2*dt))

plotname = 'ttend'
fig, ax = plt.subplots()
ax.plot(mons, ice['ttend'], '-', color='tab:blue', label='AQUAice')
ax.plot(mons, noice['ttend'], '-k', label='AQUAnoice')
ax.plot(mons, qflux['ttend'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$\partial_t T_{s}$ (K s$^{-1}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot ice thickness tendency
dt = 30*86400

# ice['sntend'] = (ice['sni'][2:] - ice['sni'][:-2])/(2*dt)
# ice['sntend'] = np.insert(ice['sntend'], 0, (ice['sni'][1] - ice['sni'][-1])/(2*dt) )
# ice['sntend'] = np.append(ice['sntend'], (ice['sni'][0] - ice['sni'][-2])/(2*dt))

# # centepurple
# ice['stend'] = (ice['siced'][2:] - ice['siced'][:-2])/(2*dt)
# ice['stend'] = np.insert(ice['stend'], 0, (ice['siced'][1] - ice['siced'][-1])/(2*dt) )
# ice['stend'] = np.append(ice['stend'], (ice['siced'][0] - ice['siced'][-2])/(2*dt))

# # forward
# ice['stend'] = (ice['siced'][1:] - ice['siced'][:-1])/(dt)
# ice['stend'] = np.append(ice['stend'], (ice['siced'][0] - ice['siced'][-1])/(dt))

# ice['twtend'] = (ice['tsw'][2:] - ice['tsw'][:-2])/(2*dt)
# ice['twtend'] = np.insert(ice['twtend'], 0, (ice['tsw'][1] - ice['tsw'][-1])/(2*dt) )
# ice['twtend'] = np.append(ice['twtend'], (ice['tsw'][0] - ice['tsw'][-2])/(2*dt))

# ice['titend'] = (ice['tsi'][2:] - ice['tsi'][:-2])/(2*dt)
# ice['titend'] = np.insert(ice['titend'], 0, (ice['tsi'][1] - ice['tsi'][-1])/(2*dt) )
# ice['titend'] = np.append(ice['titend'], (ice['tsi'][0] - ice['tsi'][-2])/(2*dt))

# noice['stend'] = (noice['siced'][2:] - noice['siced'][:-2])/(2*dt)
# noice['stend'] = np.insert(noice['stend'], 0, (noice['siced'][1] - noice['siced'][-1])/(2*dt) )
# noice['stend'] = np.append(noice['stend'], (noice['siced'][0] - noice['siced'][-2])/(2*dt))

# qflux['stend'] = (qflux['siced'][2:] - qflux['siced'][:-2])/(2*dt)
# qflux['stend'] = np.insert(qflux['stend'], 0, (qflux['siced'][1] - qflux['siced'][-1])/(2*dt) )
# qflux['stend'] = np.append(qflux['stend'], (qflux['siced'][0] - qflux['siced'][-2])/(2*dt))

plotname = 'stend'
fig, ax = plt.subplots()
ax.plot(mons, ice['stend'], '-', color='tab:blue', label='AQUAice')
ax.plot(mons, noice['stend'], '-k', label='AQUAnoice')
ax.plot(mons, qflux['stend'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$\partial_t h_{i}$ (m s$^{-1}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

plotname = 'rlstend'
fig, ax = plt.subplots()
ax.plot(mons, rhoi*Lf*ice['stend'], '-', color='tab:blue', label='AQUAice')
ax.plot(mons, rhoi*Lf*noice['stend'], '-k', label='AQUAnoice')
ax.plot(mons, rhoi*Lf*qflux['stend'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$\partial_t h_{i}$ (m s$^{-1}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot inferpurple mixed layer depth
cw = 4182; rhow = 997;
ice['d'] = (ice['fsfc'] + Lf*rhoi*ice['stend']) / (cw*rhow*ice['ttend'])
ice['d'][ice['stend']<0] = np.nan # depth is infinite where ice is melting
noice['d'] = (noice['fsfc'] + Lf*rhoi*noice['stend']) / (cw*rhow*noice['ttend'])
qflux['d'] = (qflux['fsfc'] + Lf*rhoi*qflux['stend']) / (cw*rhow*qflux['ttend'])

plotname = 'd'
fig, ax = plt.subplots()
ax.plot(mons, ice['d'], '-', color='tab:blue', label='AQUAice')
ax.plot(mons, noice['d'], '-k', label='AQUAnoice')
# ax.plot(mons, qflux['d'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$d_{\mathrm{inferpurple}}$ (m)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot seaice in aquaice
plotname = 'seaice'
fig, ax = plt.subplots()
# ax.axhline(0, linewidth=0.5, color='k')
ax.plot(mons, ice['seaice'], '-', color='k')
ax.set_title(r'Control AQUAice, $80-90^\circ$N')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'Sea ice fraction (unitless)')
fig.set_size_inches(4,3)
plt.tight_layout()
# plt.legend(loc='upper right', prop={'size':6}, frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot inferpurple net FSFCice
d = 40
di = 0.05
ci = 2106
cs = 2090
rhos = 300
rhoh2o = 1000
dinf=0.25
ki = 2.1656 # ice thermal conductivity [W/K/m]
ks = 0.31 # snow thermal conductivity [W/K/m]
T0 = 273.15-1.8 # sea ice freezing temperature
# ice['istor']=rhoi*ci*ice['siced']*ice['titend']
ceff=(rhoi*ci*di+rhoh2o*cs*ice['sni'])
ice['istor']=ceff*ice['titend']
print('Max Cond (output) is %g W m**-2' % np.max(ice['ahfcon']))
plotname = 'fsfci_budget_diag'
fig, ax = plt.subplots()
ax.axhline(0, linewidth=0.5, color='k')
ax.plot(mons, ice['fsfc'], '-', color='k', label=r'$F_{SFC}$')
ax.plot(mons, ice['ahfcon'], '-', color='tab:purple', label='Conduction')
ax.plot(mons, ice['ahfres'], '-', color='tab:blue', label='Melt')
ax.plot(mons, ice['istor'], '-', color='tab:green', label='Ice Storage')
ax.plot(mons, ice['fsfc']-ice['ahfcon']-ice['ahfres']-ice['istor'], '--', color='k', label='Residual')
# ax.plot(mons, ice['istor'], '-', color='tab:green', label='Ice Storage')
# ax.plot(mons, ice['wstor'], '--', color='tab:green', label='Water Storage')
# ax.plot(mons, -Lf*rhoi*ice['stend'], '-', color='tab:gray', label='Latent energy')
# ax.plot(mons, -Lf*rhos*ice['sntend'], '--', color='tab:gray', label='Snow Enthalpy')
# ax.plot(mons, cw*rhow*dinf*ice['ttend'] - Lf*rhoi*ice['stend'], '--', color='tab:blue', label='inferpurple')
# ax.plot(mons, qflux['d'], '-', color='tab:purple', label=qlabel)
ax.set_title(r'Control AQUAice, $80-90^\circ$N')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'Energy flux (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(loc='upper right', prop={'size':6}, frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# # plot inferpurple net FSFCice
# # heff=ice['siced']+ki*rhoi/(ks*rhos)*ice['sn']
# heff=ice['siced']+ki*rhoh2o/(ks*rhos)*ice['sni']
# cond=ki*(ice['tsi']-T0)/heff
# ice['wstor']=rhow*cw*d*ice['twtend']
# le=ice['smelt']+ice['snmelt']
# # le=-Lf*rhoi*ice['smelt']-Lf*rhoh2o*ice['snmelt']
# # le=-Lf*rhoh2o*ice['sntend']
# sn=Lf*ice['aprs']
# # le=-Lf*rhoi*ice['stend']
# print('Max Cond (diagnosed) is %g W m**-2' % np.max(cond))
# plotname = 'fsfci_budget'
# fig, ax = plt.subplots()
# ax.axhline(0, linewidth=0.5, color='k')
# ax.plot(mons, ice['fsfc'], '-', color='k', label=r'$F_{SFC}$')
# ax.plot(mons, cond, '--', color='tab:purple', label='Conduction')
# ax.plot(mons, ice['istor'], '-', color='tab:green', label='Sensible Heat Storage')
# ax.plot(mons, ice['tmelt'], '--', color='tab:blue', label='Melt')
# # ax.plot(mons, le, '-', color='tab:gray', label='Latent energy')
# # ax.plot(mons, sn, '-', color='tab:blue', label='Snowfall')
# # ax.plot(mons, -Lf*rhoh2o*ice['sntend'], '--', color='tab:gray', label='Snow Enthalpy')
# ax.plot(mons, ice['fsfc']-cond-ice['tmelt']-ice['istor'], '--', color='k', label='Residual')
# # ax.plot(mons, cw*rhow*dinf*ice['ttend'] - Lf*rhoi*ice['stend'], '--', color='tab:blue', label='inferpurple')
# # ax.plot(mons, qflux['d'], '-', color='tab:purple', label=qlabel)
# ax.set_title(r'Control AQUAice, $80-90^\circ$N')
# ax.set_xlim([0,11])
# ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
# ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
# ax.yaxis.set_minor_locator(AutoMinorLocator())
# ax.set_ylabel(r'Energy flux (W m$^{-2}$)')
# fig.set_size_inches(4,3)
# plt.tight_layout()
# plt.legend(loc='upper right', prop={'size':6}, frameon=False)
# plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
# plt.close()

# # compare output and diagnosed melt
# plotname = 'decomp_melt'
# fig, ax = plt.subplots()
# ax.axhline(0, linewidth=0.5, color='k')
# ax.plot(mons, ice['ahfres'], '-', color='tab:blue', label='Melt (output)')
# ax.plot(mons, ice['tmelt'], '--', color='tab:blue', label='Melt (diagnosed)')
# # ax.plot(mons, ice['smelt']+ice['snmelt'], '-', color='tab:red', label='Melt (total)')
# ax.plot(mons, ice['smelt'], ':', color='tab:blue', label='Melt (bare sea ice)')
# ax.plot(mons, ice['snmelt'], '-.', color='tab:blue', label='Melt (snow)')
# # ax.plot(mons, Lf*ice['aprs'], '-.', color='tab:gray', label='Snowfall')
# ax.set_title(r'Control AQUAice, $80-90^\circ$N')
# ax.set_xlim([0,11])
# ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
# ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
# ax.yaxis.set_minor_locator(AutoMinorLocator())
# ax.set_ylabel(r'Energy flux (W m$^{-2}$)')
# fig.set_size_inches(4,3)
# plt.tight_layout()
# plt.legend(loc='upper right', prop={'size':6}, frameon=False)
# plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
# plt.close()

# # compare output and diagnosed melt
# plotname = 'diag_melt'
# fig, ax = plt.subplots()
# ax.axhline(0, linewidth=0.5, color='k')
# ax.plot(mons, ice['ahfres'], '-', color='tab:blue', label='Melt (output)')
# ax.plot(mons, le, '--', color='tab:blue', label='Melt (diagnosed)')
# ax.plot(mons, ice['tmelt'], '--', color='tab:red', label='Melt (t, diagnosed)')
# ax.set_title(r'Control AQUAice, $80-90^\circ$N')
# ax.set_xlim([0,11])
# ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
# ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
# ax.yaxis.set_minor_locator(AutoMinorLocator())
# ax.set_ylabel(r'Energy flux (W m$^{-2}$)')
# fig.set_size_inches(4,3)
# plt.tight_layout()
# plt.legend(loc='upper right', prop={'size':6}, frameon=False)
# plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
# plt.close()

# # compare output and diagnosed conduction
# plotname = 'diag_cond'
# fig, ax = plt.subplots()
# ax.axhline(0, linewidth=0.5, color='k')
# ax.plot(mons, ice['ahfcon'], '-', color='tab:purple', label='Conduction (output)')
# ax.plot(mons, ice['cond'], '--', color='tab:purple', label='Conduction (diagnosed)')
# ax.set_title(r'Control AQUAice, $80-90^\circ$N')
# ax.set_xlim([0,11])
# ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
# ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
# ax.yaxis.set_minor_locator(AutoMinorLocator())
# ax.set_ylabel(r'Energy flux (W m$^{-2}$)')
# fig.set_size_inches(4,3)
# plt.tight_layout()
# plt.legend(loc='lower left', prop={'size':6}, frameon=False)
# plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
# plt.close()

# plot inferpurple net FSFCwater
plotname = 'fsfcw_infer'
fig, ax = plt.subplots()
ax.plot(mons, noice['fsfc'], '-', color='k', label='actual')
ax.plot(mons, cw*rhow*40*noice['ttend'], '--', color='k', label='inferpurple')
# ax.plot(mons, qflux['d'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$F_{SFC,w}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface sw radiation
plotname = 'srads'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['srads'], '-', color='tab:blue', label='AQUAice')
ax.plot(mons, noice['srads'], '-k', label='AQUAnoice')
ax.plot(mons, qflux['srads'], '-', color='tab:purple', label=qlabel)
print('%s wm**-2' % np.mean(qflux['srads']-ice['srads']))
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot net surface fluxes
plotname = 'fsfc'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['fsfc'], '-', color='tab:blue', label='AQUAice')
# ax.plot(mons, noice['fsfc'], '-k', label='AQUAnoice')
ax.plot(mons, qflux['fsfc'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$F_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot net surface fluxes
plotname = 'fsfc_cl'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['fsfc'], '-', color='tab:blue', label='AQUAice')
ax.plot(mons, noice['fsfc'], '-k', label='AQUAnoice')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$F_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot net surface fluxes
plotname = 'fsfc_mq'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['fsfc'], '-', color='tab:blue', label='AQUAice')
# ax.plot(mons, noice['fsfc'], '-k', label='AQUAnoice')
ax.plot(mons, qflux['fsfc']-qflux['qflux'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$F_{SFC}-Q_{SW}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface latent heat flux 
plotname = 'ahfl'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['ahfl'], '-', color='tab:blue', label='AQUAice')
# ax.plot(mons, noice['ahfl'], '-k', label='AQUAnoice')
ax.plot(mons, qflux['ahfl'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$LH$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface latent heat flux 
plotname = 'ahfl_cl'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['ahfl'], '-', color='tab:blue', label='AQUAice')
ax.plot(mons, noice['ahfl'], '-k', label='AQUAnoice')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$LH$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface sensible heat flux
plotname = 'ahfs'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['ahfs'], '-', color='tab:blue', label='AQUAice')
# ax.plot(mons, noice['ahfs'], '-k', label='AQUAnoice')
ax.plot(mons, qflux['ahfs'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SH$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface sensible heat flux
plotname = 'ahfs_cl'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['ahfs'], '-', color='tab:blue', label='AQUAice')
ax.plot(mons, noice['ahfs'], '-k', label='AQUAnoice')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SH$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface LW radiation
plotname = 'trads'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['trads'], '-', color='tab:blue', label='AQUAice')
# ax.plot(mons, noice['trads'], '-k', label='AQUAnoice')
ax.plot(mons, qflux['trads'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$LW_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface LW radiation
plotname = 'trads_cl'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['trads'], '-', color='tab:blue', label='AQUAice')
ax.plot(mons, noice['trads'], '-k', label='AQUAnoice')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$LW_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface sw radiation
plotname = 'srads'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['srads'], '-', color='tab:blue', label='AQUAice')
# ax.plot(mons, noice['srads'], '-k', label='AQUAnoice')
ax.plot(mons, qflux['srads'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface sw radiation
plotname = 'srads_cl'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['srads'], '-', color='tab:blue', label='AQUAice')
ax.plot(mons, noice['srads'], '-k', label='AQUAnoice')
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface sw radiation
plotname = 'srads_mq'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['srads'], '-', color='tab:blue', label='AQUAice')
# ax.plot(mons, noice['srads'], '-k', label='AQUAnoice')
ax.plot(mons, qflux['srads']-qflux['qflux'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW_{SFC}-Q_{SW}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot downward surface sw
ice['sradsd'] = ice['srads'] - ice['sradsu']
noice['sradsd'] = noice['srads'] - noice['sradsu']
qflux['sradsd'] = qflux['srads'] - qflux['sradsu']
plotname = 'sradsd'
fig, ax = plt.subplots()
ax.plot(mons, ice['sradsd'], '-', color='tab:blue', label='AQUAice')
ax.plot(mons, noice['sradsd'], '-k', label='AQUAnoice')
ax.plot(mons, qflux['sradsd'], '-', color='tab:purple', label=qlabel)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW^{\downarrow}_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

