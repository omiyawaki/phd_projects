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
varnames = ['albedo', 'siced', 'temp2', 'tsurf', 'fsfc', 'trads', 'srads', 'sradsu', 'ahfl', 'ahfs']

qfsim='rp000190f'
plotdir='./paper_plots/%s' % (qfsim) 
if not os.path.exists(plotdir):
    os.mkdir(plotdir)

ice_label='AQUAice'
dq_label=r'AQUAnoice with $Q_C$'
qflux_label=r'AQUAnoice'

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

dq = {}
run = 'rp000190f'
for varname in varnames:
    ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
    dsvar = getattr(ds, varname)
    dsvar = dsvar.mean(dim='lon') # take zonal mean
    # take area mean over high latitudes
    dsvar = dsvar.where(dsvar.lat >= latarc, drop=True)
    w = np.cos(np.deg2rad(dsvar.lat))
    dsvarw = dsvar.weighted(w)
    dsvar = dsvarw.mean(dim='lat')
    dq[varname] = dsvar.data

qflux = {}
run = qfsim
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
ax.plot(mons, ice['tsurf'], '-', color='tab:blue', label=ice_label)
# ax.plot(mons, dq['tsurf'], '-', color='tab:red', label=dq_label)
ax.plot(mons, qflux['tsurf'], '-', color='tab:purple', label=qflux_label)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$T_{s}$ (K)')
fig.set_size_inches(4,3)
plt.tight_layout()
# plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface temperature
plotname = 'temp2'
fig, ax = plt.subplots()
ax.plot(mons, ice['temp2'], '-', color='tab:blue', label=ice_label)
# ax.plot(mons, dq['temp2'], '-', color='tab:red', label=dq_label)
ax.plot(mons, qflux['temp2'], '-', color='tab:purple', label=qflux_label)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$T_{2\,m}$ (K)')
fig.set_size_inches(4,3)
plt.tight_layout()
# plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot temp tendency
dt = 30*86400

ice['ttend'] = (ice['tsurf'][2:] - ice['tsurf'][:-2])/(2*dt)
ice['ttend'] = np.insert(ice['ttend'], 0, (ice['tsurf'][1] - ice['tsurf'][-1])/(2*dt) )
ice['ttend'] = np.append(ice['ttend'], (ice['tsurf'][0] - ice['tsurf'][-2])/(2*dt))

dq['ttend'] = (dq['tsurf'][2:] - dq['tsurf'][:-2])/(2*dt)
dq['ttend'] = np.insert(dq['ttend'], 0, (dq['tsurf'][1] - dq['tsurf'][-1])/(2*dt) )
dq['ttend'] = np.append(dq['ttend'], (dq['tsurf'][0] - dq['tsurf'][-2])/(2*dt))

qflux['ttend'] = (qflux['tsurf'][2:] - qflux['tsurf'][:-2])/(2*dt)
qflux['ttend'] = np.insert(qflux['ttend'], 0, (qflux['tsurf'][1] - qflux['tsurf'][-1])/(2*dt) )
qflux['ttend'] = np.append(qflux['ttend'], (qflux['tsurf'][0] - qflux['tsurf'][-2])/(2*dt))

plotname = 'ttend'
fig, ax = plt.subplots()
ax.plot(mons, ice['ttend'], '-', color='tab:blue', label=ice_label)
# ax.plot(mons, dq['ttend'], '-', color='tab:red', label=dq_label)
ax.plot(mons, qflux['ttend'], '-', color='tab:purple', label=qflux_label)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$\partial_t T_{s}$ (K s$^{-1}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
# plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface sw radiation
plotname = 'srads'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['srads'], '-', color='tab:blue', label=ice_label)
# ax.plot(mons, dq['srads'], '-', color='tab:red', label=dq_label)
ax.plot(mons, qflux['srads'], '-', color='tab:purple', label=qflux_label)
print('%s wm**-2' % np.mean(qflux['srads']-ice['srads']))
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
# plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot net surface fluxes
plotname = 'fsfc'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['fsfc'], '-', color='tab:blue', label=ice_label)
# ax.plot(mons, dq['fsfc'], '-', color='tab:red', label=dq_label)
ax.plot(mons, qflux['fsfc'], '-', color='tab:purple', label=qflux_label)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$F_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
# plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface latent heat flux 
plotname = 'ahfl'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['ahfl'], '-', color='tab:blue', label=ice_label)
# ax.plot(mons, dq['ahfl'], '-', color='tab:red', label=dq_label)
ax.plot(mons, qflux['ahfl'], '-', color='tab:purple', label=qflux_label)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$LH$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
# plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface sensible heat flux
plotname = 'ahfs'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['ahfs'], '-', color='tab:blue', label=ice_label)
# ax.plot(mons, dq['ahfs'], '-', color='tab:red', label=dq_label)
ax.plot(mons, qflux['ahfs'], '-', color='tab:purple', label=qflux_label)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SH$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
# plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface LW radiation
plotname = 'trads'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(mons, ice['trads'], '-', color='tab:blue', label=ice_label)
# ax.plot(mons, dq['trads'], '-', color='tab:red', label=dq_label)
ax.plot(mons, qflux['trads'], '-', color='tab:purple', label=qflux_label)
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
ax.plot(mons, ice['srads'], '-', color='tab:blue', label=ice_label)
# ax.plot(mons, dq['srads'], '-', color='tab:red', label=dq_label)
ax.plot(mons, qflux['srads'], '-', color='tab:purple', label=qflux_label)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
# plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot downward surface sw
ice['sradsd'] = ice['srads'] - ice['sradsu']
dq['sradsd'] = dq['srads'] - dq['sradsu']
qflux['sradsd'] = qflux['srads'] - qflux['sradsu']
plotname = 'sradsd'
fig, ax = plt.subplots()
ax.plot(mons, ice['sradsd'], '-', color='tab:blue', label=ice_label)
# ax.plot(mons, dq['sradsd'], '-', color='tab:red', label=dq_label)
ax.plot(mons, qflux['sradsd'], '-', color='tab:purple', label=qflux_label)
ax.set_xlim([0,11])
ax.set_xticks(mons, ['J','F','M','A','M','J','J','A','S','O','N','D'])
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW^{\downarrow}_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
# plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()
