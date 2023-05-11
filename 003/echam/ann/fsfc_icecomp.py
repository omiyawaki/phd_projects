import os
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

latarc = 80 # equatorward latitude of Arctic domain
Lf = 3.337e5 # J kg**-1 latent heat of fusion
rhoi = 917 # kg m**-3 density of ice

prefix = '/project2/tas1/miyawaki/projects/003/data/raw/echam'
# varnames = ['srafs', 'srad0', 'srad0u', 'srad0d', 'albedo', 'siced', 'temp2', 'tsurf', 'fsfc', 'trads', 'srads', 'sradsu', 'ahfl', 'ahfs']
varnames = ['swrefu','tradsu','aclcov','twtend', 'titend', 'ttend', 'stend', 'ahfres', 'tsi', 'tsw', 'ahfcon', 'srafs', 'srad0', 'srad0d', 'albedo', 'siced', 'temp2', 'tsurf', 'fsfc', 'trads', 'srads', 'sradsu', 'ahfl', 'ahfs']

ice_label='AQUAice'
# noice_label='AQUAnoice with $Q^0$'
noice_label='AQUAnoice with no $Q$'

keepvars=['siced', 'tsi', 'tsw', 'tsurf', 'srad0', 'srad0u', 'srad0d', 'srads', 'sradsu', 'albedo']
qfsim='rp000190f'
# q_label=r'AQUAqflux'
# q_label=r'AQUAnoice with $Q_C+Q_{SW_\alpha}$'
q_label=r'AQUAnoice with $Q^3$'
print(qfsim)
plotdir='./plot/%s' % (qfsim) 
if not os.path.exists(plotdir):
    os.makedirs(plotdir)

ice = {}; icek={};
run = 'rp000134'
for varname in varnames:
    ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
    dsvar = getattr(ds, varname)
    # keep original albedo and srads
    if varname in keepvars:
        icek[varname] = dsvar.data
    # dsvar = dsvar.mean(dim='lon',skipna=True) # take zonal mean
    # # take annual mean
    # dsvar = dsvar.mean('time',skipna=True)
    ice[varname] = dsvar.data
    ice[varname] = dsvar.data
    ice[varname][ice[varname]>1e10]=np.nan
    ice[varname]=np.nanmean(ice[varname],axis=(0,2))

noice = {}; noicek={};
run = 'rp000135'
for varname in varnames:
    ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
    dsvar = getattr(ds, varname)
    if varname in keepvars:
        noicek[varname] = dsvar.data
    # dsvar = dsvar.mean(dim='lon',skipna=True) # take zonal mean
    # take annual mean
    # dsvar = dsvar.mean('time',skipna=True)
    noice[varname] = dsvar.data
    noice[varname][noice[varname]>1e10]=np.nan
    noice[varname]=np.nanmean(noice[varname],axis=(0,2))

# # compute planetary albedo
# icek['ap'] = -icek['srad0u']/icek['srad0d']
# noicek['ap'] = -noicek['srad0u']/noicek['srad0d']
# icek['ap'][np.where(icek['srad0d']==0)] = 0
# noicek['ap'][np.where(noicek['srad0d']==0)] = 0
# ice['ap'] = np.mean(icek['ap'], axis=(0,2))
# noice['ap'] = np.mean(noicek['ap'], axis=(0,2))

# # compute pswabs using ice sw down
# ice['pswabs'] = (1-icek['ap'])*icek['srad0d']
# noice['pswabs'] = (1-noicek['ap'])*icek['srad0d']
# ice['pswabs'] = np.mean(ice['pswabs'], axis=(0,2))
# noice['pswabs'] = np.mean(noice['pswabs'], axis=(0,2))

# # compute sswabs using ice sw down
# ice['sswabs'] = (1-icek['albedo'])*(icek['srads']-icek['sradsu'])
# noice['sswabs'] = (1-noicek['albedo'])*(icek['srads']-icek['sradsu'])
# ice['sswabs'] = np.mean(ice['sswabs'], axis=(0,2))
# noice['sswabs'] = np.mean(noice['sswabs'], axis=(0,2))

qflux = {}; qfluxk = {};
run = qfsim
for varname in varnames:
    # ds = xr.open_dataset('%s/%s/%s_%s_0001_0009.ymonmean-5.nc' % (prefix, run, varname, run))
    ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
    dsvar = getattr(ds, varname)
    if varname in keepvars:
        qfluxk[varname] = dsvar.data
    dsvar = dsvar.mean(dim='lon') # take zonal mean
    # take annual mean
    dsvar = dsvar.mean('time')
    qflux[varname] = dsvar.data

# also load qflux
varname='qflux'
# ds = xr.open_dataset('%s/%s/%s_%s_0001_0009.ymonmean-5.nc' % (prefix, run, varname, run))
ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
dsvar = getattr(ds, 'aflux')
if varname in keepvars:
    qfluxk[varname] = dsvar.data
dsvar = dsvar.mean(dim='lon') # take zonal mean
# take annual mean
dsvar = dsvar.mean('time')
qflux[varname] = dsvar.data

##################################################
# PLOTS
##################################################
lat = ds.lat

# plot upward sw reflectivity
plotname = 'swrefu'
fig, ax = plt.subplots()
print(noice['swrefu'])
ax.plot(lat, ice['swrefu'], '-', color='tab:purple', label=ice_label)
ax.plot(lat, noice['swrefu'], '-k', label=noice_label)
# ax.plot(lat, qflux['aclcov'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$R_{\uparrow}$ (unitless)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot total cloud cover
plotname = 'aclcov'
fig, ax = plt.subplots()
ax.plot(lat, ice['aclcov'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['aclcov'], '-k', label=noice_label)
ax.plot(lat, qflux['aclcov'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'Cloud fraction (unitless)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface temperature
plotname = 'tsurf'
fig, ax = plt.subplots()
ax.plot(lat, ice['tsurf'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['tsurf'], '-k', label=noice_label)
ax.plot(lat, qflux['tsurf'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
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
ax.plot(lat, ice['tsurf'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['tsurf'], '-k', label=noice_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
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
ax.plot(lat, ice['temp2'], '-', color='tab:purple', label=ice_label)
ax.plot(lat, noice['temp2'], '-', color='tab:blue', label=noice_label)
ax.plot(lat, qflux['temp2'], '-', color='tab:orange', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
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

plotname = 'ttend'
fig, ax = plt.subplots()
ax.plot(lat, ice['ttend'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['ttend'], '-k', label=noice_label)
ax.plot(lat, qflux['ttend'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
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

plotname = 'stend'
fig, ax = plt.subplots()
ax.plot(lat, ice['stend'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['stend'], '-k', label=noice_label)
ax.plot(lat, qflux['stend'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
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
ax.plot(lat, rhoi*Lf*ice['stend'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, rhoi*Lf*noice['stend'], '-k', label=noice_label)
ax.plot(lat, rhoi*Lf*qflux['stend'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$\partial_t h_{i}$ (m s$^{-1}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot inferred mixed layer depth
cw = 3994; rhow = 1025; d=40;
ice['d'] = (ice['fsfc'] + Lf*rhoi*ice['stend']) / (cw*rhow*ice['ttend'])
ice['d'][ice['stend']<0] = np.nan # depth is infinite where ice is melting
noice['d'] = (noice['fsfc'] + Lf*rhoi*noice['stend']) / (cw*rhow*noice['ttend'])
qflux['d'] = (qflux['fsfc'] + Lf*rhoi*qflux['stend']) / (cw*rhow*qflux['ttend'])

plotname = 'd'
fig, ax = plt.subplots()
ax.plot(lat, ice['d'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['d'], '-k', label=noice_label)
# ax.plot(lat, qflux['d'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$d_{\mathrm{inferred}}$ (m)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot inferred net FSFCice
ci = 2106
ice['istor']=rhoi*ci*ice['siced']*ice['titend']
ice['wstor']=rhow*cw*d*ice['twtend']
dinf=0.25
plotname = 'fsfci_budget'
fig, ax = plt.subplots()
ax.axhline(0, linewidth=0.5, color='k')
ax.plot(lat, ice['fsfc'], '-', color='k', label='$F_{SFC}$')
# ax.plot(lat, cw*rhow*dinf*ice['ttend'] - Lf*rhoi*ice['stend'], '--', color='tab:blue', label='inferred')
ax.plot(lat, ice['ahfcon'], '-', color='tab:purple', label=r'Conductive flux')
ax.plot(lat, ice['ahfres'], '-', color='tab:blue', label='Ice melt')
# ax.plot(lat, ice['istor'], '-', color='tab:green', label=r'Ice heat storage')
# ax.plot(lat, ice['wstor'], '--', color='tab:green', label=r'Sea heat storage')
# ax.plot(lat, -Lf*rhoi*ice['stend'], '-', color='tab:gray', label=r'$-L_f \rho_i \partial_t h$')
# res=ice['fsfc']-ice['ahfcon']-ice['ahfres']-ice['istor']-ice['wstor']+Lf*rhoi*ice['stend']
res=ice['fsfc']-ice['ahfcon']-ice['ahfres']
ax.plot(lat, res, '--', color='k', label='Residual')
ax.set_title(r'Control AQUAice, ANN')
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$F_{SFC,i}$ (W m$^{-2}$)')
fig.set_size_inches(5,3)
plt.tight_layout()
box=ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot inferred net FSFCwater
plotname = 'fsfcw_infer'
fig, ax = plt.subplots()
ax.plot(lat, noice['fsfc'], '-', color='k', label='actual')
ax.plot(lat, cw*rhow*40*noice['ttend'], '--', color='k', label='inferred')
# ax.plot(lat, qflux['d'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
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
ax.plot(lat, ice['srads'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['srads'], '-k', label=noice_label)
ax.plot(lat, qflux['srads'], '-', color='tab:purple', label=q_label)
print('%s wm**-2' % np.mean(qflux['srads']-ice['srads']))
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
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
ax.plot(lat, ice['fsfc'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['fsfc'], '-k', label=noice_label)
ax.plot(lat, qflux['fsfc'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
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
ax.plot(lat, ice['fsfc'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['fsfc'], '-k', label=noice_label)
# ax.plot(lat, qflux['fsfc'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
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
ax.plot(lat, ice['fsfc'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['fsfc'], '-k', label=noice_label)
ax.plot(lat, qflux['fsfc']-qflux['qflux'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
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
ax.plot(lat, ice['ahfl'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['ahfl'], '-k', label=noice_label)
ax.plot(lat, qflux['ahfl'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
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
ax.plot(lat, ice['ahfl'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['ahfl'], '-k', label=noice_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
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
ax.plot(lat, ice['ahfs'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['ahfs'], '-k', label=noice_label)
ax.plot(lat, qflux['ahfs'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
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
ax.plot(lat, ice['ahfs'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['ahfs'], '-k', label=noice_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
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
ax.plot(lat, ice['trads'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['trads'], '-k', label=noice_label)
ax.plot(lat, qflux['trads'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$LW_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface LW radiation up
plotname = 'tradsu'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(lat, ice['tradsu'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['tradsu'], '-k', label=noice_label)
ax.plot(lat, qflux['tradsu'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$LW^{\uparrow}_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface LW radiation
plotname = 'trads_cl'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(lat, ice['trads'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['trads'], '-k', label=noice_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
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
ax.plot(lat, ice['srads'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['srads'], '-k', label=noice_label)
ax.plot(lat, qflux['srads'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
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
ax.plot(lat, ice['srads'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['srads'], '-k', label=noice_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
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
ax.plot(lat, ice['srads'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['srads'], '-k', label=noice_label)
ax.plot(lat, qflux['srads']-qflux['qflux'], 'tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW_{SFC}-Q_{SW}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot net surface sw minus qflux
plotname = 'srads-qflux'
fig, ax = plt.subplots()
ax.plot(lat, ice['srads'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['srads'], '-k', label=noice_label)
ax.plot(lat, qflux['srads']-qflux['qflux'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW_{SFC} - Q$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot DIFF surface sw radiation
plotname = 'srads_diff'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
# ax.plot(lat, noice['srads']-ice['srads'], '-k', label='ni - i')
ax.plot(lat, np.mean(noicek['srads']-icek['srads'], axis=(0,2)), '-k', label='ni - i')
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$\Delta SW_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot DIFF surface sw radiation
plotname = 'srads_diffq'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
# ax.plot(lat, noice['srads']-ice['srads'], '-k', label='ni - i')
ax.plot(lat, np.mean(qfluxk['srads']-icek['srads'], axis=(0,2)), '-k', label='q - i')
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$\Delta SW_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot TOA sw radiation
plotname = 'srad0'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(lat, ice['srad0'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['srad0'], '-k', label=noice_label)
ax.plot(lat, qflux['srad0'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW_{TOA}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot DIFF TOA sw radiation
plotname = 'srad0_diff'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
# ax.plot(lat, noice['srad0']-ice['srad0'], '-k', label='ni - i')
ax.plot(lat, np.mean(noicek['srad0']-icek['srad0'], axis=(0,2)), '-k', label='ni - i')
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$\Delta SW_{TOA}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# # plot planetary albedo
# plotname = 'ap'
# fig, ax = plt.subplots()
# ax.axhline(0, color='k', linewidth=0.5)
# ax.plot(lat, ice['ap'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['ap'], '-k', label=noice_label)
# ax.set_xlim([-90,90])
# ax.set_xticks(np.arange(-90,91,30))
# ax.set_xlabel('Latitude (deg)')
# ax.xaxis.set_minor_locator(MultipleLocator(10))
# ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
# ax.yaxis.set_minor_locator(AutoMinorLocator())
# ax.set_ylabel(r'${SW}_p$ (unitless)')
# fig.set_size_inches(4,3)
# plt.tight_layout()
# plt.legend(frameon=False)
# plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
# plt.close()

# plot downward surface sw
ice['sradsd'] = ice['srads'] - ice['sradsu']
noice['sradsd'] = noice['srads'] - noice['sradsu']
qflux['sradsd'] = qflux['srads']- qflux['sradsu']
plotname = 'sradsd'
fig, ax = plt.subplots()
ax.plot(lat, ice['sradsd'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['sradsd'], '-k', label=noice_label)
ax.plot(lat, qflux['sradsd'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW^{\downarrow}_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot swabs
ice['swabs']=ice['srad0d']-ice['sradsd']
noice['swabs']=noice['srad0d']-noice['sradsd']
qflux['swabs']=qflux['srad0d']-qflux['sradsd']
plotname = 'swabs'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(lat, ice['swabs'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['swabs'], '-k', label=noice_label)
ax.plot(lat, qflux['swabs'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW^{\downarrow}_{TOA}-SW^{\downarrow}_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface LW radiation down
ice['tradsd'] = ice['trads'] - ice['tradsu']
noice['tradsd'] = noice['trads'] - noice['tradsu']
qflux['tradsd'] = qflux['trads'] - qflux['tradsu']
plotname = 'tradsd'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(lat, ice['tradsd'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['tradsd'], '-k', label=noice_label)
ax.plot(lat, qflux['tradsd'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$LW^{\downarrow}_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface net radiation down
plotname = 'radsd'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(lat, ice['sradsd']+ice['tradsd'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['radsd'], '-k', label=noice_label)
ax.plot(lat, qflux['sradsd']+qflux['tradsd'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$R^{\downarrow}_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface net radiation up
plotname = 'radsu'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(lat, ice['sradsu']+ice['tradsu'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['radsu'], '-k', label=noice_label)
ax.plot(lat, qflux['sradsu']+qflux['tradsu'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$R^{\uparrow}_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot surface net radiation
plotname = 'rads'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(lat, ice['srads']+ice['trads'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['radsu'], '-k', label=noice_label)
ax.plot(lat, qflux['srads']+qflux['trads'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$R_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot clear sky net surface sw
plotname = 'srafs'
fig, ax = plt.subplots()
ax.plot(lat, ice['srafs'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['srafs'], '-k', label=noice_label)
ax.plot(lat, qflux['srafs'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW_{SFC,\,clear}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot clear sky net surface sw minus qflux
plotname = 'srafs-qflux'
fig, ax = plt.subplots()
ax.plot(lat, ice['srafs'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['srafs'], '-k', label=noice_label)
ax.plot(lat, qflux['srafs']-qflux['qflux'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW_{SFC,\,clear} - Q$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot cloudy sky net surface sw
ice['sracs'] = ice['srads'] - ice['srafs']
noice['sracs'] = noice['srads'] - noice['srafs']
qflux['sracs'] = qflux['srads'] - qflux['srafs']
plotname = 'sracs'
fig, ax = plt.subplots()
ax.plot(lat, ice['sracs'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['sracs'], '-k', label=noice_label)
ax.plot(lat, qflux['sracs'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW_{SFC,\,cloudy}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot downward surface sw clear
ice['srafsd'] = ice['srafs']/(1-ice['albedo'])
noice['srafsd'] = noice['srafs']/(1-noice['albedo'])
qflux['srafsd'] = (qflux['srafs']-qflux['qflux'])/(1-qflux['albedo'])
plotname = 'srafsd'
fig, ax = plt.subplots()
ax.plot(lat, ice['srafsd'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['srafsd'], '-k', label=noice_label)
ax.plot(lat, qflux['srafsd'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW^{\downarrow}_{SFC,clear}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot swabs
ice['swabsc']=ice['srad0']-ice['srafsd']
noice['swabsc']=noice['srad0']-noice['srafsd']
qflux['swabsc']=qflux['srad0']-qflux['srafsd']
plotname = 'swabs_clear'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(lat, ice['swabsc'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['swabsc'], '-k', label=noice_label)
ax.plot(lat, qflux['swabsc'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW^{\downarrow}_{TOA}-SW^{\downarrow}_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot downward surface sw cloudy
ice['sracsd'] = ice['sracs']/(1-ice['albedo'])
noice['sracsd'] = noice['sracs']/(1-noice['albedo'])
qflux['sracsd'] = (qflux['sracs']-qflux['qflux'])/(1-qflux['albedo'])
plotname = 'sracsd'
fig, ax = plt.subplots()
ax.plot(lat, ice['sracsd'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['sracsd'], '-k', label=noice_label)
ax.plot(lat, qflux['sracsd'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW^{\downarrow}_{SFC,clear}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# # plot surf SW abs DIFF using same downward surface sw
# # ice['sswabs'] = (1-ice['albedo'])*(ice['srads'] - ice['sradsu'])
# # noice['sswabs'] = (1-noice['albedo'])*(ice['srads'] - ice['sradsu'])
# plotname = 'sswabs_diff'
# fig, ax = plt.subplots()
# ax.plot(lat, noice['sswabs']-ice['sswabs'], '-k', label='ni - i')
# ax.set_xlim([-90,90])
# ax.set_xticks(np.arange(-90,91,30))
# ax.set_xlabel('Latitude (deg)')
# ax.xaxis.set_minor_locator(MultipleLocator(10))
# ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
# ax.yaxis.set_minor_locator(AutoMinorLocator())
# ax.set_ylabel(r'$\Delta (1-{SW})SW^{\downarrow}_{SFC}$ (W m$^{-2}$)')
# fig.set_size_inches(4,3)
# plt.tight_layout()
# plt.legend(frameon=False)
# plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
# plt.close()

# # plot TOA SW abs DIFF using same downward TOA sw
# plotname = 'pswabs_diff'
# fig, ax = plt.subplots()
# ax.plot(lat, noice['pswabs']-ice['pswabs'], '-k', label='ni - i')
# ax.set_xlim([-90,90])
# ax.set_xticks(np.arange(-90,91,30))
# ax.set_xlabel('Latitude (deg)')
# ax.xaxis.set_minor_locator(MultipleLocator(10))
# ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
# ax.yaxis.set_minor_locator(AutoMinorLocator())
# ax.set_ylabel(r'$\Delta (1-{SW}_p)SW^{\downarrow}_{TOA}$ (W m$^{-2}$)')
# fig.set_size_inches(4,3)
# plt.tight_layout()
# plt.legend(frameon=False)
# plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
# plt.close()

# plot approx net shortwave
ice['asrads'] =(1-ice['albedo'])*ice['sradsd']
noice['asrads'] =(1-noice['albedo'])*noice['sradsd']
qflux['asrads'] =(1-qflux['albedo'])*qflux['sradsd']
plotname = 'srads_approx'
fig, ax = plt.subplots()
ax.plot(lat, ice['asrads'], '--', color='tab:blue', label=ice_label)
ax.plot(lat, noice['asrads'], '--k', label=noice_label)
ax.plot(lat, qflux['asrads'], '--', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$(1-{SW}_s)SW^{\downarrow}_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot approx net shortwave
plotname = 'srads_comp'
fig, ax = plt.subplots()
ax.plot(lat, ice['asrads'], '--', color='tab:blue', label=ice_label)
ax.plot(lat, noice['asrads'], '--k', label=noice_label)
ax.plot(lat, qflux['asrads'], '--', color='tab:purple', label=q_label)
ax.plot(lat, ice['srads'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['srads'], '-k', label=noice_label)
ax.plot(lat, qflux['srads'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$SW_{SFC}$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
# plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot inferred albedo
ice['infalb'] = 1+(ice['trads']+ice['ahfl']+ice['ahfs'])/ice['sradsd']
noice['infalb'] = 1+(noice['trads']+noice['ahfl']+noice['ahfs'])/noice['sradsd']
qflux['infalb'] = 1+(qflux['trads']+qflux['ahfl']+qflux['ahfs'])/qflux['sradsd']
plotname = 'infalb_comp'
fig, ax = plt.subplots()
ax.plot(lat, ice['infalb'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['infalb'], '-k', label=noice_label)
ax.plot(lat, qflux['infalb'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'Inferred ${SW}_{s}$ (unitless)')
fig.set_size_inches(4,3)
plt.tight_layout()
# plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot albedo
ice['albi'] = -ice['sradsu']/ice['sradsd']
noice['albi'] = -noice['sradsu']/noice['sradsd']
qflux['albi'] = -qflux['sradsu']/qflux['sradsd']
plotname = 'albedo_diff'
fig, ax = plt.subplots()
ax.plot(lat, ice['albedo']-ice['albi'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['albi'], '-k', label=noice_label)
ax.plot(lat, qflux['albedo']-qflux['albi'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$\alpha_{diag}-\frac{SW^{\uparrow}_{SFC}}{SW^{\downarrow}_{SFC}}$ (unitless)')
fig.set_size_inches(4,3)
plt.tight_layout()
# plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# # plot inferred SW optical thickness tau
# ice['swtau'] = np.log(ice['srad0d']/ice['sradsd'])
# noice['swtau'] = np.log(noice['srad0d']/noice['sradsd'])
# qflux['swtau'] = np.log(qflux['srad0d']/qflux['sradsd'])
# plotname = 'swtau'
# fig, ax = plt.subplots()
# ax.plot(lat, ice['swtau'], '-', color='tab:blue', label=ice_label)
# ax.plot(lat, noice['swtau'], '-k', label=noice_label)
# ax.plot(lat, qflux['swtau'], '-', color='tab:purple', label=q_label)
# ax.set_xlim([-90,90])
# ax.set_xticks(np.arange(-90,91,30))
# ax.set_xlabel('Latitude (deg)')
# ax.xaxis.set_minor_locator(MultipleLocator(10))
# ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
# ax.yaxis.set_minor_locator(AutoMinorLocator())
# ax.set_ylabel(r'$\tau_{SW}$ (unitless)')
# fig.set_size_inches(4,3)
# plt.tight_layout()
# plt.legend(frameon=False)
# plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
# plt.close()

# inferred q flux
ice['q'] = ice['fsfc'] - rhow*cw*d*ice['ttend']
noice['q'] = noice['fsfc'] - rhow*cw*d*noice['ttend']
qflux['q'] = qflux['fsfc'] - rhow*cw*d*qflux['ttend']
plotname = 'qflux'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(lat, ice['q'], '-', color='tab:blue', label=ice_label)
ax.plot(lat, noice['q'], '-k', label=noice_label)
ax.plot(lat, qflux['q'], '-', color='tab:purple', label=q_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$Q$ (W m$^{-2}$)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

print('\nice d_eff = %g m' % np.mean(ice['d']))
print('no ice d_eff = %g m' % np.mean(noice['d']))
print('amplification factor = %g' % (np.mean(noice['d'])/np.mean(ice['d'])) )
