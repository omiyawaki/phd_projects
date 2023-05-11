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
varnames = ['tradsu','aclcov','twtend', 'titend', 'ttend', 'stend', 'ahfres', 'tsi', 'tsw', 'ahfcon', 'srafs', 'srad0', 'srad0d', 'albedo', 'siced', 'temp2', 'tsurf', 'fsfc', 'trads', 'srads', 'sradsu', 'ahfl', 'ahfs']

ice_label='AQUAice'
noice_label='AQUAnoice with $Q^0$'

keepvars=['siced', 'tsi', 'tsw', 'tsurf', 'srad0', 'srad0u', 'srad0d', 'srads', 'sradsu', 'albedo']

q1sim='rp000190d'
q1_label=r'AQUAnoice with $Q^1$'
q2sim='rp000190e'
q2_label=r'AQUAnoice with $Q^2$'
q3sim='rp000190f'
q3_label=r'AQUAnoice with $Q^3$'

plotdir='./plot/q0-q3'
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
    dsvar = dsvar.mean(dim='lon') # take zonal mean
    # take annual mean
    dsvar = dsvar.mean('time')
    ice[varname] = dsvar.data

noice = {}; noicek={};
run = 'rp000190b'
for varname in varnames:
    ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
    dsvar = getattr(ds, varname)
    if varname in keepvars:
        noicek[varname] = dsvar.data
    dsvar = dsvar.mean(dim='lon') # take zonal mean
    # take annual mean
    dsvar = dsvar.mean('time')
    noice[varname] = dsvar.data

# also load qflux
varname='qflux'
ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
dsvar = getattr(ds, 'aflux')
if varname in keepvars:
    noicek[varname] = dsvar.data
dsvar = dsvar.mean(dim='lon') # take zonal mean
# take annual mean
dsvar = dsvar.mean('time')
noice[varname] = dsvar.data

q1 = {}; q1k = {};
run = q1sim
for varname in varnames:
    ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
    dsvar = getattr(ds, varname)
    if varname in keepvars:
        q1k[varname] = dsvar.data
    dsvar = dsvar.mean(dim='lon') # take zonal mean
    # take annual mean
    dsvar = dsvar.mean('time')
    q1[varname] = dsvar.data

# also load qflux
varname='qflux'
ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
dsvar = getattr(ds, 'aflux')
if varname in keepvars:
    q1k[varname] = dsvar.data
dsvar = dsvar.mean(dim='lon') # take zonal mean
# take annual mean
dsvar = dsvar.mean('time')
q1[varname] = dsvar.data

q2 = {}; q2k = {};
run = q2sim
for varname in varnames:
    ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
    dsvar = getattr(ds, varname)
    if varname in keepvars:
        q2k[varname] = dsvar.data
    dsvar = dsvar.mean(dim='lon') # take zonal mean
    # take annual mean
    dsvar = dsvar.mean('time')
    q2[varname] = dsvar.data

# also load qflux
varname='qflux'
ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
dsvar = getattr(ds, 'aflux')
if varname in keepvars:
    q2k[varname] = dsvar.data
dsvar = dsvar.mean(dim='lon') # take zonal mean
# take annual mean
dsvar = dsvar.mean('time')
q2[varname] = dsvar.data

q3 = {}; q3k = {};
run = q3sim
for varname in varnames:
    ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
    dsvar = getattr(ds, varname)
    if varname in keepvars:
        q3k[varname] = dsvar.data
    dsvar = dsvar.mean(dim='lon') # take zonal mean
    # take annual mean
    dsvar = dsvar.mean('time')
    q3[varname] = dsvar.data

# also load qflux
varname='qflux'
ds = xr.open_dataset('%s/%s/%s_%s_0020_0039.ymonmean-20.nc' % (prefix, run, varname, run))
dsvar = getattr(ds, 'aflux')
if varname in keepvars:
    q3k[varname] = dsvar.data
dsvar = dsvar.mean(dim='lon') # take zonal mean
# take annual mean
dsvar = dsvar.mean('time')
q3[varname] = dsvar.data

##################################################
# PLOTS
##################################################
lat = ds.lat

# plot total cloud cover
plotname = 'aclcov'
fig, ax = plt.subplots()
ax.plot(lat, ice['aclcov'], '-', color='tab:purple', label=ice_label)
ax.plot(lat, noice['aclcov'], '-', color='tab:blue', label=noice_label)
ax.plot(lat, q1['aclcov'], '-', color='tab:orange', label=q1_label)
ax.plot(lat, q2['aclcov'], '-', color='tab:green', label=q2_label)
ax.plot(lat, q3['aclcov'], '-', color='tab:red', label=q3_label)
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
ax.plot(lat, ice['tsurf'], '-', color='tab:purple', label=ice_label)
ax.plot(lat, noice['tsurf'], '-', color='tab:blue', label=noice_label)
ax.plot(lat, q1['tsurf'], '-', color='tab:orange', label=q1_label)
ax.plot(lat, q2['tsurf'], '-', color='tab:green', label=q2_label)
ax.plot(lat, q3['tsurf'], '-', color='tab:red', label=q3_label)
ax.set_xlim([-90,90])
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.set_title('ANN')
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylabel(r'$T_{s}$ (K)')
fig.set_size_inches(4,3)
plt.tight_layout()
plt.legend(frameon=False)
plt.savefig('./%s/%s.pdf' % (plotdir, plotname) )
plt.close()

# plot 2m temperature
plotname = 'temp2'
fig, ax = plt.subplots()
ax.plot(lat, ice['temp2'], '-', color='tab:purple', label=ice_label)
ax.plot(lat, noice['temp2'], '-', color='tab:blue', label=noice_label)
ax.plot(lat, q1['temp2'], '-', color='tab:orange', label=q1_label)
ax.plot(lat, q2['temp2'], '-', color='tab:green', label=q2_label)
ax.plot(lat, q3['temp2'], '-', color='tab:red', label=q3_label)
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
ax.plot(lat, ice['ttend'], '-', color='tab:purple', label=ice_label)
ax.plot(lat, noice['ttend'], '-', color='tab:blue', label=noice_label)
ax.plot(lat, q1['ttend'], '-', color='tab:orange', label=q1_label)
ax.plot(lat, q2['ttend'], '-', color='tab:green', label=q2_label)
ax.plot(lat, q3['ttend'], '-', color='tab:red', label=q3_label)
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

# plot surface sw radiation
plotname = 'srads'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(lat, ice['srads'], '-', color='tab:purple', label=ice_label)
ax.plot(lat, noice['srads'], '-', color='tab:blue', label=noice_label)
ax.plot(lat, q1['srads'], '-', color='tab:orange', label=q1_label)
ax.plot(lat, q2['srads'], '-', color='tab:green', label=q2_label)
ax.plot(lat, q3['srads'], '-', color='tab:red', label=q3_label)
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
ax.plot(lat, ice['fsfc'], '-', color='tab:purple', label=ice_label)
ax.plot(lat, noice['fsfc'], '-', color='tab:blue', label=noice_label)
ax.plot(lat, q1['fsfc'], '-', color='tab:orange', label=q1_label)
ax.plot(lat, q2['fsfc'], '-', color='tab:green', label=q2_label)
ax.plot(lat, q3['fsfc'], '-', color='tab:red', label=q3_label)
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

# plot surface latent heat flux 
plotname = 'ahfl'
fig, ax = plt.subplots()
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(lat, ice['ahfl'], '-', color='tab:purple', label=ice_label)
ax.plot(lat, noice['ahfl'], '-', color='tab:blue', label=noice_label)
ax.plot(lat, q1['ahfl'], '-', color='tab:orange', label=q1_label)
ax.plot(lat, q2['ahfl'], '-', color='tab:green', label=q2_label)
ax.plot(lat, q3['ahfl'], '-', color='tab:red', label=q3_label)
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
ax.plot(lat, ice['ahfs'], '-', color='tab:purple', label=ice_label)
ax.plot(lat, noice['ahfs'], '-', color='tab:blue', label=noice_label)
ax.plot(lat, q1['ahfs'], '-', color='tab:orange', label=q1_label)
ax.plot(lat, q2['ahfs'], '-', color='tab:green', label=q2_label)
ax.plot(lat, q3['ahfs'], '-', color='tab:red', label=q3_label)
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
ax.plot(lat, ice['ahfs'], '-', color='tab:purple', label=ice_label)
ax.plot(lat, noice['ahfs'], '-', color='tab:blue', label=noice_label)
ax.plot(lat, q1['ahfs'], '-', color='tab:orange', label=q1_label)
ax.plot(lat, q2['ahfs'], '-', color='tab:green', label=q2_label)
ax.plot(lat, q3['ahfs'], '-', color='tab:red', label=q3_label)
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

# plot net surface sw minus qflux
plotname = 'srads-qflux'
fig, ax = plt.subplots()
ax.plot(lat,   ice['srads'], '-', color='tab:purple', label=ice_label)
ax.plot(lat, noice['srads'] - noice['qflux'], '-', color='tab:blue', label=noice_label)
ax.plot(lat,    q1['srads'] -    q1['qflux'], '-', color='tab:orange', label=q1_label)
ax.plot(lat,    q2['srads'] -    q2['qflux'], '-', color='tab:green', label=q2_label)
ax.plot(lat,    q3['srads'] -    q3['qflux'], '-', color='tab:red', label=q3_label)
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
