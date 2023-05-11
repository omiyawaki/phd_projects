import os
import sys
import pickle
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from tools import get_modelidx, get_modelstr

nref=147
nrcp=295
model='mmm'
models=['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MPI-ESM-LR']

# specify directory for output data and plots
datadir = '%s/data' % ( os.getcwd() ) # output data directory
plotdir = '%s/plot/%s' % ( os.getcwd(), model ) # plot directory
if not os.path.exists(plotdir):
    os.mkdir(plotdir)

##################################################
# LOAD GCM DATA
##################################################
ds = xr.open_dataset('/project2/tas1/miyawaki/projects/003/echam/ghg/ghg_rcp85_1765-2500_c100203.nc')
co2 = 1e-6*ds.CO2.sel(time=slice('2006-01-01', '2300-12-31')).data

[mmm_clim, grid, indiv_clim] = pickle.load(open('../climlab/input_data/clima.pickle', 'rb'))
[mmm_forc, grid, indiv_forc] = pickle.load(open('../climlab/input_data/forcing.pickle', 'rb'))

gcm_clim = mmm_clim
gcm_forc = mmm_forc

for im in range(len(models)):
    gcm_clim_mmm = {}
    gcm_clim_mmm['prc25'] = {}
    gcm_clim_mmm['prc75'] = {}
    gcm_clim_mmm['std'] = {}
    gcm_forc_mmm = {}
    gcm_forc_mmm['prc25'] = {}
    gcm_forc_mmm['prc75'] = {}
    gcm_forc_mmm['std'] = {}
    for varname in indiv_clim:
        gcm_clim_mmm['prc25'][varname] = np.percentile(indiv_clim[varname], 25, axis=0)
        gcm_clim_mmm['prc75'][varname] = np.percentile(indiv_clim[varname], 75, axis=0)
        gcm_forc_mmm['prc25'][varname] = np.percentile(indiv_forc[varname], 25, axis=0)
        gcm_forc_mmm['prc75'][varname] = np.percentile(indiv_forc[varname], 75, axis=0)

        gcm_clim_mmm['std'][varname] = np.nanstd(indiv_clim[varname], axis=0)
        gcm_forc_mmm['std'][varname] = np.nanstd(indiv_forc[varname], axis=0)

##################################################
# LOAD RRTMG DATA
##################################################
# RRTMG full RCP8.5 response
ref = {}
[ref['ra'], ref['racs'], ref['diag']] = pickle.load(open('%s/rad.hist.%s.pickle' % (datadir, model), 'rb'))

ref_indiv = {}
ref_indiv['racs'] = np.empty([len(models), len(ref['racs'])])
for im in range(len(models)):
    [ra, racs, diag] = pickle.load(open('%s/rad.hist.%s.pickle' % (datadir, models[im]), 'rb'))
    ref_indiv['racs'][im,:] = racs
ref_mmm = {}
ref_mmm['prc25'] = {}
ref_mmm['prc75'] = {}
ref_mmm['std'] = {}
ref_mmm['prc25'] = np.percentile(ref_indiv['racs'], 25, axis=0)
ref_mmm['prc75'] = np.percentile(ref_indiv['racs'], 75, axis=0)
ref_mmm['std'] = np.nanstd(ref_indiv['racs'], axis=0)

# RRTMG full RCP8.5 response
full = {}
[full['ra'], full['racs'], full['diag']] = pickle.load(open('%s/rad.rcp.%s.pickle' % (datadir, model), 'rb'))

full_indiv = {}
full_indiv['racs'] = np.empty([len(models), len(full['racs'])])
for im in range(len(models)):
    [ra, racs, diag] = pickle.load(open('%s/rad.rcp.%s.pickle' % (datadir, models[im]), 'rb'))
    full_indiv['racs'][im,:] = racs
full_mmm = {}
full_mmm['prc25'] = {}
full_mmm['prc75'] = {}
full_mmm['std'] = {}
full_mmm['prc25'] = np.percentile(full_indiv['racs'], 25, axis=0)
full_mmm['prc75'] = np.percentile(full_indiv['racs'], 75, axis=0)
full_mmm['std'] = np.nanstd(full_indiv['racs'], axis=0)

# RRTMG fixed T and q RCP8.5 response
fixqt = {}
[fixqt['ra'], fixqt['racs'], fixqt['diag']] = pickle.load(open('%s/rad.rcp_fixq_fixT.%s.pickle' % (datadir, model), 'rb'))

fixqt_indiv = {}
fixqt_indiv['racs'] = np.empty([len(models), len(fixqt['racs'])])
for im in range(len(models)):
    [ra, racs, diag] = pickle.load(open('%s/rad.rcp_fixq_fixT.%s.pickle' % (datadir, models[im]), 'rb'))
    fixqt_indiv['racs'][im,:] = racs
fixqt_mmm = {}
fixqt_mmm['prc25'] = {}
fixqt_mmm['prc75'] = {}
fixqt_mmm['std'] = {}
fixqt_mmm['prc25'] = np.percentile(fixqt_indiv['racs'], 25, axis=0)
fixqt_mmm['prc75'] = np.percentile(fixqt_indiv['racs'], 75, axis=0)
fixqt_mmm['std'] = np.nanstd(fixqt_indiv['racs'], axis=0)

# RRTMG fixed q and CO2 RCP8.5 response
fixqco2 = {}
[fixqco2['ra'], fixqco2['racs'], fixqco2['diag']] = pickle.load(open('%s/rad.rcp_fixq_fixco2.%s.pickle' % (datadir, model), 'rb'))

fixqco2_indiv = {}
fixqco2_indiv['racs'] = np.empty([len(models), len(fixqco2['racs'])])
for im in range(len(models)):
    [ra, racs, diag] = pickle.load(open('%s/rad.rcp_fixq_fixco2.%s.pickle' % (datadir, models[im]), 'rb'))
    fixqco2_indiv['racs'][im,:] = racs
fixqco2_mmm = {}
fixqco2_mmm['prc25'] = {}
fixqco2_mmm['prc75'] = {}
fixqco2_mmm['std'] = {}
fixqco2_mmm['prc25'] = np.percentile(fixqco2_indiv['racs'], 25, axis=0)
fixqco2_mmm['prc75'] = np.percentile(fixqco2_indiv['racs'], 75, axis=0)
fixqco2_mmm['std'] = np.nanstd(fixqco2_indiv['racs'], axis=0)

# RRTMG fixed T and CO2 RCP8.5 response
fixtco2 = {}
[fixtco2['ra'], fixtco2['racs'], fixtco2['diag']] = pickle.load(open('%s/rad.rcp_fixT_fixco2.%s.pickle' % (datadir, model), 'rb'))

fixtco2_indiv = {}
fixtco2_indiv['racs'] = np.empty([len(models), len(fixtco2['racs'])])
for im in range(len(models)):
    [ra, racs, diag] = pickle.load(open('%s/rad.rcp_fixT_fixco2.%s.pickle' % (datadir, models[im]), 'rb'))
    fixtco2_indiv['racs'][im,:] = racs
fixtco2_mmm = {}
fixtco2_mmm['prc25'] = {}
fixtco2_mmm['prc75'] = {}
fixtco2_mmm['std'] = {}
fixtco2_mmm['prc25'] = np.percentile(fixtco2_indiv['racs'], 25, axis=0)
fixtco2_mmm['prc75'] = np.percentile(fixtco2_indiv['racs'], 75, axis=0)
fixtco2_mmm['std'] = np.nanstd(fixtco2_indiv['racs'], axis=0)

# RRTMG fixed relative humidity and CO2 RCP8.5 response
fixrhco2 = {}
[fixrhco2['ra'], fixrhco2['racs'], fixrhco2['diag']] = pickle.load(open('%s/rad.rcp_fixrh_fixco2.%s.pickle' % (datadir, model), 'rb'))

fixrhco2_indiv = {}
fixrhco2_indiv['racs'] = np.empty([len(models), len(fixrhco2['racs'])])
for im in range(len(models)):
    [ra, racs, diag] = pickle.load(open('%s/rad.rcp_fixrh_fixco2.%s.pickle' % (datadir, models[im]), 'rb'))
    fixrhco2_indiv['racs'][im,:] = racs
fixrhco2_mmm = {}
fixrhco2_mmm['prc25'] = {}
fixrhco2_mmm['prc75'] = {}
fixrhco2_mmm['std'] = {}
fixrhco2_mmm['prc25'] = np.percentile(fixrhco2_indiv['racs'], 25, axis=0)
fixrhco2_mmm['prc75'] = np.percentile(fixrhco2_indiv['racs'], 75, axis=0)
fixrhco2_mmm['std'] = np.nanstd(fixrhco2_indiv['racs'], axis=0)

# RRTMG fixed T and CO2, devRH RCP8.5 response
fixtco2devrh = {}
[fixtco2devrh['ra'], fixtco2devrh['racs'], fixtco2devrh['diag']] = pickle.load(open('%s/rad.rcp_fixT_fixco2_devrh.%s.pickle' % (datadir, model), 'rb'))

fixtco2devrh_indiv = {}
fixtco2devrh_indiv['racs'] = np.empty([len(models), len(fixtco2devrh['racs'])])
for im in range(len(models)):
    [ra, racs, diag] = pickle.load(open('%s/rad.rcp_fixT_fixco2_devrh.%s.pickle' % (datadir, models[im]), 'rb'))
    fixtco2devrh_indiv['racs'][im,:] = racs
fixtco2devrh_mmm = {}
fixtco2devrh_mmm['prc25'] = {}
fixtco2devrh_mmm['prc75'] = {}
fixtco2devrh_mmm['std'] = {}
fixtco2devrh_mmm['prc25'] = np.percentile(fixtco2devrh_indiv['racs'], 25, axis=0)
fixtco2devrh_mmm['prc75'] = np.percentile(fixtco2devrh_indiv['racs'], 75, axis=0)
fixtco2devrh_mmm['std'] = np.nanstd(fixtco2devrh_indiv['racs'], axis=0)

# # RRTMG fixed CO2 RCP8.5 response
# fixco2 = {}
# [fixco2['ra'], fixco2['racs'], fixco2['diag']] = pickle.load(open('%s/rad.rcp_fixco2.%s.pickle' % (datadir, model), 'rb'))

# # RRTMG fixed specific humidity RCP8.5 response
# fixq = {}
# [fixq['ra'], fixq['racs'], fixq['diag']] = pickle.load(open('%s/rad.rcp_fixq.%s.pickle' % (datadir, model), 'rb'))

# # RRTMG fixed relative humidity RCP8.5 response
# fixrh = {}
# [fixrh['ra'], fixrh['racs'], fixrh['diag']] = pickle.load(open('%s/rad.rcp_fixrh.%s.pickle' % (datadir, model), 'rb'))

##################################################
# MISC data
##################################################
plev = 1e-2*grid['lev']
modelstr = get_modelstr(model)

##################################################
# RCP
##################################################
time = np.arange(2006,2301)

for varname in gcm_clim:
    gcm_clim[varname] = np.nanmean(gcm_clim[varname][-30:,...], axis=0)

ref['ra'] = np.nanmean(ref['ra'][-30:])
ref['racs'] = np.nanmean(ref['racs'][-30:])

for varname in gcm_clim:
    gcm_clim_mmm['prc25'][varname] = np.nanmean(gcm_clim_mmm['prc25'][varname][-30:,...], axis=0)
    gcm_clim_mmm['prc75'][varname] = np.nanmean(gcm_clim_mmm['prc75'][varname][-30:,...], axis=0)
    gcm_clim_mmm['std'][varname] = np.nanmean(gcm_clim_mmm['std'][varname][-30:,...], axis=0)

for stat in ref_mmm:
    ref_mmm[stat] = np.nanmean(ref_mmm[stat][-30:])

##################################################
# FIX q and T
##################################################
# compare change in radiative cooling
plotname='%s/dracs_fixq_fixT.pdf' % (plotdir)
fig, ax = plt.subplots()
ax.axhline(0,color='k', linewidth=0.5)
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.plot(time, full['racs']-ref['racs'], color='black', label='RRTMG')
ax.plot(time, fixqt['racs']-ref['racs'], '-', color='tab:green', label='RRTMG (fixq+fixT)')
ax.set_xlim([time[0], time[-1]])
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('Time (yr)')
ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
fig.set_size_inches(4,3.3)
plt.legend()
plt.tight_layout()
plt.savefig(plotname, format='pdf', dpi=300)
plt.close()

##################################################
# FIX q and CO2
##################################################
# compare change in radiative cooling
plotname='%s/dracs_fixq_fixco2.pdf' % (plotdir)
fig, ax = plt.subplots()
ax.axhline(0,color='k', linewidth=0.5)
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.plot(time, full['racs']-ref['racs'], color='black', label='RRTMG')
ax.plot(time, fixqco2['racs']-ref['racs'], '-', color='tab:orange', label='RRTMG (fixq+fixCO2)')
ax.set_xlim([time[0], time[-1]])
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('Time (yr)')
ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
fig.set_size_inches(4,3.3)
plt.legend()
plt.tight_layout()
plt.savefig(plotname, format='pdf', dpi=300)
plt.close()

##################################################
# FIX T and CO2
##################################################
# compare change in radiative cooling
plotname='%s/dracs_fixT_fixco2.pdf' % (plotdir)
fig, ax = plt.subplots()
ax.axhline(0,color='k', linewidth=0.5)
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.plot(time, full['racs']-ref['racs'], color='black', label='RRTMG')
ax.plot(time, fixtco2['racs']-ref['racs'], '-', color='tab:blue', label='RRTMG (fixT+fixCO2)')
ax.set_xlim([time[0], time[-1]])
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('Time (yr)')
ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
fig.set_size_inches(4,3.3)
plt.legend()
plt.tight_layout()
plt.savefig(plotname, format='pdf', dpi=300)
plt.close()

##################################################
# FIX RH and CO2
##################################################
# compare change in radiative cooling
plotname='%s/dracs_fixrh_fixco2.pdf' % (plotdir)
fig, ax = plt.subplots()
ax.axhline(0,color='k', linewidth=0.5)
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.plot(time, full['racs']-ref['racs'], color='black', label='RRTMG')
ax.plot(time, fixrhco2['racs']-ref['racs'], '--', color='tab:blue', label='RRTMG (fixRH+fixCO2)')
ax.set_xlim([time[0], time[-1]])
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('Time (yr)')
ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
fig.set_size_inches(4,3.3)
plt.legend()
plt.tight_layout()
plt.savefig(plotname, format='pdf', dpi=300)
plt.close()

##################################################
# DECOMP CO2, T, q (no legend)
##################################################
# compare change in radiative cooling
plotname='%s/dracs_decomp_noleg.pdf' % (plotdir)
fig, ax = plt.subplots()
ax.axhline(0,color='k', linewidth=0.5)
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
lres,=ax.plot(time, full['racs']-fixqt['racs']-fixqco2['racs']-fixtco2['racs']+2*ref['racs'], '-.', color='k', label='Residual')
lfull,=ax.plot(time, full['racs']-ref['racs'], color='black', label='All')
lco2,=ax.plot(time, fixqt['racs']-ref['racs'], '-', color='tab:green', label='CO$_2$')
lt,=ax.plot(time, fixqco2['racs']-ref['racs'], '-', color='tab:orange', label='T')
lq,=ax.plot(time, fixtco2['racs']-ref['racs'], '-', color='tab:blue', label='Q')
ax.set_xlim([time[0], time[-1]])
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('Time (yr)')
ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
if model == 'mmm':
    ax.set_title('RRTMG')
else:
    ax.set_title('RRTMG(%s)' % (modelstr))
fig.set_size_inches(4,3.3)
plt.tight_layout()
plt.savefig(plotname, format='pdf', dpi=300)
plt.close()

##################################################
# DECOMP CO2, T, q (inside legend)
##################################################
# compare change in radiative cooling
plotname='%s/dracs_decomp.pdf' % (plotdir)
fig, ax = plt.subplots()
ax.axhline(0,color='k', linewidth=0.5)
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)

rrtmgstd = full_mmm['std'] - ref_mmm['std']
ax.fill_between(time, full['racs']-ref['racs']-rrtmgstd, full['racs']-ref['racs']+rrtmgstd, alpha=0.2, color='k', edgecolor=None)
rrtmgstd = fixqt_mmm['std'] - ref_mmm['std']
ax.fill_between(time, fixqt['racs']-ref['racs']-rrtmgstd, fixqt['racs']-ref['racs']+rrtmgstd, alpha=0.2, color='tab:green', edgecolor=None)
rrtmgstd = fixqco2_mmm['std'] - ref_mmm['std']
ax.fill_between(time, fixqco2['racs']-ref['racs']-rrtmgstd, fixqco2['racs']-ref['racs']+rrtmgstd, alpha=0.2, color='tab:orange', edgecolor=None)
rrtmgstd = fixtco2_mmm['std'] - ref_mmm['std']
ax.fill_between(time, fixtco2['racs']-ref['racs']-rrtmgstd, fixtco2['racs']-ref['racs']+rrtmgstd, alpha=0.2, color='tab:blue', edgecolor=None)
res=full['racs']-fixqt['racs']-fixqco2['racs']-fixtco2['racs']+2*ref['racs']
rrtmgstd=full_mmm['std']-fixqt_mmm['std']-fixqco2_mmm['std']-fixtco2_mmm['std']+2*ref_mmm['std']
ax.fill_between(time, res-rrtmgstd, res+rrtmgstd, alpha=0.2, color='k', edgecolor=None)

lres,=ax.plot(time, res, '-.', color='k', label='Residual')
lfull,=ax.plot(time, full['racs']-ref['racs'], color='black', label='All')
lco2,=ax.plot(time, fixqt['racs']-ref['racs'], '-', color='tab:green', label=r'$\Delta CO_2$')
lt,=ax.plot(time, fixqco2['racs']-ref['racs'], '-', color='tab:orange', label=r'$\Delta T$')
lq,=ax.plot(time, fixtco2['racs']-ref['racs'], '-', color='tab:blue', label=r'$\Delta q$')
ax.set_xlim([time[0], time[-1]])
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('Time (yr)')
ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
if model == 'mmm':
    ax.set_title('RRTMG')
else:
    ax.set_title('RRTMG(%s)' % (modelstr))
fig.set_size_inches(4,3.3)
plt.tight_layout()
leg=plt.legend(handles=[lfull,lco2,lt,lq,lres], loc='lower left', prop={'size':8})
leg.get_frame().set_linewidth(0)
leg.get_frame().set_facecolor('none')
plt.savefig(plotname, format='pdf', dpi=300)
plt.close()

##################################################
# DECOMP CO2, T, q (w/ legend)
##################################################
# compare change in radiative cooling
plotname='%s/dracs_decomp_leg.pdf' % (plotdir)
fig, ax = plt.subplots()
ax.axhline(0,color='k', linewidth=0.5)
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
lres,=ax.plot(time, full['racs']-fixqt['racs']-fixqco2['racs']-fixtco2['racs']+2*ref['racs'], '-.', color='k', label='Residual')
lfull,=ax.plot(time, full['racs']-ref['racs'], color='black', label='All')
lco2,=ax.plot(time, fixqt['racs']-ref['racs'], '-', color='tab:green', label='CO$_2$')
lt,=ax.plot(time, fixqco2['racs']-ref['racs'], '-', color='tab:orange', label='T')
lq,=ax.plot(time, fixtco2['racs']-ref['racs'], '-', color='tab:blue', label='Q')
ax.set_xlim([time[0], time[-1]])
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('Time (yr)')
ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
ax.set_title(modelstr)
fig.set_size_inches(6,5)
plt.tight_layout()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.15,
                     box.width, box.height * 0.85])
ax.legend(handles=[lfull,lco2,lt,lq,lres], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
plt.savefig(plotname, format='pdf', dpi=300)
plt.close()

##################################################
# DECOMP CO2, T, RH (no legend)
##################################################
# compare change in radiative cooling
plotname='%s/dracs_decomp_rh_noleg.pdf' % (plotdir)
fig, ax = plt.subplots()
ax.axhline(0,color='k', linewidth=0.5)
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
lres,=ax.plot(time, full['racs']-fixqt['racs']-fixrhco2['racs']-fixtco2devrh['racs']+2*ref['racs'], '-.', color='k', label='Residual')
lfull,=ax.plot(time, full['racs']-ref['racs'], color='black', label='All')
lco2,=ax.plot(time, fixqt['racs']-ref['racs'], '-', color='tab:green', label='CO$_2$')
lcc,=ax.plot(time, fixrhco2['racs']-ref['racs'], '-', color='tab:orange', label=r'T+CC')
ldcc,=ax.plot(time, fixtco2devrh['racs']-ref['racs'], '-', color='tab:blue', label=r'RH')
ax.set_xlim([time[0], time[-1]])
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('Time (yr)')
ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
if model == 'mmm':
    ax.set_title('RRTMG')
else:
    ax.set_title('RRTMG(%s)' % (modelstr))
fig.set_size_inches(4,3.3)
plt.tight_layout()
plt.savefig(plotname, format='pdf', dpi=300)
plt.close()

##################################################
# DECOMP CO2, T, RH (inside legend)
##################################################
# compare change in radiative cooling
plotname='%s/dracs_decomp_rh.pdf' % (plotdir)
fig, ax = plt.subplots()
ax.axhline(0,color='k', linewidth=0.5)
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)

rrtmgstd = full_mmm['std'] - ref_mmm['std']
ax.fill_between(time, full['racs']-ref['racs']-rrtmgstd, full['racs']-ref['racs']+rrtmgstd, alpha=0.2, color='k', edgecolor=None)
rrtmgstd = fixqt_mmm['std'] - ref_mmm['std']
ax.fill_between(time, fixqt['racs']-ref['racs']-rrtmgstd, fixqt['racs']-ref['racs']+rrtmgstd, alpha=0.2, color='tab:green', edgecolor=None)
rrtmgstd = fixrhco2_mmm['std'] - ref_mmm['std']
ax.fill_between(time, fixrhco2['racs']-ref['racs']-rrtmgstd, fixrhco2['racs']-ref['racs']+rrtmgstd, alpha=0.2, color='tab:orange', edgecolor=None)
rrtmgstd = fixtco2devrh_mmm['std'] - ref_mmm['std']
ax.fill_between(time, fixtco2devrh['racs']-ref['racs']-rrtmgstd, fixtco2devrh['racs']-ref['racs']+rrtmgstd, alpha=0.2, color='tab:blue', edgecolor=None)
res=full['racs']-fixqt['racs']-fixrhco2['racs']-fixtco2devrh['racs']+2*ref['racs']
rrtmgstd=full_mmm['std']-fixqt_mmm['std']-fixrhco2_mmm['std']-fixtco2devrh_mmm['std']+2*ref_mmm['std']
ax.fill_between(time, res-rrtmgstd, res+rrtmgstd, alpha=0.2, color='k', edgecolor=None)

lres,=ax.plot(time, res, '-.', color='k', label='Residual')
lfull,=ax.plot(time, full['racs']-ref['racs'], color='black', label='All')
lco2,=ax.plot(time, fixqt['racs']-ref['racs'], '-', color='tab:green', label=r'$\Delta CO_2$')
lcc,=ax.plot(time, fixrhco2['racs']-ref['racs'], '-', color='tab:orange', label=r'$\Delta T,\Delta q|_\mathrm{fix\,RH}$')
ldcc,=ax.plot(time, fixtco2devrh['racs']-ref['racs'], '-', color='tab:blue', label=r'$\Delta RH q^*|_\mathrm{clim}$')
ax.set_xlim([time[0], time[-1]])
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('Time (yr)')
ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
if model == 'mmm':
    ax.set_title('RRTMG')
else:
    ax.set_title('RRTMG(%s)' % (modelstr))
fig.set_size_inches(4,3.3)
plt.tight_layout()
plt.tight_layout()
leg=plt.legend(handles=[lfull,lco2,lcc,ldcc,lres], prop={'size':8})
leg.get_frame().set_linewidth(0)
leg.get_frame().set_facecolor('none')
plt.savefig(plotname, format='pdf', dpi=300)
plt.close()

##################################################
# DECOMP CO2, T, RH (w/ legend)
##################################################
# compare change in radiative cooling
plotname='%s/dracs_decomp_rh_leg.pdf' % (plotdir)
fig, ax = plt.subplots()
ax.axhline(0,color='k', linewidth=0.5)
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
lres,=ax.plot(time, full['racs']-fixqt['racs']-fixrhco2['racs']-fixtco2devrh['racs']+2*ref['racs'], '-.', color='k', label='Residual')
lfull,=ax.plot(time, full['racs']-ref['racs'], color='black', label='All')
lco2,=ax.plot(time, fixqt['racs']-ref['racs'], '-', color='tab:green', label='CO$_2$')
lcc,=ax.plot(time, fixrhco2['racs']-ref['racs'], '-', color='tab:orange', label=r'T+CC')
ldcc,=ax.plot(time, fixtco2devrh['racs']-ref['racs'], '-', color='tab:blue', label=r'RH')
ax.set_xlim([time[0], time[-1]])
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('Time (yr)')
ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
ax.set_title(modelstr)
fig.set_size_inches(5.5,4.5)
plt.tight_layout()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.15,
                     box.width, box.height * 0.85])
ax.legend(handles=[lfull,lco2,lcc,ldcc,lres], loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)
plt.savefig(plotname, format='pdf', dpi=300)
plt.close()

###################################################
## FIX q and T (compared to CO2)
###################################################
#vmin = -111; vmax=-108;
#ra_f30 = np.mean(fixqt['racs'][:30])
#ra_l30 = np.mean(fixqt['racs'][-30:])
#co2_f30 = np.mean(co2[:30])
#co2_l30 = np.mean(co2[-30:])
#maxis = (co2_l30 - co2_f30) / (ra_l30 - ra_f30)
#vmax_alg = co2_f30 + maxis*(vmax-ra_f30)
#vmin_alg = co2_l30 + maxis*(vmin-ra_l30)
## compare absolute change in radiative cooling
#plotname='%s/co2_racs_fixq_fixT.pdf' % (plotdir)
#fig, ax = plt.subplots()
#ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
#ax.plot(time, fixqt['racs'], '-', color='tab:green', label='RRTMG (fixq+T)')
#ax.set_xlim([time[0], time[-1]])
#ax.xaxis.set_minor_locator(AutoMinorLocator())
#ax.yaxis.set_minor_locator(AutoMinorLocator())
#ax.set_xlabel('Time (yr)')
#ax.set_ylabel(r'$R_{a,\,cs}$ (W m$^{-2}$)')
#ax.set_ylim(vmin, vmax)
#sax = ax.twinx()
#sax.plot(time, co2, '-', color='tab:purple', label='p_{CO_2}')
#sax.set_ylim(vmin_alg, vmax_alg)
#fig.set_size_inches(4,3.3)
#plt.legend()
#plt.tight_layout()
#plt.savefig(plotname, format='pdf', dpi=300)
#plt.close()

###################################################
## FIX Q
###################################################
## compare absolute change in radiative cooling
#plotname='%s/racs_fixq.pdf' % (plotdir)
#fig, ax = plt.subplots()
#ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
#ax.plot(time, gcm_forc['racs'], color='tab:red', label=modelstr)
#ax.plot(time, full['racs'], color='black', label='RRTMG')
#ax.plot(time, fixq['racs'], color='tab:orange', label='RRTMG (fixq)')
#ax.set_xlim([time[0], time[-1]])
#ax.xaxis.set_minor_locator(AutoMinorLocator())
#ax.yaxis.set_minor_locator(AutoMinorLocator())
#ax.set_xlabel('Time (yr)')
#ax.set_ylabel(r'$R_{a,\,cs}$ (W m$^{-2}$)')
#fig.set_size_inches(4,3.3)
#plt.legend()
#plt.tight_layout()
#plt.savefig(plotname, format='pdf', dpi=300)
#plt.close()

###################################################
## FIX RH
###################################################
## compare absolute change in radiative cooling
#plotname='%s/racs_fixrh.pdf' % (plotdir)
#fig, ax = plt.subplots()
#ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
#ax.plot(time, gcm_forc['racs'], color='tab:red', label=modelstr)
#ax.plot(time, full['racs'], color='black', label='RRTMG')
#ax.plot(time, fixrh['racs'], '--', color='tab:blue', label='RRTMG (fixrh)')
#ax.set_xlim([time[0], time[-1]])
#ax.xaxis.set_minor_locator(AutoMinorLocator())
#ax.yaxis.set_minor_locator(AutoMinorLocator())
#ax.set_xlabel('Time (yr)')
#ax.set_ylabel(r'$R_{a,\,cs}$ (W m$^{-2}$)')
#fig.set_size_inches(4,3.3)
#plt.legend()
#plt.tight_layout()
#plt.savefig(plotname, format='pdf', dpi=300)
#plt.close()

###################################################
## FIX RH and CO2
###################################################
## compare absolute change in radiative cooling
#plotname='%s/racs_fixrh_fixco2.pdf' % (plotdir)
#fig, ax = plt.subplots()
#ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
#ax.plot(time, gcm_forc['racs'], color='tab:red', label=modelstr)
#ax.plot(time, full['racs'], color='black', label='RRTMG')
#ax.plot(time, fixrhco2['racs'], '--', color='tab:blue', label='RRTMG (fixrh+co2)')
#ax.set_xlim([time[0], time[-1]])
#ax.xaxis.set_minor_locator(AutoMinorLocator())
#ax.yaxis.set_minor_locator(AutoMinorLocator())
#ax.set_xlabel('Time (yr)')
#ax.set_ylabel(r'$R_{a,\,cs}$ (W m$^{-2}$)')
#fig.set_size_inches(4,3.3)
#plt.legend()
#plt.tight_layout()
#plt.savefig(plotname, format='pdf', dpi=300)
#plt.close()

###################################################
## FIX CO2
###################################################
## compare absolute change in radiative cooling
#plotname='%s/racs_fixco2.pdf' % (plotdir)
#fig, ax = plt.subplots()
#ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
#ax.plot(time, gcm_forc['racs'], color='tab:red', label=modelstr)
#ax.plot(time, full['racs'], color='black', label='RRTMG')
#ax.plot(time, fixco2['racs'], color='tab:blue', label='RRTMG (fixco2)')
#ax.set_xlim([time[0], time[-1]])
#ax.xaxis.set_minor_locator(AutoMinorLocator())
#ax.yaxis.set_minor_locator(AutoMinorLocator())
#ax.set_xlabel('Time (yr)')
#ax.set_ylabel(r'$R_{a,\,cs}$ (W m$^{-2}$)')
#fig.set_size_inches(4,3.3)
#plt.legend()
#plt.tight_layout()
#plt.savefig(plotname, format='pdf', dpi=300)
#plt.close()
