import os
import sys
import pickle
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from tools import get_model6idx, get_model6str

# allmodels = ['ACCESS-CM2', 'ACCESS-ESM1-5', 'CanESM5', 'CESM2-WACCM', 'IPSL-CM6A-LR', 'MRI-ESM2-0','MIROC-ES2L','GISS-E2-1-G','GISS-E2-1-H','UKESM1-0-LL']
# models = ['mmm']
# from scipy.stats import t
# alpha=0.1 # significance
# df=len(allmodels)
# tscore=t.ppf(1-alpha/2,df)

models = ['ACCESS-CM2', 'ACCESS-ESM1-5', 'CanESM5', 'CESM2-WACCM', 'IPSL-CM6A-LR', 'MRI-ESM2-0','MIROC-ES2L','GISS-E2-1-G','GISS-E2-1-H','UKESM1-0-LL']

for model in models:

    # specify directory for output data and plots
    datadir = '%s/data' % ( os.getcwd() ) # output data directory
    plotdir = '%s/plot/%s' % ( os.getcwd(), model ) # plot directory
    if not os.path.exists(plotdir):
        os.mkdir(plotdir)

    ##################################################
    # LOAD GCM DATA
    ##################################################
    ds = xr.open_dataset('/project2/tas1/miyawaki/projects/003/echam/ghg/greenhouse_hist+ssp585.nc')
    co2 = 1e-6*ds.CO2.sel(time=slice('2005-12-31', '2299-12-31')).data

    [mmm_clim, grid, indiv_clim] = pickle.load(open('../climlab/input_data/clima.cmip6.pickle', 'rb'))
    [mmm_forc, grid, indiv_forc] = pickle.load(open('../climlab/input_data/forcing.cmip6.pickle', 'rb'))

    if model == 'mmm':
        gcm_clim = mmm_clim
        gcm_forc = mmm_forc

        indiv_resp={}
        gcm_resp = {}
        gcm_resp_mmm = {}
        gcm_resp_mmm['std'] = {}

        gcm_clim_mmm = {}
        gcm_clim_mmm['std'] = {}
        gcm_forc_mmm = {}
        gcm_forc_mmm['std'] = {}
        for varname in indiv_clim:
            gcm_clim_mmm['std'][varname] = np.nanstd(indiv_clim[varname], axis=0)
            gcm_forc_mmm['std'][varname] = np.nanstd(indiv_forc[varname], axis=0)

            gcm_resp[varname]=gcm_forc[varname]-np.nanmean(gcm_clim[varname][-30:,...],axis=0,keepdims=True)
            indiv_resp[varname]=indiv_forc[varname]-np.nanmean(indiv_clim[varname][:,-30:,...],axis=1,keepdims=True)
            gcm_resp_mmm['std'][varname] = np.nanstd(indiv_resp[varname], axis=0)
    else:
        imod = get_model6idx(model)
        gcm_clim = {}
        gcm_forc = {}
        gcm_resp = {}
        indiv_resp={}
        for varname in indiv_clim:
            gcm_clim[varname] = indiv_clim[varname][imod,...]
            gcm_forc[varname] = indiv_forc[varname][imod,...]
            gcm_resp[varname]=gcm_forc[varname]-np.nanmean(gcm_clim[varname][-30:,...],axis=0,keepdims=True)
            indiv_resp[varname]=indiv_forc[varname]-np.nanmean(indiv_clim[varname][:,-30:,...],axis=1,keepdims=True)

    ##################################################
    # LOAD RRTMG DATA
    ##################################################
    # RRTMG climatology
    ref = {}
    [ref['ra'], ref['racs'], ref['diag']] = pickle.load(open('%s/rad.hist.%s.cmip6.pickle' % (datadir, model), 'rb'))

    # RRTMG full RCP8.5 response
    full = {}
    [full['ra'], full['racs'], full['diag']] = pickle.load(open('%s/rad.hist+ssp585.%s.cmip6.pickle' % (datadir, model), 'rb'))

    print(ref['ra'].shape)
    
    if model == 'mmm':
        ref_indiv = {}
        ref_indiv['racs'] = np.empty([len(allmodels), len(ref['racs'])])
        for im in range(len(allmodels)):
            [ra, racs, diag] = pickle.load(open('%s/rad.hist.%s.cmip6.pickle' % (datadir, allmodels[im]), 'rb'))
            ref_indiv['racs'][im,:] = racs
        ref_mmm = {}
        ref_mmm['std'] = {}
        ref_mmm['std'] = np.nanstd(ref_indiv['racs'], axis=0)

        full_indiv = {}
        full_indiv['racs'] = np.empty([len(allmodels), len(full['racs'])])
        for im in range(len(allmodels)):
            [ra, racs, diag] = pickle.load(open('%s/rad.hist+ssp585.%s.cmip6.pickle' % (datadir, allmodels[im]), 'rb'))
            full_indiv['racs'][im,:] = racs
        full_mmm = {}
        full_mmm['std'] = {}
        full_mmm['std'] = np.nanstd(full_indiv['racs'], axis=0)

    ##################################################
    # MISC data
    ##################################################
    plev = 1e-2*grid['lev']
    modelstr = get_model6str(model)

    ##################################################
    # HISTORICAL
    ##################################################
    time0 = np.arange(1860,2015+1)

    ##################################################
    # TOTAL RA
    ##################################################
    # compare absolute change in radiative cooling
    plotname='%s/ra_hist.cmip6.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time0, gcm_clim['ra'], color='tab:gray', label=modelstr)
    ax.plot(time0, ref['ra'], color='black', label='RRTMG')
    ax.set_xlim([time0[0], time0[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$R_a$ (W m$^{-2}$)')
    fig.set_size_inches(4,3.3)
    if model=='mmm':
        plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # CLEAR SKY RA
    ##################################################
    # compare absolute change in radiative cooling
    plotname='%s/racs_hist.cmip6.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time0, gcm_clim['racs'], color='tab:red', label=modelstr)
    ax.plot(time0, ref['racs'], color='black', label='RRTMG')
    ax.set_xlim([time0[0], time0[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$R_{a,\,cs}$ (W m$^{-2}$)')
    fig.set_size_inches(4,3.3)
    if model=='mmm':
        plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # RCP
    ##################################################
    # time = np.arange(2006,2301)
    time = np.arange(1860,2300+1)

    for varname in gcm_clim:
        gcm_clim[varname] = np.nanmean(gcm_clim[varname][-30:,...], axis=0)

    ref['ra'] = np.nanmean(ref['ra'][-30:])
    ref['racs'] = np.nanmean(ref['racs'][-30:])

    if model == 'mmm':
        for varname in gcm_clim:
            gcm_clim_mmm['std'][varname] = np.nanmean(gcm_clim_mmm['std'][varname][-30:,...], axis=0)

        for stat in ref_mmm:
            ref_mmm[stat] = np.nanmean(ref_mmm[stat][-30:])

    ##################################################
    # ALL COMP DRA
    ##################################################
    # compare relative change in radiative cooling
    plotname='%s/dra_comp.cmip6.pdf' % (plotdir)

    sw_forc = gcm_forc['rsdt'] - gcm_forc['rsut'] + gcm_forc['rsus'] - gcm_forc['rsds']
    lw_forc = - gcm_forc['rlut'] + gcm_forc['rlus'] - gcm_forc['rlds']
    lwcs_forc = - gcm_forc['rlutcs'] + gcm_forc['rlus'] - gcm_forc['rldscs']
    lwcld_forc = lw_forc - lwcs_forc
    sw_clim = gcm_clim['rsdt'] - gcm_clim['rsut'] + gcm_clim['rsus'] - gcm_clim['rsds']
    lw_clim = - gcm_clim['rlut'] + gcm_clim['rlus'] - gcm_clim['rlds']
    lwcs_clim = - gcm_clim['rlutcs'] + gcm_clim['rlus'] - gcm_clim['rldscs']
    lwcld_clim = lw_clim - lwcs_clim

    dsw = gcm_resp['rsdt'] - gcm_resp['rsut'] + gcm_resp['rsus'] - gcm_resp['rsds']
    dlw = - gcm_resp['rlut'] + gcm_resp['rlus'] - gcm_resp['rlds']
    dlwcs = - gcm_resp['rlutcs'] + gcm_resp['rlus'] - gcm_resp['rldscs']
    dlwcld = dlw - dlwcs

    # dsw = sw_forc - sw_clim
    # dlw = lw_forc - lw_clim
    # dlwcs = lwcs_forc - lwcs_clim
    # dlwcld = lwcld_forc - lwcld_clim

    if model == 'mmm':
        sw_forcstd = gcm_forc_mmm['std']['rsdt'] - gcm_forc_mmm['std']['rsut'] + gcm_forc_mmm['std']['rsus'] - gcm_forc_mmm['std']['rsds']
        lw_forcstd = - gcm_forc_mmm['std']['rlut'] + gcm_forc_mmm['std']['rlus'] - gcm_forc_mmm['std']['rlds']
        lwcs_forcstd = - gcm_forc_mmm['std']['rlutcs'] + gcm_forc_mmm['std']['rlus'] - gcm_forc_mmm['std']['rldscs']
        lwcld_forcstd = lw_forcstd - lwcs_forcstd
        sw_climstd = gcm_clim_mmm['std']['rsdt'] - gcm_clim_mmm['std']['rsut'] + gcm_clim_mmm['std']['rsus'] - gcm_clim_mmm['std']['rsds']
        lw_climstd = - gcm_clim_mmm['std']['rlut'] + gcm_clim_mmm['std']['rlus'] - gcm_clim_mmm['std']['rlds']
        lwcs_climstd = - gcm_clim_mmm['std']['rlutcs'] + gcm_clim_mmm['std']['rlus'] - gcm_clim_mmm['std']['rldscs']
        lwcld_climstd = lw_climstd - lwcs_climstd

        dswstd = np.sqrt(1/4*(gcm_resp_mmm['std']['rsdt']**2 + gcm_resp_mmm['std']['rsut']**2 + gcm_resp_mmm['std']['rsus']**2 + gcm_resp_mmm['std']['rsds']**2))
        dlwstd = np.sqrt(1/3*(gcm_resp_mmm['std']['rlut']**2 + gcm_resp_mmm['std']['rlus']**2 + gcm_resp_mmm['std']['rlds']**2))
        dlwcsstd = np.sqrt(1/3*(gcm_resp_mmm['std']['rlutcs']**2 + gcm_resp_mmm['std']['rlus']**2 + gcm_resp_mmm['std']['rldscs']**2))
        dlwcldstd = np.sqrt(1/2*(dlwstd**2 + dlwcsstd**2))

        # dswstd = sw_forcstd - sw_climstd
        # dlwstd = lw_forcstd - lw_climstd
        # dlwcsstd = lwcs_forcstd - lwcs_climstd
        # dlwcldstd = lwcld_forcstd - lwcld_climstd

    fig, ax = plt.subplots()
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if model == 'mmm':
        sqn=2*np.sqrt(len(allmodels))/tscore
        # rastd = np.sqrt(1/2*(gcm_forc_mmm['std']['ra']**2+gcm_clim_mmm['std']['ra']**2))
        rastd = gcm_resp_mmm['std']['ra']
        ax.fill_between(time, gcm_forc['ra']-gcm_clim['ra'] - rastd/sqn, gcm_forc['ra']-gcm_clim['ra'] + rastd/sqn, alpha=0.2, color='tab:gray', edgecolor=None)
        ax.fill_between(time, dsw-dswstd/sqn, dsw+dswstd/sqn, alpha=0.2, color='cyan', edgecolor=None)
        ax.fill_between(time, dlwcs-dlwcsstd/sqn, dlwcs+dlwcsstd/sqn, alpha=0.2, color='tab:red', edgecolor=None)
        ax.fill_between(time, dlwcld-dlwcldstd/sqn, dlwcld+dlwcldstd/sqn, alpha=0.2, color='tab:purple', edgecolor=None)
        # ax.fill_between(time, gcm_forc['ra']-gcm_clim['ra'] - rastd/sqn, gcm_forc['ra']-gcm_clim['ra'] + rastd/sqn, alpha=0.2, color='k', edgecolor=None)

        # rrtmgstd = np.sqrt(1/2*(full_mmm['std']**2 + ref_mmm['std']**2))
        # ax.fill_between(time, full['racs']-ref['racs']-rrtmgstd/sqn, full['racs']-ref['racs']+rrtmgstd/sqn, alpha=0.2, color='k', edgecolor=None)

    ax.plot(time, gcm_forc['ra']-gcm_clim['ra'], color='tab:gray', label='$\Delta R_a$')
    ax.plot(time, dsw, color='cyan', label='$\Delta SW$')
    ax.plot(time, dlwcs, '-', color='tab:red', label='$\Delta LW_{clear}$')
    ax.plot(time, dlwcld, '-', color='tab:purple', label='$\Delta LW_{cloud}$')
    ax.plot(time, full['racs']-ref['racs'], color='black', label='$\Delta$RRTMG$_{clear}$')
    # ax.set_xlim([time[0], time[-1]])
    ax.set_xlim([1987, 2200])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta$ Energy flux (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,2.5)
    leg=plt.legend(loc='lower left', prop={'size':8})
    leg.get_frame().set_linewidth(0)
    leg.get_frame().set_facecolor('none')
    if model=='mmm':
        fig.set_size_inches(4,3.3)
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    # no leg
    plotname='%s/dra_comp_noleg.cmip6.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if model == 'mmm':
        rastd = gcm_forc_mmm['std']['ra']-gcm_clim_mmm['std']['ra']
        ax.fill_between(time, gcm_forc['ra']-gcm_clim['ra'] - rastd, gcm_forc['ra']-gcm_clim['ra'] + rastd, alpha=0.2, color='k', edgecolor=None)
        ax.fill_between(time, dsw-dswstd, dsw+dswstd, alpha=0.2, color='cyan', edgecolor=None)
        ax.fill_between(time, dlwcs-dlwcsstd, dlwcs+dlwcsstd, alpha=0.2, color='tab:red', edgecolor=None)
        ax.fill_between(time, dlwcld-dlwcldstd, dlwcld+dlwcldstd, alpha=0.2, color='tab:purple', edgecolor=None)
        ax.fill_between(time, gcm_forc['ra']-gcm_clim['ra'] - rastd, gcm_forc['ra']-gcm_clim['ra'] + rastd, alpha=0.2, color='k', edgecolor=None)

        rrtmgstd = full_mmm['std'] - ref_mmm['std']
        ax.fill_between(time, full['racs']-ref['racs']-rrtmgstd, full['racs']-ref['racs']+rrtmgstd, alpha=0.2, color='k', edgecolor=None)

    ax.plot(time, gcm_forc['ra']-gcm_clim['ra'], color='tab:gray', label='$\Delta R_a$')
    ax.plot(time, dsw, color='cyan', label='$\Delta SW$')
    ax.plot(time, dlwcs, '-', color='tab:red', label='$\Delta LW_{clear}$')
    ax.plot(time, dlwcld, '-', color='tab:purple', label='$\Delta LW_{cloud}$')
    ax.plot(time, full['racs']-ref['racs'], color='black', label='$\Delta$RRTMG$_{clear}$')
    # ax.set_xlim([time[0], time[-1]])
    ax.set_xlim([1987, 2200])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta$ Energy flux (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,2.5)
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    # no leg
    plotname='%s/dra_comp_noleg_sw.cmip6.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if model == 'mmm':
        rastd = gcm_forc_mmm['std']['ra']-gcm_clim_mmm['std']['ra']
        ax.fill_between(time, gcm_forc['ra']-gcm_clim['ra'] - rastd, gcm_forc['ra']-gcm_clim['ra'] + rastd, alpha=0.2, color='k', edgecolor=None)
        # ax.fill_between(time, dsw-dswstd, dsw+dswstd, alpha=0.2, color='cyan', edgecolor=None)
        # ax.fill_between(time, dlwcs-dlwcsstd, dlwcs+dlwcsstd, alpha=0.2, color='tab:red', edgecolor=None)
        # ax.fill_between(time, dlwcld-dlwcldstd, dlwcld+dlwcldstd, alpha=0.2, color='tab:purple', edgecolor=None)
        # ax.fill_between(time, gcm_forc['ra']-gcm_clim['ra'] - rastd, gcm_forc['ra']-gcm_clim['ra'] + rastd, alpha=0.2, color='k', edgecolor=None)

        # rrtmgstd = full_mmm['std'] - ref_mmm['std']
        # ax.fill_between(time, full['racs']-ref['racs']-rrtmgstd, full['racs']-ref['racs']+rrtmgstd, alpha=0.2, color='k', edgecolor=None)

    ax.plot(time, gcm_forc['ra']-gcm_clim['ra'], color='tab:gray', label='$\Delta R_a$')
    # ax.plot(time, dsw, color='cyan', label='$\Delta SW$')
    # ax.plot(time, dlwcs, '-', color='tab:red', label='$\Delta LW_{clear}$')
    # ax.plot(time, dlwcld, '-', color='tab:purple', label='$\Delta LW_{cloud}$')
    # ax.plot(time, full['racs']-ref['racs'], color='black', label='$\Delta$RRTMG$_{clear}$')
    # ax.set_xlim([time[0], time[-1]])
    ax.set_xlim([1987, 2200])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta$ Energy flux (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,3.3)
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()


    ##################################################
    # TOTAL RA
    ##################################################
    # compare absolute change in radiative cooling
    plotname='%s/ra.cmip6.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, gcm_forc['ra'], color='tab:gray', label=modelstr)
    ax.plot(time, full['ra'], color='black', label='RRTMG')
    # ax.set_xlim([time[0], time[-1]])
    ax.set_xlim([1987, 2200])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$R_a$ (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,3.3)
    if model=='mmm':
        plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    # compare relative (to historical) change in radiative cooling
    plotname='%s/dra.cmip6.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.axhline(0,color='k', linewidth=0.5)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, gcm_forc['ra'] - gcm_clim['ra'], color='tab:gray', label=modelstr)
    ax.plot(time, full['ra'] - ref['ra'], color='black', label='RRTMG')
    # ax.set_xlim([time[0], time[-1]])
    ax.set_xlim([1987, 2200])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta R_a$ (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,3.3)
    if model=='mmm':
        plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # CLEAR SKY RA
    ##################################################
    # compare absolute change in radiative cooling
    plotname='%s/racs.cmip6.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, gcm_forc['racs'], color='tab:red', label=modelstr)
    ax.plot(time, full['racs'], color='black', label='RRTMG')
    # ax.set_xlim([time[0], time[-1]])
    ax.set_xlim([1987, 2200])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$R_{a,\,cs}$ (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,3.3)
    if model=='mmm':
        plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()


    # compare relative (to historical) change in clear sky radiative cooling
    plotname='%s/dracs.cmip6.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.axhline(0,color='k', linewidth=0.5)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, gcm_forc['racs'] - gcm_clim['racs'], color='tab:red', label=modelstr)
    # ax.plot(time, full['racs'] - gcm_clim['racs'], color='black', label='RRTMG')
    ax.plot(time, full['racs'] - ref['racs'], color='black', label='RRTMG')
    # ax.set_xlim([time[0], time[-1]])
    ax.set_xlim([1987, 2200])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,3.3)
    if model=='mmm':
        plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

