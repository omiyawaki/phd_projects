import os
import sys
import pickle
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from tools import get_modelidx, get_modelstr

# models = ['mmm']
# models=['CCSM4']
# models=['bcc-csm1-1', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MPI-ESM-LR']
models=['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MPI-ESM-LR']

for model in models:

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

    if model == 'mmm':
        gcm_clim = mmm_clim
        gcm_forc = mmm_forc
    else:
        imod = get_modelidx(model)
        gcm_clim = {}
        gcm_forc = {}
        for varname in indiv_clim:
            gcm_clim[varname] = indiv_clim[varname][imod,...]
            gcm_forc[varname] = indiv_forc[varname][imod,...]

    ##################################################
    # LOAD RRTMG DATA
    ##################################################
    # RRTMG climatology
    ref = {}
    [ref['ra'], ref['racs'], ref['diag']] = pickle.load(open('%s/rad.hist.%s.pickle' % (datadir, model), 'rb'))

    # RRTMG full RCP8.5 response
    full = {}
    [full['ra'], full['racs'], full['diag']] = pickle.load(open('%s/rad.rcp.%s.pickle' % (datadir, model), 'rb'))

    ##################################################
    # MISC data
    ##################################################
    plev = 1e-2*grid['lev']
    modelstr = get_modelstr(model)

    ##################################################
    # HISTORICAL
    ##################################################
    time0 = np.arange(1860,2007)

    ##################################################
    # TOTAL RA
    ##################################################
    # compare absolute change in radiative cooling
    plotname='%s/ra_hist.pdf' % (plotdir)
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
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # CLEAR SKY RA
    ##################################################
    # compare absolute change in radiative cooling
    plotname='%s/racs_hist.pdf' % (plotdir)
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
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # RCP
    ##################################################
    time = np.arange(2006,2301)

    for varname in gcm_clim:
        gcm_clim[varname] = np.nanmean(gcm_clim[varname][-30:,...], axis=0)

    ref['ra'] = np.nanmean(ref['ra'][-30:])
    ref['racs'] = np.nanmean(ref['racs'][-30:])

    ##################################################
    # ALL COMP DRA
    ##################################################
    # compare relative change in radiative cooling
    plotname='%s/dra_comp.pdf' % (plotdir)
    sw_forc = gcm_forc['rsdt'] - gcm_forc['rsut'] + gcm_forc['rsus'] - gcm_forc['rsds']
    lw_forc = - gcm_forc['rlut'] + gcm_forc['rlus'] - gcm_forc['rlds']
    lwcs_forc = - gcm_forc['rlutcs'] + gcm_forc['rlus'] - gcm_forc['rldscs']
    lwcld_forc = lw_forc - lwcs_forc
    sw_clim = gcm_clim['rsdt'] - gcm_clim['rsut'] + gcm_clim['rsus'] - gcm_clim['rsds']
    lw_clim = - gcm_clim['rlut'] + gcm_clim['rlus'] - gcm_clim['rlds']
    lwcs_clim = - gcm_clim['rlutcs'] + gcm_clim['rlus'] - gcm_clim['rldscs']
    lwcld_clim = lw_clim - lwcs_clim
    dsw = sw_forc - sw_clim
    dlw = lw_forc - lw_clim
    dlwcs = lwcs_forc - lwcs_clim
    dlwcld = lwcld_forc - lwcld_clim
    fig, ax = plt.subplots()
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, gcm_forc['ra']-gcm_clim['ra'], color='tab:gray', label='$\Delta R_a$')
    ax.plot(time, dsw, color='cyan', label='$\Delta SW$')
    ax.plot(time, dlwcs, '--', color='tab:red', label='$\Delta LW_{clear}$')
    ax.plot(time, dlwcld, ':', color='tab:red', label='$\Delta LW_{cloud}$')
    ax.plot(time, full['racs']-ref['racs'], color='black', label='$\Delta$RRTMG$_{clear}$')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta$ Energy flux (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,3.3)
    leg=plt.legend()
    leg.get_frame().set_linewidth(0)
    leg.get_frame().set_facecolor('none')
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # TOTAL RA
    ##################################################
    # compare absolute change in radiative cooling
    plotname='%s/ra.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, gcm_forc['ra'], color='tab:gray', label=modelstr)
    ax.plot(time, full['ra'], color='black', label='RRTMG')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$R_a$ (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,3.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    # compare relative (to historical) change in radiative cooling
    plotname='%s/dra.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.axhline(0,color='k', linewidth=0.5)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, gcm_forc['ra'] - gcm_clim['ra'], color='tab:gray', label=modelstr)
    ax.plot(time, full['ra'] - ref['ra'], color='black', label='RRTMG')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta R_a$ (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,3.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

    ##################################################
    # CLEAR SKY RA
    ##################################################
    # compare absolute change in radiative cooling
    plotname='%s/racs.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, gcm_forc['racs'], color='tab:red', label=modelstr)
    ax.plot(time, full['racs'], color='black', label='RRTMG')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$R_{a,\,cs}$ (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,3.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()


    # compare relative (to historical) change in clear sky radiative cooling
    plotname='%s/dracs.pdf' % (plotdir)
    fig, ax = plt.subplots()
    ax.axhline(0,color='k', linewidth=0.5)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.plot(time, gcm_forc['racs'] - gcm_clim['racs'], color='tab:red', label=modelstr)
    # ax.plot(time, full['racs'] - gcm_clim['racs'], color='black', label='RRTMG')
    ax.plot(time, full['racs'] - ref['racs'], color='black', label='RRTMG')
    ax.set_xlim([time[0], time[-1]])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('Time (yr)')
    ax.set_ylabel(r'$\Delta R_{a,\,cs}$ (W m$^{-2}$)')
    ax.set_title(modelstr)
    fig.set_size_inches(4,3.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotname, format='pdf', dpi=300)
    plt.close()

