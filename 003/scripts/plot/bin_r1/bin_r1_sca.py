import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import *
from misc.filenames import *
# from misc.translate import translate_varname
from misc.means import lat_mean, global_int
from misc.load_data import *
from misc import par
from proc.r1 import save_r1
from plot.titles import make_title_sim_time_lat
from yaxis_bin_r1_sca import yaxis_bin_r1_sca
import os
import pickle
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.ndimage import uniform_filter
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

def bin_r1_sca(sim, **kwargs):

    categ = 'bin_r1'

    annmean = kwargs.get('annmean', 0)
    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    latbnd = kwargs.get('latbnd', (80,90))
    latstep = kwargs.get('latstep', 0.25) # latitude step size used for interpolation
    plotvar = kwargs.get('plotvar', None) # plot overlay (sic for sea ice, ga_dev for lapse rate deviation, pr for precip, decomp for linear decomposition)?
    refclim = kwargs.get('refclim', 'hist-30') # reference climate from which to compute deviations (init is first time step, hist-30 is the last 30 years of the historical run)
    legend = kwargs.get('legend', 0)

    if sim == 'longrun':
        model = kwargs.get('model', 'MPIESM12_abrupt4x')
        yr_span = kwargs.get('yr_span', '1000')
        yr_base = 0
    elif sim == 'rcp85':
        model = kwargs.get('model', 'MPI-ESM-LR')
        yr_span = kwargs.get('yr_span', '200601-230012')
        if 'ymonmean' not in timemean:
            yr_base = 2006
        else:
            yr_base = 0
    elif sim == 'historical':
        model = kwargs.get('model', 'MPI-ESM-LR')
        yr_span = kwargs.get('yr_span', '186001-200512')
        if 'ymonmean' not in timemean:
            yr_base = 1860
        else:
            yr_base = 0
    elif sim == 'echam':
        model = kwargs.get('model', 'rp000140')
        yr_span = kwargs.get('yr_span', '0001_0039')
        yr_base = 0
    elif sim == 'era5':
        model = None
        yr_span = kwargs.get('yr_span', '1979_2019')
        if 'ymonmean' not in timemean:
            yr_base = int(yr_span[0:4])
        else:
            yr_base = 0

    ##########################################
    ## X AXIS SPECIFICATIONS
    ##########################################
    vmin, vmax, ylabel = yaxis_bin_r1_sca(sim, timemean=timemean, latbnd=latbnd, plotvar=plotvar)

    ##################################
    # LOAD DATA
    ##################################
    if isinstance(model, str) or model is None:
        [r1, grid, datadir, plotdir, modelstr] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 
    else:
        [r1, grid, datadir, plotdir, modelstr, r1_mmm] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span, refclim=refclim) 

    # load total precipitation
    if plotvar == 'pr':
        if isinstance(model, str) or model is None:
            [hydro, grid, datadir, plotdir, modelstr] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
        else:
            [hydro, grid, datadir, plotdir, modelstr, hydro_mmm] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

        # save pr
        binvar = hydro['pr']

        # if multimodel mean save statistics
        if not ( isinstance(model, str) or model is None ):
            binvar_mmm = hydro_mmm['pr']

    # load convection precipitation
    if plotvar == 'prc':
        if isinstance(model, str) or model is None:
            [hydro, grid, datadir, plotdir, modelstr] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
        else:
            [hydro, grid, datadir, plotdir, modelstr, hydro_mmm] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

        # save prc
        binvar = hydro['prc']

        # if multimodel mean save statistics
        if not ( isinstance(model, str) or model is None ):
            binvar_mmm = hydro_mmm['prc']

    # convective precipitation fraction
    if plotvar == 'prfrac':
        if isinstance(model, str) or model is None:
            [hydro, grid, datadir, plotdir, modelstr] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)
        else:
            [hydro, grid, datadir, plotdir, modelstr, hydro_mmm] = load_hydro(sim, categ, zonmean=zonmean, timemean=timemean, try_load=try_load, model=model, yr_span=yr_span)

        # save prfrac
        binvar = hydro['prfrac']

        # if multimodel mean save statistics
        if not ( isinstance(model, str) or model is None ):
            binvar_mmm = hydro_mmm['prfrac']

    # take annual mean?
    if annmean:
        r1 = np.mean(r1, axis=0, keepdims=True)
        binvar = np.mean(binvar, axis=0, keepdims=True)

    # limit data to select region?
    lat_sel = np.where( np.logical_and(grid['lat'] >= latbnd[0], grid['lat'] <= latbnd[1]) )
    r1 = r1[:,lat_sel[0]]
    binvar = binvar[:,lat_sel[0]]

    # flatten data
    r1_flat = np.ndarray.flatten(r1)
    binvar_flat = np.ndarray.flatten(binvar)

    ############################################
    # SCATTER PLOT
    ############################################

    if plotvar == 'prfrac':
        # print(binvar_avg[par.r1_bins < 1])
        # r1_sel = np.squeeze(r1_avg[par.r1_bins < 1])
        # binvar_sel = np.squeeze(binvar_avg[par.r1_bins < 1])
        A = np.vstack([np.ndarray.flatten(r1), np.ones(r1.size)]).T
        m, c = np.linalg.lstsq(A, np.ndarray.flatten(binvar), rcond=None)[0]
        r1_reg = np.linspace(r1.min(),r1.max(),21)

    if annmean:
        plotname = '%s/%s.scat.annmean.%.0f.%.0f.%s' % (plotdir, plotvar, latbnd[0], landbnd[1], timemean)
    else:
        plotname = '%s/%s.scat.%.0f.%.0f.%s' % (plotdir, plotvar, latbnd[0], latbnd[1], timemean)

    fig, ax = plt.subplots()
    ax.axhline(0,0.8,1.1, color='k', linewidth=0.5)
    # ax.plot(r1_flat, binvar_flat, '.k')
    # hsv = cm.get_cmap('hsv', 12)
    hsv = cm.get_cmap('twilight', 12)
    monstr = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
    # hleg = {}
    for mon in range(12):
        labelstr = monstr[mon]
        # hleg[mon] = ax.plot(r1[mon,0], binvar[mon,0], '.', color=hsv(mon), label=labelstr)
        ax.plot(r1[mon,0], binvar[mon,0], '.', color=hsv(mon), label=labelstr)
    for mon in range(12):
        ax.plot(r1[mon,:], binvar[mon,:], '.', color=hsv(mon))
    ax.set_title('%.0f$^\circ$ to %.0f$^\circ$' % (latbnd[0], latbnd[1]))
    ax.axhline(0,0,1, color='k', linewidth=0.5)
    if plotvar == 'prfrac':
        ax.plot(r1_reg, c + m*r1_reg, '--k', label='$P_c/P = %.2f R_1 + %.2f$' % (m, c))
    ax.legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('$R_1$ (unitless)')
    ax.set_ylabel(ylabel)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300, bbox_inches="tight")

    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # SCATTER PLOT (SKIP months 5 and 6)
    ############################################

    # compute linear regression through prfrac=0, r1=1

    skip56 = [0,1,2,3,6,7,8,9,10,11]
    r1_skip56 = r1[skip56,:]
    binvar_skip56 = binvar[skip56,:]
    
    if plotvar == 'prfrac':
        # print(binvar_avg[par.r1_bins < 1])
        # r1_sel = np.squeeze(r1_avg[par.r1_bins < 1])
        # binvar_sel = np.squeeze(binvar_avg[par.r1_bins < 1])
        A = np.vstack([np.ndarray.flatten(r1_skip56), np.ones(r1_skip56.size)]).T
        m, c = np.linalg.lstsq(A, np.ndarray.flatten(binvar_skip56), rcond=None)[0]
        r1_reg = np.linspace(r1.min(),r1.max(),21)
        
    if annmean:
        plotname = '%s/%s.scat.skip56.annmean.%.0f.%.0f.%s' % (plotdir, plotvar, latbnd[0], landbnd[1], timemean)
    else:
        plotname = '%s/%s.scat.skip56.%.0f.%.0f.%s' % (plotdir, plotvar, latbnd[0], latbnd[1], timemean)

    fig, ax = plt.subplots()
    hsv = cm.get_cmap('twilight', 12)
    monstr = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
    for mon in range(12):
        if mon in [4,5]:
            continue
        labelstr = monstr[mon]
        ax.plot(r1[mon,0], binvar[mon,0], '.', color=hsv(mon), label=labelstr)
    for mon in range(12):
        if mon in [4,5]:
            continue
        ax.plot(r1[mon,:], binvar[mon,:], '.', color=hsv(mon))
        
    ax.set_title('%.0f$^\circ$ to %.0f$^\circ$' % (latbnd[0], latbnd[1]))
    ax.axhline(0,0,1, color='k', linewidth=0.5)
    if plotvar == 'prfrac':
        ax.plot(r1_reg, c + m*r1_reg, '--k', label='$P_c/P = %.2f R_1 + %.2f$' % (m, c))
    ax.legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xlabel('$R_1$ (unitless)')
    ax.set_ylabel(ylabel)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300, bbox_inches="tight")

    if viewplt:
        plt.show()
    plt.close()

    ############################################
    # BIN VARIABLE ACCORDING TO R1
    ############################################

    # categorize r1 data into bins
    id_bins = np.digitize(r1, par.r1_bins)

    # create array to save binned averages
    r1_avg = np.empty([len(par.r1_bins), 1])
    binvar_avg = np.empty([len(par.r1_bins), 1])

    # create cos(lat) data for area averaging
    clat = np.cos( np.deg2rad( np.tile(grid['lat'][lat_sel[0]], (r1.shape[0], 1)) ) )

    for bin in range(len(par.r1_bins)):
        idx_bins = np.where( id_bins == bin)
        r1_avg[bin] = np.sum( clat[idx_bins]*r1[idx_bins] ) / np.sum(clat[idx_bins])
        binvar_avg[bin] = np.sum( clat[idx_bins] * binvar[idx_bins] ) / np.sum(clat[idx_bins])

    # for multimodel mean compute the spread
    if not ( isinstance(model, str) or model is None ):
        binvar_p25 = np.empty([len(par.r1_bins), 1])
        binvar_p75 = np.empty([len(par.r1_bins), 1])
        binvar_std = np.empty([len(par.r1_bins), 1])

        for bin in range(len(par.r1_bins)):
            idx_bins = np.where( id_bins == bin)
            # r1_avg[bin] = np.sum( clat[idx_bins]*r1[idx_bins] ) / np.sum(clat[idx_bins])
            binvar_p25[bin] = np.sum( clat[idx_bins] * binvar_mmm['prc25'][idx_bins] ) / np.sum(clat[idx_bins])
            binvar_p75[bin] = np.sum( clat[idx_bins] * binvar_mmm['prc75'][idx_bins] ) / np.sum(clat[idx_bins])
            binvar_std[bin] = np.sum( clat[idx_bins] * binvar_mmm['std'][idx_bins] ) / np.sum(clat[idx_bins])

    ############################################
    # for prfrac, fit linear regression through prfrac=0, r1=1
    ############################################
    
    if plotvar == 'prfrac':
        # print(binvar_avg[par.r1_bins < 1])
        r1_sel = np.squeeze(r1_avg[par.r1_bins < 1])
        binvar_sel = np.squeeze(binvar_avg[par.r1_bins < 1])
        A = np.vstack([r1_sel, np.ones(len(r1_sel))]).T
        m, c = np.linalg.lstsq(A, binvar_sel, rcond=None)[0]
        
    ############################################
    # DISTRIBUTION PLOT
    ############################################

    if annmean:
        plotname = '%s/%s.dist.annmean.%.0f.%.0f.%s' % (plotdir, plotvar, latbnd[0], latbnd[1], timemean)
    else:
        plotname = '%s/%s.dist.%.0f.%.0f.%s' % (plotdir, plotvar, latbnd[0], latbnd[1], timemean)

    fig, ax = plt.subplots()
    rae = patches.Rectangle((0.9,vmin),1.5-0.9,vmax-vmin, alpha=0.2, facecolor='tab:blue', edgecolor=None)
    rce = patches.Rectangle((-0.6,vmin),0.1-(-0.6),vmax-vmin, alpha=0.2, facecolor='tab:orange', edgecolor=None)
    ax.add_patch(rae)
    ax.add_patch(rce)
    ax.axhline(0,0,1, color='k', linewidth=0.5)
    ax.plot(r1_avg, binvar_avg, '-k')
    if not ( isinstance(model, str) or model is None ):
        # ax.fill_between(np.squeeze(r1_avg), np.squeeze(binvar_p25), np.squeeze(binvar_p75), facecolor='k', alpha=0.2, edgecolor=None)
        ax.fill_between(np.squeeze(r1_avg), np.squeeze(binvar_avg-binvar_std), np.squeeze(binvar_avg+binvar_std), facecolor='k', alpha=0.2, edgecolor=None)
    if plotvar == 'prfrac':
        ax.plot(r1_avg, r1_avg*m + c, '--k', label='$P_c/P = %.2f R_1 + %.2f$' % (m, c))
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    ax.set_xticks(np.arange(-70,150,20)/100)
    ax.set_xlim(par.r1_bins.min()-0.05, par.r1_bins.max()-0.05)
    ax.set_ylim(vmin,vmax)
    ax.set_xlabel('$R_1$ (unitless)')
    ax.set_ylabel(ylabel)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    plt.legend()
    fig.set_size_inches(5, 4)
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300, bbox_inches="tight")

    if viewplt:
        plt.show()
    plt.close()

