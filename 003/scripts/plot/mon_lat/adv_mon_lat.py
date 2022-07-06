import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import get_datadir, get_plotdir
from misc.filenames import *
from proc.r1 import save_r1
from plot.titles import make_title_sim_time
import os
import pickle
import numpy as np
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
# import tikzplotlib

def adv_mon_lat(sim, **kwargs):

    categ = 'mon_lat'

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    refclim = kwargs.get('refclim', 'hist-30') # reference climate from which to compute deviations (init is first time step, hist-30 is the last 30 years of the historical run)
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)
    sim_ref = kwargs.get('sim_ref', 'historical')
    timemean_ref = kwargs.get('timemean_ref', 'ymonmean-30')
    yr_span_ref = kwargs.get('yr_span_ref', '186001-200512')

    if sim == 'longrun':
        model = kwargs.get('model', 'MPIESM12_abrupt4x')
        yr_span = kwargs.get('yr_span', '1000')
        yr_base = 0
    elif sim == 'rcp85':
        model = kwargs.get('model', 'MPI-ESM-LR')
        yr_span = kwargs.get('yr_span', '200601-230012')
        yr_base = 2006
    elif sim == 'echam':
        model = kwargs.get('model', 'rp000140')
        yr_span = kwargs.get('yr_span', '0001_0039')
        yr_base = 0
    elif sim == 'era5':
        model = None
        yr_span = kwargs.get('yr_span', '1980_2005')
        yr_base = 1980

# load data and plot directories
    datadir = get_datadir(sim, model=model, yr_span=yr_span)
    plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

# location of pickled advection data
    stg_adv_file = remove_repdots('%s/stg_adv.%s.%s.pickle' % (datadir, zonmean, timemean))

    if not (os.path.isfile(stg_adv_file) and try_load):
        save_r1(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

    [stg_adv, grid] = pickle.load(open(stg_adv_file, 'rb'))

    [mesh_lat, mesh_time] = np.meshgrid(grid['lat'], np.arange(stg_adv.shape[0])) # create mesh

# compute deviation from control climate
    if refclim == 'hist-30':
        datadir_ref = get_datadir(sim_ref, model=model, yr_span=yr_span_ref)

        # location of pickled historical stg_adv data
        stg_adv_file_ref = remove_repdots('%s/stg_adv.%s.%s.pickle' % (datadir_ref, zonmean, timemean_ref))

        if not (os.path.isfile(stg_adv_file_ref) and try_load):
            save_r1(sim_ref, model=model, zonmean=zonmean, timemean=timemean_ref, yr_span=yr_span_ref, refclim='init')

        [stg_adv_ref, grid_ref] = pickle.load(open(stg_adv_file_ref, 'rb'))
        
        if timemean == 'djfmean':
            stg_adv_ref = np.mean(np.roll(stg_adv_ref,1,axis=0)[0:3], 0)
        elif timemean == 'jjamean':
            stg_adv_ref = np.mean(stg_adv_ref[5:8], 0)

    stg_adv_dev = stg_adv - stg_adv_ref

    ############################################
    # PLOT (ABS VALUE)
    ############################################
    plotname = remove_repdots('%s/stg_adv_mon_lat.%s' % (plotdir, timemean))

    fig, ax = plt.subplots()
    vmin = -150
    vmax = 150 
    csf = ax.contourf(mesh_time, mesh_lat, stg_adv, np.arange(vmin,vmax,10), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    make_title_sim_time(ax, sim, model=model, timemean=timemean)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Year')
    ax.set_ylabel('Latitude (deg)')
    ax.set_yticks(np.arange(-90,91,30))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    cbar = plt.colorbar(csf)
    cbar.set_label('$\partial_t m + \partial_y (vm)$ (Wm$^{-2}$)')
# plt.savefig('%s.png' % (plotname), dpi=300)
    plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    if viewplt:
        plt.show()

    ############################################
    # PLOT (DEVIATION FROM INITIAL)
    ############################################
    plotname = remove_repdots('%s/stg_adv_dev_mon_lat.%s' % (plotdir, timemean))

    fig, ax = plt.subplots()
    vmin = -50
    vmax = 50 
    csf = ax.contourf(mesh_time, mesh_lat, stg_adv_dev, np.arange(vmin,vmax,10), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    make_title_sim_time(ax, sim, model=model, timemean=timemean)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Year')
    ax.set_ylabel('Latitude (deg)')
    ax.set_yticks(np.arange(-90,91,30))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    cbar = plt.colorbar(csf)
    cbar.set_label('$\Delta(\partial_t m + \partial_y (vm))$ (Wm$^{-2}$)')
# plt.savefig('%s.png' % (plotname), dpi=300)
    plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    if viewplt:
        plt.show()
