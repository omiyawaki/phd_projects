import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import get_datadir, get_plotdir
from misc.filenames import *
from proc.hydro import save_hydro
from plot.titles import make_title_sim_time
import os
import pickle
import numpy as np
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
# import tikzplotlib

def hydro_mon_lat(sim, **kwargs):

    categ = 'mon_lat'

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    vertcoord = kwargs.get('vertcoord', 'si') # vertical coordinate (si for sigma, pa for pressure, z for height)
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute hydro
    vertbnd = kwargs.get('vertbnd', (0.7, 0.3)) # sigma bounds of vertical integprl

    if sim == 'longrun':
        model = kwargs.get('model', 'MPIESM12_abrupt4x')
        yr_span = kwargs.get('yr_span', '1000')
        yr_base = 0
    elif sim == 'rcp85':
        model = kwargs.get('model', 'MPI-ESM-LR')
        yr_span = kwargs.get('yr_span', '200601-230012')
        yr_base = 0 if 'ymonmean' in timemean else 2006
    elif sim == 'echam':
        model = kwargs.get('model', 'rp000140')
        yr_span = kwargs.get('yr_span', '0001_0039')
        yr_base = 0
    elif sim == 'era5':
        model = None
        yr_span = kwargs.get('yr_span', '1980_2005')
        yr_base = 0 if 'ymonmean' in timemean else 1980

    vmin =-10 
    vmax = 10
    vint = 0.2

    # load data and plot directories
    datadir = get_datadir(sim, model=model, yr_span=yr_span)
    plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

    # location of pickled hydro data
    hydro_file = '%s/hydro.%s.%s.pickle' % (datadir, zonmean, timemean)

    if not (os.path.isfile(hydro_file) and try_load):
        save_hydro(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

    [hydro, grid] = pickle.load(open(hydro_file, 'rb'))
    pr = hydro['pr']
    prc = hydro['prc']

    rolling_mean = 0; # smooth data using a rolling mean? (units: yr)
    pr_filt = uniform_filter(pr, [rolling_mean,0]) # apply rolling mean

    [mesh_lat, mesh_time] = np.meshgrid(grid['lat'], np.arange(pr.shape[0])) # create mesh

    plotname = '%s/pr_mon_lat%s' % (plotdir, timemean)
    fig, ax = plt.subplots()
    csf = ax.contourf(mesh_time, mesh_lat, pr, np.arange(vmin,vmax,vint), cmap='RdBu_r', vmin=vmin, vmax=vmax, extend='both')
    # cs_0 = ax.contour(mesh_time, mesh_lat, pr, levels=[0], colors='gray', linewidths=3)
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
    ax.yaxis.set_minor_locator(MultipleLocator(10))
    cbar = plt.colorbar(csf)
    cbar.set_label('$P$ (mm d$^{-1}$)')
    # plt.savefig('%s.png' % (plotname), dpi=300)
    fig.set_size_inches(5, 4)
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    plt.show()

    plotname = '%s/prc_mon_lat%s' % (plotdir, timemean)
    fig, ax = plt.subplots()
    csf = ax.contourf(mesh_time, mesh_lat, prc, np.arange(vmin,vmax,vint), cmap='RdBu_r', vmin=vmin, vmax=vmax, extend='both')
    # cs_0 = ax.contour(mesh_time, mesh_lat, prc, levels=[0], colors='gray', linewidths=3)
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
    ax.yaxis.set_minor_locator(MultipleLocator(10))
    cbar = plt.colorbar(csf)
    cbar.set_label('$P_c$ (mm d$^{-1}$)')
    # plt.savefig('%s.png' % (plotname), dpi=300)
    fig.set_size_inches(5, 4)
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    plt.show()

    # tikzplotlib.save('%s.tex' % (plotname))
