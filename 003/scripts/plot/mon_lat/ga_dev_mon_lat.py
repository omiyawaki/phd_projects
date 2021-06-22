import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import get_datadir, get_plotdir
from misc.filenames import remove_repdots
from proc.ga import make_ga_dev_vint
from plot.titles import make_title_sim_time
import os
import pickle
import numpy as np
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
# import tikzplotlib

def ga_dev_mon_lat(sim, **kwargs):

    categ = 'mon_lat'

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (.yearmean, .jjamean, .djfmean, .ymonmean-30)
    vertcoord = kwargs.get('vertcoord', 'si') # vertical coordinate (si for sigma, pa for pressure, z for height)
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    vertbnd = kwargs.get('vertbnd', (0.7, 0.3)) # sigma bounds of vertical integral

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

    ga_dev_vint = make_ga_dev_vint(sim, vertbnd, model=model, vertcoord = vertcoord, zonmean=zonmean, timemean=timemean, yr_span=yr_span, try_load=try_load)

    # print(np.reshape(ga_dev, (-1,96,12)).shape)
    if timemean == '':
        ga_dev_vint['ga_dev_vint'] = np.mean(np.reshape(ga_dev_vint['ga_dev_vint'], (-1,12,ga_dev_vint['ga_dev_vint'].shape[1])),1)

    rolling_mean = 0; # smooth data using a rolling mean? (units: yr)
    ga_dev_vint_filt = uniform_filter(ga_dev_vint['ga_dev_vint'], [rolling_mean,0]) # apply rolling mean

    [mesh_lat, mesh_time] = np.meshgrid(ga_dev_vint['grid']['lat'], yr_base + np.arange(ga_dev_vint['ga_dev_vint'].shape[0])) # create mesh

    ##################################
    # PLOT
    ##################################
    plotname = '%s/ga_dev_mon_lat.%g.%g.%s' % (plotdir, vertbnd[0], vertbnd[1], timemean)
    fig, ax = plt.subplots()
    vmin = -50
    vmax = 50
    csf = ax.contourf(mesh_time, mesh_lat, ga_dev_vint['ga_dev_vint'], np.arange(vmin,vmax,1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Year')
    # cbar = plt.colorbar(csf)
    # cbar.set_label('$R_1$ (unitless)')
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    plt.show()
