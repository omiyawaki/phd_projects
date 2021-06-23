import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import get_datadir, get_plotdir
from misc.filenames import remove_repdots
from proc.dTdt import make_dTdt_lev
from plot.titles import make_title_sim_time
import os
import pickle
import numpy as np
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
# import tikzplotlib

def dTdt_mon_lat(sim, **kwargs):
    # plot time evolution of temperature tendency evaluated at some level

    categ = 'mon_lat'

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    vertcoord = kwargs.get('vertcoord', 'si') # vertical coordinate (si for sigma, pa for pressure, z for height)
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    vertlev = kwargs.get('vertlev', 1) # vertical level to evaluate dTdt

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

    if timemean == 'yearmean' or timemean == '':
        vmin = -20
        vmax = 20

    # load data and plot directories
    datadir = get_datadir(sim, model=model, yr_span=yr_span)
    plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

    dTdt_lev = make_dTdt_lev(sim, vertlev, model=model, vertcoord = vertcoord, zonmean=zonmean, timemean=timemean, yr_span=yr_span, try_load=try_load)

    # print(np.reshape(dTdt, (-1,96,12)).shape)
    if timemean == '':
        dTdt_lev['dTdt_lev'] = np.mean(np.reshape(dTdt_lev['dTdt_lev'], (-1,12,dTdt_lev['dTdt_lev'].shape[1])),1)

    rolling_mean = 10; # smooth data using a rolling mean? (units: yr)
    dTdt_lev['dTdt_lev'] = uniform_filter(dTdt_lev['dTdt_lev'], [rolling_mean,0]) # apply rolling mean

    [mesh_lat, mesh_time] = np.meshgrid(dTdt_lev['grid']['lat'], yr_base + np.arange(dTdt_lev['dTdt_lev'].shape[0])) # create mesh

    ##################################
    # PLOT
    ##################################
    plotname = '%s/dTdt_mon_lat.%g.%s' % (plotdir, vertlev, timemean)
    fig, ax = plt.subplots()
    csf = ax.contourf(mesh_time, mesh_lat, dTdt_lev['dTdt_lev'], np.linspace(vmin,vmax,50), cmap='RdBu_r', vmin=vmin, vmax=vmax, extend='both')
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Year')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    make_title_sim_time(ax, sim, model=model, timemean=timemean)
    ax.set_ylabel('Latitude (deg)')
    ax.set_yticks(np.arange(-90,91,30))
    ax.yaxis.set_minor_locator(MultipleLocator(10))
    cbar = plt.colorbar(csf)
    # cbar.set_label(r'$\langle(\Gamma_m - \Gamma)/\Gamma_m\rangle_{%0.1f}^{%0.1f}$ (%%)' % (vertbnd[0], vertbnd[1]))
    plt.savefig(remove_repdots('%s.pdf' % (plotname)), format='pdf', dpi=300)
    plt.show()
