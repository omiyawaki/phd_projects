import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import *
from misc.filenames import *
from misc.means import *
from misc.load_data import load_r1
import misc.par as par
from proc.r1 import save_r1
from plot.titles import make_title_sim_time
import os
import pickle
import numpy as np
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
# import tikzplotlib

def r1_mon_lat(sim, **kwargs):

    categ = 'mon_lat'

    zonmean = kwargs.get('zonmean', 'zonmean') # zonal mean?
    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    domain = kwargs.get('domain', '')
    try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
    viewplt = kwargs.get('viewplt', 0) # view plot? (plt.show)

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

    if domain == '':
        vmin = -90
        vmax = 90
    elif domain == 'rce_nh':
        if timemean == 'yearmean':
            vmin = 40
            vmax = 45
        elif timemean == 'djfmean':
            vmin = 38 
            vmax = 43
        elif timemean == 'jjamean':
            vmin = 40
            vmax = 60
    elif domain == 'rce_sh':
        if timemean == 'yearmean':
            vmin = -38
            vmax = -43
        elif timemean == 'djfmean':
            vmin = -38
            vmax = -43
        elif timemean == 'jjamean':
            vmin = -38
            vmax = -43
    elif domain == 'rae_nh':
        if timemean == 'yearmean':
            vmin = 70
            vmax = 90
        elif timemean == 'djfmean':
            vmin = 70
            vmax = 90
        elif timemean == 'jjamean':
            vmin = 70
            vmax = 90
    elif domain == 'rae_sh':
        if timemean == 'yearmean':
            vmin = -60
            vmax = -80
        elif timemean == 'djfmean':
            vmin = -60
            vmax = -80
        elif timemean == 'jjamean':
            vmin = -60
            vmax = -80

    ##################################
    # LOAD DATA
    ##################################
    if isinstance(model, str):
        [r1, grid, datadir, plotdir, modelstr] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, domain=domain, try_load=try_load, model=model, yr_span=yr_span) 
    else:
        [r1, grid, datadir, plotdir, modelstr, r1_mmm] = load_r1(sim, categ, zonmean=zonmean, timemean=timemean, domain=domain, try_load=try_load, model=model, yr_span=yr_span) 

    # print(np.reshape(r1, (-1,96,12)).shape)
    if timemean == '':
        r1 = np.mean(np.reshape(r1, (-1,12,r1.shape[1])),1)
        
    [mesh_lat, mesh_time] = np.meshgrid(grid['lat'], yr_base + np.arange(r1.shape[0])) # create mesh

    ##################################
    # PLOT
    ##################################
    plotname = remove_repdots('%s/r1_mon_lat.%s.%s' % (plotdir, domain, timemean))
    fig, ax = plt.subplots()
    vmin = -1.7
    vmax = 1.7
    csf = ax.contourf(mesh_time, mesh_lat, r1, np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
    cs_rae = ax.contour(mesh_time, mesh_lat, r1, levels=[0.9], colors='royalblue', linewidths=3)
    cs_rce = ax.contour(mesh_time, mesh_lat, r1, levels=[0.1], colors='sandybrown', linewidths=3)
    make_title_sim_time(ax, sim, model=modelstr, timemean=timemean)
    ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
    if 'ymonmean' in timemean:
        ax.set_xticks(np.arange(0,12,1))
        ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    else:
        ax.set_xlabel('Year')
    ax.set_ylim([vmin, vmax])
    ax.set_ylabel('Latitude (deg)')
    ax.set_yticks(np.arange(-90,91,30))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(MultipleLocator(10))
    cbar = plt.colorbar(csf)
    cbar.set_label('$R_1$ (unitless)')
    # plt.savefig('%s.png' % (plotname), dpi=300)
    plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
    if viewplt:
        plt.show()
    plt.close()