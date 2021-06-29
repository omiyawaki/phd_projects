import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import get_datadir, get_plotdir
from misc.filenames import filenames_raw
from misc.translate import *
from proc.r1 import save_r1
from plot.titles import make_title_sim_time_seas
import os
import pickle
import numpy as np
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from netCDF4 import Dataset
import tikzplotlib

def plot_plev3d_lat(sim, varname, **kwargs):
	
	categ = 'lat'

	timemean = kwargs.get('timemean', 'ymonmean-30') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
	seasmean = kwargs.get('seasmean', '')
	domain = kwargs.get('domain', '')
	try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1
	vertlev = kwargs.get('vertlev', 1) # vertical level to evaluate 

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
		
	if translate_varname(varname) == 'ta':
		ylabel = '$T$ (K)'
	
	# load data and plot directories
	plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)
	
	# location of raw 3d data
	file = filenames_raw(sim, varname, model=model, timemean=timemean, yr_span=yr_span)
	
	# select level
	var2d = np.squeeze(file.variables[varname][:][:,vertlev,:,:])
	plev = file.variables[latetrans_grid(sim, 'lev')][:][vertlev]
	
	# select seasons
	if seasmean=='djf':
		var2d = np.roll(var2d,1,axis=0)[:3,:,:]
	elif seasmean=='jja':
		var2d = var2d[5:8,:,:]
	elif seasmean=='':
		var2d = var2d[:]

	# take time and zonal mean
	var2d = np.mean(np.squeeze(var2d),(0,2))
		
	grid = {}
	grid['lat'] = file.variables['lat'][:]
		
	plotname = '%s/%s_lat.lev%g.%s.%s' % (plotdir, varname, vertlev, timemean, seasmean)
	fig, ax = plt.subplots()
	pvar = ax.plot(grid['lat'], var2d, color='black')
	# make_title_sim_time(ax, sim, model=model, timemean=timemean)
	ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
	ax.set_xticks(np.arange(-90,91,30))
	ax.set_xlabel('Latitude (deg)')
	ax.set_ylabel(ylabel)
	# ax.set_yticks(np.arange(-90,91,30))
	ax.xaxis.set_minor_locator(MultipleLocator(10))
	ax.yaxis.set_minor_locator(AutoMinorLocator())
	ax.set_xlim([-90,90])
	make_title_sim_time_seas(ax, sim, model=model, timemean=timemean, seasmean=seasmean, levstr='p=%g' % (plev))
	plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
	plt.show()
	
	# tikzplotlib.save('%s.tex' % (plotname))
