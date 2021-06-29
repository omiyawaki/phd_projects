import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import get_datadir, get_plotdir
from misc.filenames import filenames_raw
from misc.translate import *
from plot.titles import make_title_sim_time
import os
import pickle
import numpy as np
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from netCDF4 import Dataset
import tikzplotlib

def plot_raw2d_mon_lat(sim, varname, **kwargs):
	
	categ = 'lat'

	timemean = kwargs.get('timemean', 'yearmean') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
	domain = kwargs.get('domain', '')
	try_load = kwargs.get('try_load', 1) # try to load data if available; otherwise, compute R1

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
		
	if translate_varname(varname) == 'tas':
		vmin = 200
		vmax = 310
		vint = 2
		ylabel = '$T_{2\,m}$ (K)'
	elif translate_varname(varname) == 'pr':
		vmin = 0
		vmax = 20
		vint = 1
		ylabel = '$P$ (mm d$^{-1}$)'
	
	# load data and plot directories
	plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)
	
	# location of raw 2d data
	file = filenames_raw(sim, varname, model=model, timemean=timemean, yr_span=yr_span)
	
	# take zonal mean
	var2d = np.mean(np.squeeze(file.variables[varname][:]),2)
	print(var2d.shape)
	
	# if precip, convert to mm/d
	if translate_varname(varname) == 'pr':
		var2d = 86400 * var2d
		
	grid = {}
	grid['lat'] = file.variables['lat'][:]
	
	[mesh_lat, mesh_time] = np.meshgrid(grid['lat'], yr_base + np.arange(var2d.shape[0])) # create mesh

	##################################
	# PLOT
	##################################
	plotname = '%s/%s_mon_lat.%s' % (plotdir, varname, timemean)
	fig, ax = plt.subplots()
	csf = ax.contourf(mesh_time, mesh_lat, var2d, np.arange(vmin,vmax,vint), cmap='viridis', vmin=vmin, vmax=vmax, extend='both')
	make_title_sim_time(ax, sim, model=model, timemean=timemean)
	ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
	if 'ymonmean' in timemean:
		ax.set_xticks(np.arange(0,12,1))
		ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
	else:
		ax.set_xlabel('Year')
	ax.set_ylabel('Latitude (deg)')
	ax.set_yticks(np.arange(-90,91,30))
	ax.xaxis.set_minor_locator(AutoMinorLocator())
	ax.yaxis.set_minor_locator(MultipleLocator(10))
	cbar = plt.colorbar(csf)
	cbar.set_label(ylabel)
	# plt.savefig('%s.png' % (plotname), dpi=300)
	plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
	plt.show()
	
	# tikzplotlib.save('%s.tex' % (plotname))
