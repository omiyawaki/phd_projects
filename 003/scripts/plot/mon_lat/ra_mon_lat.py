import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import get_datadir, get_plotdir
from proc.r1 import save_r1
from plot.titles import make_title_sim_time
import os
import pickle
import numpy as np
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import tikzplotlib

# sim = 'longrun'
# model = 'MPIESM12_abrupt16x'
# yr_span = '103'

# sim = 'rcp85'
# model = 'MPI-ESM-LR'
# yr_span = '200601-230012'

sim = 'echam'
model = 'rp000140'
yr_span = '0001_0039'

zonmean = '.zonmean' # zonal mean?
timemean = '.yearmean' # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
categ = 'mon_lat'

try_load = 0 # try to load data if available; otherwise, compute R1

# load data and plot directories
datadir = get_datadir(sim, model=model, yr_span=yr_span)
plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

# location of pickled R1 data
ra_file = '%s/ra%s%s.pickle' % (datadir, zonmean, timemean)

if not (os.path.isfile(ra_file) and try_load):
    save_r1(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

[ra, grid] = pickle.load(open(ra_file, 'rb'))

rolling_mean = 0; # smooth data using a rolling mean? (units: yr)
ra_filt = uniform_filter(ra, [rolling_mean,0]) # apply rolling mean

[mesh_lat, mesh_time] = np.meshgrid(grid['lat'], np.arange(ra.shape[0])) # create mesh

plotname = '%s/ra_mon_lat%s' % (plotdir, timemean)
fig, ax = plt.subplots()
vmin = -150
vmax = 150 
csf = ax.contourf(mesh_time, mesh_lat, ra_filt, np.arange(vmin,vmax,10), cmap='RdBu_r', vmin=vmin, vmax=vmax, extend='both')
cs_0 = ax.contour(mesh_time, mesh_lat, ra_filt, levels=[0], colors='gray', linewidths=3)
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
cbar.set_label('$R_a$ (Wm$^{-2}$)')
# plt.savefig('%s.png' % (plotname), dpi=300)
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)

# tikzplotlib.save('%s.tex' % (plotname))