import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
from misc.dirnames import get_datadir, get_plotdir
from misc.filenames import filenames_raw
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
# varname = 'hfls'

# sim = 'rcp85'
# model = 'MPI-ESM-LR'
# yr_span = '200601-230012'
# varname = 'hfls'

sim = 'echam'
model = 'rp000140'
yr_span = '0001_0039'
varname = 'ahfl'

zonmean = '.zonmean' # zonal mean?
timemean = '.djfmean' # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
categ = 'mon_lat'

try_load = 0 # try to load data if available; otherwise, compute R1

# load data and plot directories
plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

# location of pickled R1 data
hfls_file = filenames_raw(sim, varname, model=model, timemean=timemean, yr_span=yr_span)
if not zonmean:
    hfls = hfls_file.variables[varname][:]
else:
    hfls = np.mean(hfls_file.variables[varname][:], 2)

if sim == 'echam':
    hfls = -hfls

grid = {}
grid['lat'] = hfls_file.variables['lat'][:]
grid['lon'] = hfls_file.variables['lon'][:]

print(hfls.shape)

rolling_mean = 0; # smooth data using a rolling mean? (units: yr)
hfls_filt = uniform_filter(hfls, [rolling_mean,0]) # apply rolling mean

[mesh_lat, mesh_time] = np.meshgrid(grid['lat'], np.arange(hfls.shape[0])) # create mesh

plotname = '%s/hfls_mon_lat%s' % (plotdir, timemean)
fig, ax = plt.subplots()
vmin = -10
vmax = 200
csf = ax.contourf(mesh_time, mesh_lat, hfls_filt, np.arange(vmin,vmax,10), cmap='Blues', vmin=vmin, vmax=vmax)
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
cbar.set_label('LH (Wm$^{-2}$)')
# plt.savefig('%s.png' % (plotname), dpi=300)
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)

# tikzplotlib.save('%s.tex' % (plotname))