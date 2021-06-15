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
# varname = 'hfss'

# sim = 'rcp85'
# model = 'MPI-ESM-LR'
# yr_span = '200601-230012'
# varname = 'hfss'

sim = 'echam'
model = 'rp000140'
yr_span = '0001_0039'
varname = 'ahfs'

zonmean = '.zonmean' # zonal mean?
timemean = '.djfmean' # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
categ = 'mon_lat'

try_load = 0 # try to load data if available; otherwise, compute R1

# load data and plot directories
plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

# location of pickled R1 data
hfss_file = filenames_raw(sim, varname, model=model, timemean=timemean, yr_span=yr_span)
if not zonmean:
    hfss = hfss_file.variables[varname][:]
else:
    hfss = np.mean(hfss_file.variables[varname][:], 2)

if sim == 'echam':
    hfss = -hfss

grid = {}
grid['lat'] = hfss_file.variables['lat'][:]
grid['lon'] = hfss_file.variables['lon'][:]

print(hfss.shape)

rolling_mean = 0; # smooth data using a rolling mean? (units: yr)
hfss_filt = uniform_filter(hfss, [rolling_mean,0]) # apply rolling mean

[mesh_lat, mesh_time] = np.meshgrid(grid['lat'], np.arange(hfss.shape[0])) # create mesh

plotname = '%s/hfss_mon_lat%s' % (plotdir, timemean)
fig, ax = plt.subplots()
vmin = -70
vmax = 70
csf = ax.contourf(mesh_time, mesh_lat, hfss_filt, np.arange(vmin,vmax,10), cmap='RdBu_r', vmin=vmin, vmax=vmax, extend='both')
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
cbar.set_label('SH (Wm$^{-2}$)')
# plt.savefig('%s.png' % (plotname), dpi=300)
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)

# tikzplotlib.save('%s.tex' % (plotname))