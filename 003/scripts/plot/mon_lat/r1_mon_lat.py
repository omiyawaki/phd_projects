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
# yr_span = '999'

sim = 'rcp85'
model = 'MPI-ESM-LR'
yr_span = '200601-230012'

# sim = 'echam'
# model = 'rp000140'
# yr_span = '0001_0039'

zonmean = '.zonmean' # zonal mean?
timemean = '' # type of time mean (.yearmean, .jjamean, .djfmean, .ymonmean-30)
categ = 'mon_lat'

try_load = 1 # try to load data if available; otherwise, compute R1

# load data and plot directories
datadir = get_datadir(sim, model=model, yr_span=yr_span)
plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

# location of pickled R1 data
r1_file = '%s/r1%s%s.pickle' % (datadir, zonmean, timemean)

if not (os.path.isfile(r1_file) and try_load):
    save_r1(sim, model=model, zonmean=zonmean, timemean=timemean, yr_span=yr_span)

[r1, grid] = pickle.load(open(r1_file, 'rb'))

# print(np.reshape(r1, (-1,96,12)).shape)
if timemean == '':
    r1 = np.mean(np.reshape(r1, (-1,12,r1.shape[1])),1)

rolling_mean = 0; # smooth data using a rolling mean? (units: yr)
r1_filt = uniform_filter(r1, [rolling_mean,0]) # apply rolling mean

[mesh_lat, mesh_time] = np.meshgrid(grid['lat'], np.arange(r1.shape[0])) # create mesh

##################################
# REGULAR
##################################
plotname = '%s/r1_mon_lat%s' % (plotdir, timemean)
fig, ax = plt.subplots()
vmin = -1.7
vmax = 1.7
csf = ax.contourf(mesh_time, mesh_lat, r1_filt, np.arange(vmin,vmax,0.1), cmap='RdBu', vmin=vmin, vmax=vmax, extend='both')
cs_rae = ax.contour(mesh_time, mesh_lat, r1_filt, levels=[0.9], colors='royalblue', linewidths=3)
cs_rce = ax.contour(mesh_time, mesh_lat, r1_filt, levels=[0.1], colors='sandybrown', linewidths=3)
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
cbar.set_label('$R_1$ (unitless)')
# plt.savefig('%s.png' % (plotname), dpi=300)
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)

##################################
# NH RCE 
##################################
plotname = '%s/rce_mon_lat%s' % (plotdir, timemean)
fig, ax = plt.subplots()
cs_rce = ax.contour(mesh_time, mesh_lat, r1_filt, levels=[0.1], colors='sandybrown', linewidths=1)
make_title_sim_time(ax, sim, model=model, timemean=timemean)
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
if 'ymonmean' in timemean:
    ax.set_xticks(np.arange(0,12,1))
    ax.set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
else:
    ax.set_xlabel('Year')
ax.set_ylabel('Latitude (deg)')
ax.set_ylim([40,45])
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(MultipleLocator(10))
# plt.savefig('%s.png' % (plotname), dpi=300)
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)