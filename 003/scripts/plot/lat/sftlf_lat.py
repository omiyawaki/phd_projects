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
from netCDF4 import Dataset
import tikzplotlib

# sim = 'longrun'
# model = 'MPIESM12_abrupt16x'
# yr_span = '103'
# varname = 'sftlf'

sim = 'rcp85'
model = 'MPI-ESM-LR'
yr_span = '200601-230012'

# sim = 'echam'
# model = 'rp000140'
# yr_span = '0001_0039'
# varname = 'ahfs'

zonmean = '.zonmean' # zonal mean?
timemean = '' # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
categ = 'lat'

try_load = 0 # try to load data if available; otherwise, compute R1

# load data and plot directories
plotdir = get_plotdir(sim, model=model, yr_span=yr_span, categ=categ)

# location of pickled R1 data
sftlf_file = Dataset('/project2/tas1/ockham/data9/tas/CMIP5_RAW/MPI-ESM-LR/historical/atmos/fx/sftlf/r0i0p0/sftlf_fx_MPI-ESM-LR_historical_r0i0p0.nc', 'r')
sftlf = np.mean(sftlf_file.variables['sftlf'][:],1)

grid = {}
grid['lat'] = sftlf_file.variables['lat'][:]
grid['lon'] = sftlf_file.variables['lon'][:]

print(sftlf.shape)

plotname = '%s/sftlf_lat%s' % (plotdir, timemean)
fig, ax = plt.subplots()
lsft = ax.plot(grid['lat'], sftlf, color='black')
# make_title_sim_time(ax, sim, model=model, timemean=timemean)
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_xticks(np.arange(-90,91,30))
ax.set_xlabel('Latitude (deg)')
ax.set_ylabel('Land fraction (%)')
# ax.set_yticks(np.arange(-90,91,30))
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(MultipleLocator(10))
ax.set_xlim([-90,90])
ax.set_title('MPI-ESM-LR')
# plt.savefig('%s.png' % (plotname), dpi=300)
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)

# tikzplotlib.save('%s.tex' % (plotname))
