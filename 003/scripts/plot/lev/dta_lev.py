import os
import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/scripts')
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import interp1d, interp2d
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

categ = 'lev'

sim = 'hist+rcp85'
time = '186001-229912'
model = 'mmm'
models = ['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MPI-ESM-LR']
plotdir = '/project2/tas1/miyawaki/projects/003/plot/%s/%s/%s/%s' % (sim, model, time, categ)

if not os.path.exists(plotdir):
    os.mkdir(plotdir)

##################################
# LOAD DATA
##################################
if model == 'mmm':
    plev0 = 1e2*np.linspace(100,1000,101)
    lat0 = np.linspace(80,90,11)
    clat0 = np.cos(np.radians(lat0))
    ntime = 442

ta0 = np.empty([len(models), ntime, len(plev0)])
for im in range(len(models)):
    m = models[im]
    datadir = '/project2/tas1/miyawaki/projects/003/data/raw/%s/%s' % (sim, m)
    raw = Dataset('%s/ta_sm_Amon_%s_%s_r1i1p1_186001-229912.zonmean.djfmean.nc' % (datadir, m, sim))
    ta = np.squeeze(raw.variables['ta_sm'][:])
    lat = raw.variables['lat'][:]
    plev = raw.variables['plev'][:]

    # average over arctic
    flat = interp1d(lat, ta, axis=2, bounds_error=False, fill_value='extrapolate')
    ta_cap = flat(lat0)
    ta_acap = np.nansum(ta_cap*clat0[None,None,:], axis=2) / np.nansum(clat0) 

    # interpolate to standard plev
    flev = interp1d(plev, ta_acap, axis=1, bounds_error=False, fill_value='extrapolate')
    ta0[im,:,:] = flev(plev0)

tammm = np.nanmean(ta0, axis=0)

############################################
# PLOT (T)
############################################
plotname = '%s/ta' % (plotdir)
fig, ax = plt.subplots(figsize=(4,3))
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(tammm[145+50,:], 1e-2*plev0, '-k', label='2050')
ax.plot(tammm[145+100,:], 1e-2*plev0, '-r', label='2100')
ax.plot(tammm[145+150,:], 1e-2*plev0, '-b', label='2150')
ax.plot(tammm[145+200,:], 1e-2*plev0, '-g', label='2200')
ax.plot(tammm[145+250,:], 1e-2*plev0, '-m', label='2250')
# make_title_sim_time_seas(ax, sim, model=model, timemean=timemean, seasmean=seas)
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_xlabel('T (K)')
ax.set_ylabel('$p$ (hPa)')
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylim(100,1000)
ax.invert_yaxis()
ax.legend()
plt.tight_layout()
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
plt.close()

############################################
# PLOT (DELTA T)
############################################
plotname = '%s/dta' % (plotdir)
fig, ax = plt.subplots(figsize=(4,3))
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(tammm[145+50,:]-tammm[145,:], 1e-2*plev0, '-k', label='2050')
ax.plot(tammm[145+100,:]-tammm[145,:], 1e-2*plev0, '-r', label='2100')
ax.plot(tammm[145+150,:]-tammm[145,:], 1e-2*plev0, '-b', label='2150')
ax.plot(tammm[145+200,:]-tammm[145,:], 1e-2*plev0, '-g', label='2200')
ax.plot(tammm[145+250,:]-tammm[145,:], 1e-2*plev0, '-m', label='2250')
# make_title_sim_time_seas(ax, sim, model=model, timemean=timemean, seasmean=seas)
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_xlabel('$\Delta T$ (K)')
ax.set_ylabel('$p$ (hPa)')
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylim(100,1000)
ax.invert_yaxis()
ax.legend()
plt.tight_layout()
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
plt.close()

############################################
# PLOT (DELTA T SEQ)
############################################
plotname = '%s/dta_seq' % (plotdir)
fig, ax = plt.subplots(figsize=(4,3))
ax.axhline(0, color='k', linewidth=0.5)
ax.plot(tammm[140+50,:]-tammm[140,:], 1e-2*plev0, '-k', label='2050-2000')
ax.plot(tammm[140+100,:]-tammm[140+50,:], 1e-2*plev0, '-r', label='2100-2050')
ax.plot(tammm[140+150,:]-tammm[140+100,:], 1e-2*plev0, '-b', label='2150-2100')
ax.plot(tammm[140+200,:]-tammm[140+150,:], 1e-2*plev0, '-g', label='2200-2150')
ax.plot(tammm[140+250,:]-tammm[140+200,:], 1e-2*plev0, '-m', label='2250-2200')
# make_title_sim_time_seas(ax, sim, model=model, timemean=timemean, seasmean=seas)
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.set_xlabel('$\Delta T$ (K)')
ax.set_ylabel('$p$ (hPa)')
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_ylim(100,1000)
ax.invert_yaxis()
ax.legend()
plt.tight_layout()
plt.savefig('%s.pdf' % (plotname), format='pdf', dpi=300)
plt.close()

