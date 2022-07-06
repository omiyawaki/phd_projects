import os
import sys
import pickle
import climlab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# specify directory for output data
datadir = '%s/data' % ( os.getcwd() ) # output data directory

lat=83.3
day=0
# co2_lev=280e-6 # volume mixing ratio
co2_lev=348e-6 # volume mixing ratio
# co2_lev=4*280e-6 # volume mixing ratio
n_yr = 1 # number of years to run

# specify plot directory
plotdir = '%s/plots/%g' % ( os.getcwd(), 1e6*co2_lev ) # plot directory
if not os.path.exists(plotdir):
    os.mkdir(plotdir)

par={}
par['n']=2
par['b']=1
par['g']=9.81 # m s**-2
par['cp']=1004 # J kg**-1 K**-1
par['ps']=1000 # hPa

# load simulation data
[rcm, Ra, Rsfc, Fadv_lev, gcm_clim] = pickle.load(open('%s/rcm.%g.%g.pickle' % (datadir, 1e6*co2_lev, n_yr), 'rb'))

# PLOT prescribed advective heat flux vertical structure
fig, ax = plt.subplots()
ax.plot(Fadv_lev, rcm.lev, '-k')
fig.set_size_inches(4,3)
ax.set_ylim([1e2,1e3])
ax.invert_yaxis()
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel(r'$\partial_t T_{\mathrm{adv}}$ (K s$^{-1}$)')
ax.set_ylabel('p (hPa)')
plt.tight_layout()
plt.savefig('%s/Fadv_lev.pdf' % (plotdir), format='pdf', dpi=300)

# PLOT T
fig, ax = plt.subplots()
ax.plot(gcm_clim['ta'], gcm_clim['plev'], '-r')
ax.plot(rcm.Ts, 1e3, '.k')
ax.plot(rcm.Tatm, rcm.lev, '-k')
fig.set_size_inches(4,3)
ax.set_ylim([1e2,1e3])
ax.invert_yaxis()
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('T (K)')
ax.set_ylabel('p (hPa)')
plt.tight_layout()
plt.savefig('%s/tprofile.pdf' % (plotdir), format='pdf', dpi=300)

# PLOT T (ZOOM IN)
fig, ax = plt.subplots()
ax.plot(gcm_clim['ta'], gcm_clim['plev'], '-r')
ax.plot(rcm.Ts, 1e3, '.k')
ax.plot(rcm.Tatm, rcm.lev, '-k')
fig.set_size_inches(4,3)
ax.set_ylim([9e2,1e3])
ax.invert_yaxis()
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('T (K)')
ax.set_ylabel('p (hPa)')
plt.tight_layout()
plt.savefig('%s/tprofile_zoom.pdf' % (plotdir), format='pdf', dpi=300)
