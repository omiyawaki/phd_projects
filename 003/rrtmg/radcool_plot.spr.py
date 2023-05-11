import os
import sys
import pickle
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from tools import get_modelidx, get_modelstr

models=['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'IPSL-CM5A-LR', 'HadGEM2-ES', 'MPI-ESM-LR']

nlev0=101
plev0=np.logspace(3,5,nlev0)

nm = len(models)
nt = 147

rc=np.empty([nm,nt,nlev0])

plotdir = '%s/plot/mmm' % ( os.getcwd()) # plot directory
if not os.path.exists(plotdir):
    os.mkdir(plotdir)

c=0 # counter
for model in models:

    # specify directory for output data and plots
    datadir = '%s/data' % ( os.getcwd() ) # output data directory

    ##################################################
    # LOAD GCM DATA
    ##################################################
    ds = xr.open_dataset('/project2/tas1/miyawaki/projects/003/echam/ghg/ghg_rcp85_1765-2500_c100203.nc')
    co2 = 1e-6*ds.CO2.sel(time=slice('2006-01-01', '2300-12-31')).data

    [mmm_clim, grid, indiv_clim] = pickle.load(open('../climlab/input_data/clima.pickle', 'rb'))
    [mmm_forc, grid, indiv_forc] = pickle.load(open('../climlab/input_data/forcing.pickle', 'rb'))

    if model == 'mmm':
        gcm_clim = mmm_clim
        gcm_forc = mmm_forc
    else:
        imod = get_modelidx(model)
        gcm_clim = {}
        gcm_forc = {}
        for varname in indiv_clim:
            gcm_clim[varname] = indiv_clim[varname][imod,...]
            gcm_forc[varname] = indiv_forc[varname][imod,...]

    ##################################################
    # LOAD RRTMG DATA
    ##################################################
    # RRTMG climatology
    ref = {}
    [ref['ra'], ref['racs'], ref['diag']] = pickle.load(open('%s/rad.hist.%s.pickle' % (datadir, model), 'rb'))

    # RRTMG full RCP8.5 response
    full = {}
    [full['ra'], full['racs'], full['diag']] = pickle.load(open('%s/rad.rcp.%s.pickle' % (datadir, model), 'rb'))

    ##################################################
    # MISC data
    ##################################################
    plev = 1e-2*grid['lev']
    modelstr = get_modelstr(model)

    ##################################################
    # HISTORICAL
    ##################################################
    time0 = np.arange(1860,2007)

    ##################################################
    # COMPUTE RADCOOL
    ##################################################
    nlev=len(plev)
    nt_ref=len(ref['diag'])
    for i in range(nt_ref):
        rc[c,i,:]=ref['diag'][i].TdotSW + ref['diag'][i].TdotLW

    ##################################################
    # END OF LOOP, UPDATE COUNTER
    ##################################################
    c=c+1
    
##################################################
# RADCOOL
##################################################

# take time mean
rc=np.mean(rc, axis=1)

# take multimodel statistics
mean=np.mean(rc, axis=0)
print(mean.shape)
p25=np.percentile(rc, 25, axis=0)
p75=np.percentile(rc, 75, axis=0)

# plot radiative cooling rate
plotname='%s/radcool.pdf' % (plotdir)
fig, ax = plt.subplots()
ax.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax.plot(mean, plev, 'k')
ax.fill_betweenx(plev, p25, p75, plev, alpha=0.2, color='k', edgecolor=None)
# ax.set_xlim([time0[0], time0[-1]])
ax.set_ylim([100,1000])
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('Radiative Cooling (K d$^{-1}$)')
ax.set_ylabel(r'$p$ (hPa)')
ax.set_title('Historical CMIP5 mean, $80-90^\circ$N')
plt.gca().invert_yaxis()
fig.set_size_inches(4,3.3)
plt.tight_layout()
plt.savefig(plotname, format='pdf', dpi=300)
plt.close()


