from r1_mon_hl import *
from deriv_mon_hl import *
from r1_mon_hl_scat import *
from pr_mon_hl import *
from sens_mon_hl import *
from dta_mon_hl import *
from flux_mon_hl import *
from rad_mon_hl import *
from eke_mon_hl import *
from dyn_mon_hl import *
from circ_mon_hl import *

####################################
## CMIP5
####################################
# models = ['IPSL-CM5A-LR']
# sim = 'rcp85'
# yr_span = '200601-229912'
# # sim = 'historical'
# # yr_span='186001-200512'
# refclim = 'hist-30'

# STANDARD RCP8.5 RUNS, OUT TO 2100
# yr_span = '200601-210012'
# models = ['ACCESS1-0', 'ACCESS1-3', 'bcc-csm1-1', 'bcc-csm1-1-m', 'BNU-ESM', 'CanESM2', 'CCSM4', 'CESM1-BGC', 'CESM1-CAM5', 'CMCC-CESM', 'CMCC-CM', 'CMCC-CMS', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'FGOALS-g2', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M', 'GISS-E2-H', 'GISS-E2-H-CC', 'GISS-E2-R', 'GISS-E2-R-CC', 'HadGEM2-CC', 'HadGEM2-ES', 'inmcm4', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR', 'MIROC5', 'MIROC-ESM', 'MIROC-ESM-CHEM', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3', 'MRI-ESM1', 'NorESM1-M', 'NorESM1-ME']

# EXTENDED RCP8.5 RUNS, TO 2300
yr_span = '200601-229912'
sim = 'rcp85'
# refclim = 'init'
refclim = 'hist-30'
# # models = ['HadGEM2-ES', 'GISS-E2-H', 'GISS-E2-R', 'bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
# models = ['HadGEM2-ES', 'bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
models = ['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'HadGEM2-ES', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
# models = ['bcc-csm1-1']

# # MERGED HISTORICAL AND EXTENDED RCP8.5 RUNS, TO 2300
# yr_span = '186001-229912'
# sim = 'hist+rcp85'
# refclim = 'hist-30'
# # refclim = 'init'
# # models = ['HadGEM2-ES', 'GISS-E2-H', 'GISS-E2-R', 'bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
# # models = ['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
# models = ['bcc-csm1-1', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
# # models = ['bcc-csm1-1']

####################################
## ERA
####################################
# models = [None]
# sim = 'era5'
# yr_span='1979_2019'
# refclim = 'init'

# mmm = 0
try_load =0
mmm = 1
# try_load = 1
plotover = 'prfrac'
# timemeans = ['djfmean', 'mammean', 'jjamean', 'sonmean']
timemeans = ['djfmean']
# latbnd = (60,90)
# latbnd = (70,90)
latbnd = (80,90)
# latbnd = (75,85)
viewplt = 0
vertlev = 1
vertcoord = 'si'
legend = 1
spread = 'std' # show spread as percentiles (prc) or stdev (std)?

for timemean in timemeans:
    if mmm:
        r1_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # deriv_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # r1_mon_hl_scat(sim, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # pr_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # flux_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # sens_mon_hl(sim, vertlev=vertlev, vertcoord=vertcoord, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # rad_mon_hl(sim, vertlev=vertlev, vertcoord=vertcoord, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # dyn_mon_hl(sim, vertlev=vertlev, vertcoord=vertcoord, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # circ_mon_hl(sim, vertlev=vertlev, vertcoord=vertcoord, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

    else:
        for model in models:
            print(model)

            # r1_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # deriv_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # pr_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # dta_mon_hl(sim, vertlev=vertlev, vertcoord=vertcoord, viewplt=viewplt, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # sens_mon_hl(sim, vertlev=vertlev, vertcoord=vertcoord, viewplt=viewplt, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # flux_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            rad_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # dyn_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # circ_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # eke_mon_hl(sim, vertlev=vertlev, vertcoord=vertcoord, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

