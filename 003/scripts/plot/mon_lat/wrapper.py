from r1_mon_lat import *
from adv_mon_lat import *
from hydro_mon_lat import *
from gmse_mon_lat import *

# STANDARD RCP8.5 RUNS, OUT TO 2100
# sim = 'historical'
# refclim = 'hist-30'
# yr_span = '200601-210012'
# models = ['ACCESS1-0', 'ACCESS1-3', 'bcc-csm1-1', 'bcc-csm1-1-m', 'BNU-ESM', 'CanESM2', 'CCSM4', 'CESM1-BGC', 'CESM1-CAM5',  'CMCC-CESM', 'CMCC-CM', 'CMCC-CMS', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'FGOALS-g2', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M', 'GISS-E2-H', 'GISS-E2-H-CC', 'GISS-E2-R', 'GISS-E2-R-CC', 'HadGEM2-CC', 'HadGEM2-ES', 'inmcm4', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR', 'MIROC5', 'MIROC-ESM', 'MIROC-ESM-CHEM', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3', 'MRI-ESM1', 'NorESM1-M', 'NorESM1-ME']

# # HIST+EXTENDED RCP8.5 RUNS, TO 2300
# sim = 'hist+rcp85'
# timemean = 'djfmean'
# refclim = 'hist-30'
# yr_span = '186001-229912'
# # models = ['GISS-E2-H', 'GISS-E2-R', 'bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
# models = ['HadGEM2-ES', 'bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
# # models = ['bcc-csm1-1']

# # EXTENDED RCP8.5 RUNS, TO 2300
# sim = 'rcp85'
# timemean = 'yearmean'
# refclim = 'hist-30'
# yr_span = '200601-229912'
# # models = ['GISS-E2-H', 'GISS-E2-R', 'bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
# models = ['HadGEM2-ES', 'bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']

# ECHAM
sim = 'echam'
timemean = 'djfmean'
yr_span = '0040_0239'

# models = ['rp000188']
# refclim = 'rp000134'

models = ['rp000191b']
refclim = 'rp000190b'

####################################
## ERA
####################################
# models = [None]
# sim = 'era5'
# yr_span='1979_2019'
# refclim = 'init'
# timemean = 'ymonmean'

mmm = 0
try_load = 1
# try_load = 0
viewplt = 0
zonmean = 0
legend = 1

if mmm:
    r1_mon_lat(sim, viewplt=viewplt, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load)

else:
    for model in models:
            
        print(model)

        r1_mon_lat(sim, viewplt=viewplt, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load)
        
        # adv_mon_lat(sim, viewplt=viewplt, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load)

        # hydro_mon_lat(sim, viewplt=viewplt, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load)

        # gmse_mon_lat(sim, zonmean=zonmean, viewplt=viewplt, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, legend=legend)
