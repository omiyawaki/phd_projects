from ta_lev import *

####################################
## CMIP5
####################################

# # HISTORICAL
# sim = 'historical'
# yr_span = '186001-200512'
# # models = ['GISS-E2-H', 'GISS-E2-R', 'bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
# # models = ['HadGEM2-ES', 'bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
# models = ['bcc-csm1-1']

# EXTENDED RCP8.5 RUNS, TO 2300
sim = 'hist+rcp85'
yr_span = '186001-229912'
# models = ['HadGEM2-ES', 'bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
# models = ['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'HadGEM2-ES', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
models = ['bcc-csm1-1']

# # EXTENDED RCP8.5 RUNS, TO 2300
# sim = 'rcp85'
# yr_span = '200601-229912'
# # models = ['HadGEM2-ES', 'bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
# models = ['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'HadGEM2-ES', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
# # models = ['bcc-csm1-1']

timemean = 'ymonmean-30'
refclim = 'hist-30'


# mmm = 1
# try_load = 1
mmm = 0
try_load = 0
seas = ''
viewplt = 0
legend = 1

if mmm:
    ta_lev(sim, viewplt=viewplt, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, legend=legend, seas=seas)
else:
    for model in models:
        print(model)

        ta_lev(sim, viewplt=viewplt, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, legend=legend, seas=seas)
