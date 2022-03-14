from circ_lat_lev import *

# models = ['MPI-ESM-LR']
# yr_span = '186001-200512'

# STANDARD RCP8.5 RUNS, OUT TO 2100
# yr_span = '200601-210012'
# models = ['ACCESS1-0', 'ACCESS1-3', 'bcc-csm1-1', 'bcc-csm1-1-m', 'BNU-ESM', 'CanESM2', 'CCSM4', 'CESM1-BGC', 'CESM1-CAM5',  'CMCC-CESM', 'CMCC-CM', 'CMCC-CMS', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'FGOALS-g2', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M', 'GISS-E2-H', 'GISS-E2-H-CC', 'GISS-E2-R', 'GISS-E2-R-CC', 'HadGEM2-CC', 'HadGEM2-ES', 'inmcm4', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR', 'MIROC5', 'MIROC-ESM', 'MIROC-ESM-CHEM', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3', 'MRI-ESM1', 'NorESM1-M', 'NorESM1-ME']

# HISTORICAL
sim = 'historical'
yr_span = '186001-200512'
# models = ['GISS-E2-H', 'GISS-E2-R', 'bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
models = ['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'HadGEM2-ES', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
# models = ['bcc-csm1-1']

# # EXTENDED RCP8.5 RUNS, TO 2300
# sim = 'rcp85'
# yr_span = '200601-229912'
# models = ['HadGEM2-ES', 'bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
# # models = ['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']

timemean = 'ymonmean-30'
refclim = 'hist-30'

#####################################
### ERA
#####################################
#models = [None]
#sim = 'era5'
#yr_span='1979_2019'
#refclim = 'init'
#timemean = 'ymonmean'

mmm = 0
# try_load = 0
# mmm = 1
try_load = 1
viewplt = 0
seas = 'djf' # choose from ann, djf, mam, jja, son

if mmm:
    circ_lat_lev(sim, viewplt=viewplt, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, seas=seas)

else:
    for model in models:
            
        print(model)

        circ_lat_lev(sim, viewplt=viewplt, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, seas=seas)
        
