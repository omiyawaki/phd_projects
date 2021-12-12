from r1_mon_ll import *

####################################
## CMIP5
####################################
# models = ['MPI-ESM-LR']
# sim = 'rcp85'
# yr_span = '200601-230012'
# # sim = 'historical'
# # yr_span='186001-200512'
# refclim = 'hist-30'

# STANDARD RCP8.5 RUNS, OUT TO 2100
# yr_span = '200601-210012'
# models = ['ACCESS1-0', 'ACCESS1-3', 'bcc-csm1-1', 'bcc-csm1-1-m', 'BNU-ESM', 'CanESM2', 'CCSM4', 'CESM1-BGC', 'CESM1-CAM5', 'CMCC-CESM', 'CMCC-CM', 'CMCC-CMS', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'FGOALS-g2', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M', 'GISS-E2-H', 'GISS-E2-H-CC', 'GISS-E2-R', 'GISS-E2-R-CC', 'HadGEM2-CC', 'HadGEM2-ES', 'inmcm4', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR', 'MIROC5', 'MIROC-ESM', 'MIROC-ESM-CHEM', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3', 'MRI-ESM1', 'NorESM1-M', 'NorESM1-ME']

# EXTENDED RCP8.5 RUNS, TO 2300
# yr_span = '200601-229912'
# sim = 'rcp85'
# refclim = 'hist-30'
# models = ['HadGEM2-ES', 'GISS-E2-H', 'GISS-E2-R', 'bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
# models = ['HadGEM2-ES', 'bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
# models = ['MPI-ESM-LR']

####################################
## ERA
####################################
models = [None]
sim = 'era5'
yr_span='1979_2019'
refclim = 'init'

mmm = 0
plotover = 'pr'
# timemeans = ['djfmean', 'mammean', 'jjamean', 'sonmean']
timemeans = ['ymonmean']
viewplt = 0
try_load = 0
vertlev = 1
vertcoord = 'si'
legend = 1

latbnd_list = [
        (-20,-5),
        (-20,-5),
        (-20,-5),
        (5,15),
        (5,20),
        (10,30)
        ]

lonbnd_list = [
        (10,50),
        (120,150),
        (285,315),
        (0,30),
        (250,270),
        (70,100)
        ]

for region in range(len(latbnd_list)):

    latbnd = latbnd_list[region]
    lonbnd = lonbnd_list[region]

    for timemean in timemeans:
        if mmm:
            r1_mon_ll(sim, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, lonbnd=lonbnd, legend=legend)

        else:
            for model in models:
                print(model)

                r1_mon_ll(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, lonbnd=lonbnd, legend=legend)

