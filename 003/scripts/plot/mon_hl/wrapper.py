from r1_mon_hl import *
from deriv_mon_hl import *
from r1_mon_hl_scat import *
from pr_mon_hl import *
# from sens_mon_hl import *
# from dta_mon_hl import *
from flux_mon_hl import *
from rad_mon_hl import *
from co2_rad_mon_hl import *
from stf_mon_hl import *
from sic_mon_hl import *
from sic_stf_mon_hl import *
from eke_mon_hl import *
from dyn_mon_hl import *
from circ_mon_hl import *
from cld_mon_hl import *

# utype = 'rcp' # echam, ssp, rcp, ssp+hist, rcp+hist, hist
# utype = 'rcp+hist' # echam, ssp, rcp, ssp+hist, rcp+hist, hist
# utype = 'ssp+hist' # echam, ssp, rcp, ssp+hist, rcp+hist, hist
utype = 'echam' # echam, ssp, rcp, ssp+hist, rcp+hist, hist

mmm = 0
# try_load =0
# mmm = 1
try_load = 1
plotover = 'decomp'
# timemeans = ['djfmean', 'mammean', 'jjamean', 'sonmean']
timemeans = ['djfmean']
# latbnd = (60,90)
# latbnd = (70,90)
latbnd = (80,90)
# latbnd = (-90,90)
# latbnd = (75,85)
viewplt = 0
vertlev = 1
vertcoord = 'si'
legend =0
spread = 'std' # show spread as percentiles (prc) or stdev (std)?

if utype == 'echam':
    ####################################
    ## ECHAM
    ####################################
    sim = 'echam'

    ####################
    # yr_span = '0001_0009'

    # models = ['rp000190']
    # refclim = 'rp000134'

    ####################
    # yr_span = '0040_0339'

    # models = ['rp000191b']
    # refclim = 'rp000190b'

    # models = ['rp000188']
    # refclim = 'rp000134'

    ####################
    # yr_span = '0001_0079'

    # models = ['rp000130']
    # refclim = 'rp000130'

    ###################
    yr_span = '0040_0252'

    models = ['rp000191f']
    refclim = 'rp000190f'

    # models = ['rp000188']
    # refclim = 'rp000134'

    # models = ['rp000189']
    # refclim = 'rp000135'

    # models = ['rp000184']
    # refclim = 'rp000130'

    ###################
    # yr_span = '0020_0039'

    # # models = ['rp000190']
    # # refclim = 'rp000134'

    # # models = ['rp000130']
    # # refclim = 'rp000130'

    # models = ['rp000134']
    # refclim = 'rp000134'

    sim_ref = 'echam'
    timemean_ref = 'ymonmean-20'
    yr_span_ref = '0020_0039'

elif utype == 'rcpstd':
    ###################################
    # CMIP5
    ###################################
    # models = ['IPSL-CM5A-LR']
    sim = 'rcp85'
    sim_ref = 'historical'
    yr_span = '200601-229912'
    # sim = 'historical'
    # yr_span='186001-200512'
    refclim = 'hist-30'
    timemean_ref = 'ymonmean-30'
    yr_span_ref = '186001-200512'

    # STANDARD RCP8.5 RUNS, OUT TO 2100
    yr_span = '200601-210012'
    models = ['ACCESS1-0', 'ACCESS1-3', 'bcc-csm1-1', 'bcc-csm1-1-m', 'BNU-ESM', 'CanESM2', 'CCSM4', 'CESM1-BGC', 'CESM1-CAM5', 'CMCC-CESM', 'CMCC-CM', 'CMCC-CMS', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'FGOALS-g2', 'GFDL-CM3', 'GFDL-ESM2G', 'GFDL-ESM2M', 'GISS-E2-H', 'GISS-E2-H-CC', 'GISS-E2-R', 'GISS-E2-R-CC', 'HadGEM2-CC', 'HadGEM2-ES', 'inmcm4', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'IPSL-CM5B-LR', 'MIROC5', 'MIROC-ESM', 'MIROC-ESM-CHEM', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3', 'MRI-ESM1', 'NorESM1-M', 'NorESM1-ME']

elif utype == 'ssp':
    # EXTENDED SSP585 RUNS, TO 2300
    yr_span = '201501-230012'
    sim = 'ssp585'
    sim_ref = 'historical'
    # refclim = 'init'
    refclim = 'hist-30'
    timemean_ref = 'ymonmean-30'
    yr_span_ref = '186001-201412'
    models=['ACCESS-CM2', 'ACCESS-ESM1-5', 'CanESM5', 'IPSL-CM6A-LR', 'MRI-ESM2-0']

elif utype == 'rcp':
    # EXTENDED RCP8.5 RUNS, TO 2300
    yr_span = '200601-229912'
    sim = 'rcp85'
    sim_ref = 'historical'
    # refclim = 'init'
    refclim = 'hist-30'
    timemean_ref = 'ymonmean-30'
    yr_span_ref = '186001-200512'
    # # models = ['HadGEM2-ES', 'GISS-E2-H', 'GISS-E2-R', 'bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
    models = ['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'HadGEM2-ES', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
    # models = ['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'HadGEM2-ES', 'MPI-ESM-LR']
    # models = ['bcc-csm1-1', 'CCSM4', 'CSIRO-Mk3-6-0', 'HadGEM2-ES', 'MPI-ESM-LR', 'IPSL-CM5A-LR']
    # models = ['IPSL-CM5A-LR']

elif utype == 'ssp+hist':
    # MERGED HISTORICAL AND EXTENDED SSP585 RUNS, TO 2300
    yr_span = '186001-229912'
    sim = 'hist+ssp585'
    sim_ref = 'historical'
    refclim = 'hist-30'
    timemean_ref = 'ymonmean-30'
    yr_span_ref = '186001-201412'
    # refclim = 'init'
    models=['ACCESS-CM2', 'ACCESS-ESM1-5', 'CanESM5', 'CESM2-WACCM', 'IPSL-CM6A-LR', 'MRI-ESM2-0','MIROC-ES2L','GISS-E2-1-G','GISS-E2-1-H','UKESM1-0-LL']

elif utype == 'rcp+hist':
    # MERGED HISTORICAL AND EXTENDED RCP8.5 RUNS, TO 2300
    yr_span = '186001-229912'
    sim = 'hist+rcp85'
    sim_ref = 'historical'
    refclim = 'hist-30'
    timemean_ref = 'ymonmean-30'
    yr_span_ref = '186001-200512'
    # refclim = 'init'
    models = ['bcc-csm1-1', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'HadGEM2-ES', 'MPI-ESM-LR', 'IPSL-CM5A-LR','GISS-E2-H','GISS-E2-R']

####################################
## ERA
####################################
# models = [None]
# sim = 'era5'
# yr_span='1979_2019'
# refclim = 'init'

for timemean in timemeans:
    if mmm:
        r1_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # cld_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)
        
        # stf_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # sic_stf_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # deriv_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # r1_mon_hl_scat(sim, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # pr_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # flux_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # sens_mon_hl(sim, vertlev=vertlev, vertcoord=vertcoord, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # rad_mon_hl(sim, vertlev=vertlev, vertcoord=vertcoord, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # co2_rad_mon_hl(sim, vertlev=vertlev, vertcoord=vertcoord, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # dyn_mon_hl(sim, vertlev=vertlev, vertcoord=vertcoord, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

        # circ_mon_hl(sim, vertlev=vertlev, vertcoord=vertcoord, viewplt=viewplt, plotover=plotover, model=models, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, spread=spread)

    else:
        for model in models:
            print(model)

            r1_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # cld_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # stf_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # sic_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # sic_stf_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # deriv_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # pr_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # dta_mon_hl(sim, vertlev=vertlev, vertcoord=vertcoord, viewplt=viewplt, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # sens_mon_hl(sim, vertlev=vertlev, vertcoord=vertcoord, viewplt=viewplt, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # flux_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, yr_span_ref=yr_span_ref, timemean_ref=timemean_ref)

            # rad_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, sim_ref=sim_ref, timemean_ref=timemean_ref, yr_span_ref=yr_span_ref)

            # co2_rad_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend, sim_ref=sim_ref, timemean_ref=timemean_ref, yr_span_ref=yr_span_ref)

            # dyn_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # circ_mon_hl(sim, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

            # eke_mon_hl(sim, vertlev=vertlev, vertcoord=vertcoord, viewplt=viewplt, plotover=plotover, model=model, yr_span=yr_span, timemean=timemean, refclim=refclim, try_load=try_load, latbnd=latbnd, legend=legend)

