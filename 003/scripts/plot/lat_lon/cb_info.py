import numpy as np

def cb_info(plotvar, **kwargs):

    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    latbnd = kwargs.get('latbnd', (80,90))

    if plotvar == 'prfrac':
        clevs = np.arange(0,1.1,0.1)
        clabel = '$P_c/P$ (unitless)'

    elif plotvar == 'prc':
        clevs = np.arange(0,1.1,0.1)
        clabel = '$P_c$ (mm d$^{-1}$)'

    elif plotvar == 'sic':
        clevs = np.arange(0,110,10)
        clabel = 'Sea ice fraction (%)'

    return clevs, clabel
