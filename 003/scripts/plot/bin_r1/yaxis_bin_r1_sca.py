def yaxis_bin_r1_sca(sim, **kwargs):

    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    latbnd = kwargs.get('latbnd', (80,90))
    plotvar = kwargs.get('plotvar', None)

    if plotvar == 'pr':
        vmin = -0.2
        vmax = 6.6
        ylabel = '$P$ (mm d$^{-1}$)'
    if plotvar == 'prc':
        vmin = -0.2
        vmax = 6.6
        ylabel = '$P_c$ (mm d$^{-1}$)'
    elif plotvar == 'prfrac':
        vmin = -0.2
        vmax = 1.2
        ylabel = '$P_c / P$ (unitless)'

    return vmin, vmax, ylabel
