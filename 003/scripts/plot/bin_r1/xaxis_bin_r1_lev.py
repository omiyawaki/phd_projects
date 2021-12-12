def xaxis_bin_r1_lev(sim, **kwargs):

    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    latbnd = kwargs.get('latbnd', (80,90))
    plotvar = kwargs.get('plotvar', None)

    if plotvar == 'ga_dev':
        vmin = -50
        vmax = 150
        xlabel = '$(\Gamma_m - \Gamma)/\Gamma_m$ (%)'

    return vmin, vmax, xlabel
