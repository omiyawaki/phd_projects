def yaxis_r1_mon_hl(sim, **kwargs):

    timemean = kwargs.get('timemean', '') # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
    latbnd = kwargs.get('latbnd', (80,90))
    plotover = kwargs.get('plotover', None) # plot overlay (sic for sea ice, ga_dev for lapse rate deviation, pr for precip, decomp for linear decomposition)?

    vmin = {}
    vmax = {}

    if sim == 'era5':
        if plotover == 'decomp':
            vmin['r1'] = 0.5
            vmax['r1'] = 1.2
        else:
            vmin['r1'] = 0.7
            vmax['r1'] = 1.05

        vmin['prfrac'] = 0.70
        vmax['prfrac'] = 1.05
    else:
        if latbnd[0] > 0: # NH
            if timemean in ['djfmean', 'sonmean', 'jjamean', 'mammean', 'allseas']: # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
                if plotover == 'decomp':
                    # vmin['r1'] = 0.65 # previously 0.5
                    # vmax['r1'] = 1.15 # previously 1.3
                    vmin['r1'] = 0.5 # previously 0.5
                    vmax['r1'] = 1.3 # previously 1.3
                else:
                    if sim == 'echam':
                        vmin['r1'] = 0.5
                        vmax['r1'] = 1.2
                    else:
                        # vmin['r1'] = 0.65
                        # vmax['r1'] = 1.1
                        vmin['r1'] = 0.35 # previously 0.5
                        vmax['r1'] = 1.25 # previously 1.3
                vmin['dr1dt'] = -0.3e-2
                vmax['dr1dt'] = 0.2e-2
                vmin['sic'] = -60
                vmax['sic'] = 120
                # vmin['ga_dev'] = -200
                # vmax['ga_dev'] = 250
                vmin['ga_dev'] = -300
                vmax['ga_dev'] = 300
                vmin['vhur'] = 50
                vmax['vhur'] = 80
                vmin['pr'] = -0.2
                vmax['pr'] = 2.2
                vmin['prc'] = -0.2
                vmax['prc'] = 1
                vmin['prl'] = 0.5
                vmax['prl'] = 1.8
                vmin['prfrac'] = -0.10
                vmax['prfrac'] = 0.50
                vmin['clt'] = 0
                vmax['clt'] = 100
                vmin['clwvi'] = -0.25e-1
                vmax['clwvi'] = 1.25e-1
                vmin['clivi'] = -0.25e-1
                vmax['clivi'] = 1.25e-1
            # elif timemean == 'jjamean':
            #     if plotover == 'decomp':
            #         vmin['r1'] = 0.5
            #         vmax['r1'] = 1.1
            #     else:
            #         vmin['r1'] = 0.8
            #         vmax['r1'] = 1.05
            #     vmin['sic'] = 0
            #     vmax['sic'] = 120
            #     vmin['pr'] = 0.6
            #     vmax['pr'] = 2.2
            #     vmin['prc'] = -0.3
            #     vmax['prc'] = 0.3
            #     vmin['prl'] = 0.5
            #     vmax['prl'] = 1.8
            elif timemean == 'yearmean' or timemean == '':
                vmin['r1'] = 0.7
                vmax['r1'] = 0.95
                vmin['sic'] = -40
                vmax['sic'] = 100
                vmin['ga_dev'] = -200
                vmax['ga_dev'] = 175
                vmin['pr'] = 0.6
                vmax['pr'] = 2.2
                vmin['prc'] = -0.1
                vmax['prc'] = 0.4
            elif timemean == 'ymonmean-30':
                # vmin['r1'] = 0.5
                # vmax['r1'] = 1.4
                vmin['r1'] = 0.7
                vmax['r1'] = 1.2
                vmin['sic'] = -40
                vmax['sic'] = 100
                vmin['ga_dev'] = -200
                vmax['ga_dev'] = 175
                vmin['pr'] = 0
                vmax['pr'] = 4
                vmin['prc'] = -0.1
                vmax['prc'] = 1
            elif timemean == 'ymonmean':
                vmin['r1'] = 0.5
                vmax['r1'] = 1.4
                vmin['sic'] = -40
                vmax['sic'] = 100
                vmin['ga_dev'] = -200
                vmax['ga_dev'] = 175
                vmin['pr'] = 0
                vmax['pr'] = 4
                vmin['prc'] = -0.5
                vmax['prc'] = 0.4
        else: # SH
            if timemean == 'djfmean': # type of time mean (yearmean, jjamean, djfmean, ymonmean-30)
                vmin['r1'] = 0.8 
                vmax['r1'] = 1.4
                vmin['sic'] = -40
                vmax['sic'] = 100
                vmin['ga_dev'] = 50
                vmax['ga_dev'] = 350
            elif timemean == 'jjamean':
                vmin['r1'] = 0.8 
                vmax['r1'] = 1.4
                vmin['sic'] = -40
                vmax['sic'] = 100
                vmin['ga_dev'] = 50
                vmax['ga_dev'] = 350
            elif timemean == 'yearmean' or timemean == '':
                vmin['r1'] = 0.8 
                vmax['r1'] = 1.4
                vmin['sic'] = -40
                vmax['sic'] = 100
                vmin['ga_dev'] = 50
                vmax['ga_dev'] = 350

        # vmin['prfrac'] = 100*vmin['r1']+10
        # vmax['prfrac'] = 100*vmax['r1']+10

    return vmin, vmax
