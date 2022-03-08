def make_title_sim(ax, sim, **kwargs):

    model = kwargs.get('model')

    ax.set_title('%s' % (base_title_str(sim, model=model)))

def make_title_sim_time(ax, sim, **kwargs):

    model = kwargs.get('model')
    timemean = kwargs.get('timemean')

    ax.set_title('%s, %s' % (base_title_str(sim, model=model), time_str(timemean)))

def make_title_sim_time_seas(ax, sim, **kwargs):

    model = kwargs.get('model')
    timemean = kwargs.get('timemean')
    seasmean = kwargs.get('seasmean')
    levstr = kwargs.get('levstr', None)

    if levstr is None:
        ax.set_title('%s, %s, %s' % (base_title_str(sim, model=model), time_str(timemean), seas_str(seasmean)))
    else:
        ax.set_title('%s, %s, %s, %s' % (base_title_str(sim, model=model), time_str(timemean), seas_str(seasmean), levstr))
    
def make_title_sim_time_lat(ax, sim, **kwargs):

    model = kwargs.get('model')
    timemean = kwargs.get('timemean')
    lat1 = kwargs.get('lat1')
    lat2 = kwargs.get('lat2')

    ax.set_title('%s, %s, $\phi=%g^\circ$ to $%g^\circ$' % (base_title_str(sim, model=model), time_str(timemean), lat1, lat2))
    
def make_title_sim_time_lat_lon(ax, sim, **kwargs):

    model = kwargs.get('model')
    timemean = kwargs.get('timemean')
    lat1 = kwargs.get('lat1')
    lat2 = kwargs.get('lat2')
    lon1 = kwargs.get('lon1')
    lon2 = kwargs.get('lon2')

    ax.set_title('%s, %s\n $\phi=%g^\circ$%s to $%g^\circ$%s, $\lambda=%g^\circ$%s to $%g^\circ$%s' % (base_title_str(sim, model=model), time_str(timemean), abs(lat1), lat_ns(lat1), abs(lat2), lat_ns(lat2), lon1, lon_ew(lon1), lon2, lon_ew(lon2)))
    
def base_title_str(sim, **kwargs):
    
    model = kwargs.get('model')

    if sim == 'rcp85':
        base_str = '%s, %s' % ('RCP8.5', model)
    elif sim == 'historical':
        base_str = '%s, %s' % ('historical', model)
    elif sim == 'hist+rcp85':
        base_str = '%s, %s' % ('RCP8.5', model)
    elif sim == 'longrun':
        if model == 'MPIESM12_control':
            base_str = 'MPI-ESM-LR piControl'
        elif model == 'MPIESM12_abrupt2x':
            base_str = 'LongRunMIP, MPI-ESM-LR abrupt2x'
        elif model == 'MPIESM12_abrupt4x':
            base_str = 'LongRunMIP, MPI-ESM-LR abrupt4x'
        elif model == 'MPIESM12_abrupt8x':
            base_str = 'LongRunMIP, MPI-ESM-LR abrupt8x'
        elif model == 'MPIESM12_abrupt16x':
            base_str = 'LongRunMIP, MPI-ESM-LR abrupt16x'
        elif model == 'MPIESM12_abrupt32x':
            base_str = 'LongRunMIP, MPI-ESM-LR abrupt32x'
    elif sim == 'echam':
        if model == 'rp000140':
            base_str = 'AQUA 15 m'
    elif sim == 'era5':
        base_str = 'ERA5'

    return base_str

def time_str(timemean):
    if timemean == 'yearmean' or timemean == '':
        time_str = 'ANN'
    elif timemean == 'allseas':
        time_str = 'Seasonal Comparison'
    elif timemean == 'djfmean':
        time_str = 'DJF'
    elif timemean == 'mammean':
        time_str = 'MAM'
    elif timemean == 'jjamean':
        time_str = 'JJA'
    elif timemean == 'sonmean':
        time_str = 'SON'
    elif timemean == 'ymonmean-30':
        time_str = 'Last 30 years'
    elif timemean == 'ymonmean-10':
        time_str = 'Last 10 years'
    elif timemean == 'ymonmean':
        time_str = '1979-2019 Climatology'
        
    return time_str

def seas_str(seasmean):
    if seasmean == 'djf':
        seas_str = 'DJF'
    elif seasmean == 'mam':
        seas_str = 'MAM'
    elif seasmean == 'jja':
        seas_str = 'JJA'
    elif seasmean == 'son':
        seas_str = 'SON'
    elif seasmean == '':
        seas_str = 'ANN'

    return seas_str

def lat_ns(lat):

    if lat >= 0:
        latstr = 'N'
    else:
        latstr = 'S'

    return latstr

def lon_ew(lon):

    if lon >= 0:
        lonstr = 'E'
    else:
        lonstr = 'W'

    return lonstr
