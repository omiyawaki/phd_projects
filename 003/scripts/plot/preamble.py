def get_predata(sim, timemean, kwargs):

    if sim == 'longrun':
        model = kwargs.get('model', 'MPIESM12_abrupt4x')
        yr_span = kwargs.get('yr_span', '1000')
        yr_base = 0
        yr_base_show = yr_base
        yr_end_show = 2200

    elif sim == 'ssp585':
        model = kwargs.get('model', 'MPI-ESM-LR')
        yr_span = kwargs.get('yr_span', '201501-230012')
        yr_span_ref = kwargs.get('yr_span_ref', '186001-201412')
        if 'ymonmean' not in timemean:
            yr_base = 2015
            yr_base_show = yr_base
            yr_end_show = 2200
        else:
            yr_base = 0
            yr_base_show = yr_base
            yr_end_show = 2200

    elif sim == 'rcp85':
        model = kwargs.get('model', 'MPI-ESM-LR')
        yr_span = kwargs.get('yr_span', '200601-230012')
        yr_span_ref = kwargs.get('yr_span_ref', '186001-200512')
        if 'ymonmean' not in timemean:
            yr_base = 2006
            yr_base_show = yr_base
            yr_end_show = 2200
        else:
            yr_base = 0
            yr_base_show = yr_base
            yr_end_show = 2200

    elif sim == 'hist+ssp585':
        model = kwargs.get('model', 'MPI-ESM-LR')
        yr_span = kwargs.get('yr_span', '186001-230012')
        yr_span_ref = kwargs.get('yr_span_ref', '186001-201412')
        if 'ymonmean' not in timemean:
            yr_base = 1860
            yr_base_show = 1990
            yr_end_show = 2200
        else:
            yr_base = 0
            yr_base_show = yr_base
            yr_end_show = 2200

    elif sim == 'hist+rcp85':
        model = kwargs.get('model', 'MPI-ESM-LR')
        yr_span = kwargs.get('yr_span', '186001-229912')
        yr_span_ref = kwargs.get('yr_span_ref', '186001-200512')
        if 'ymonmean' not in timemean:
            yr_base = 1860
            yr_base_show = 1990
            yr_end_show = 2200
        else:
            yr_base = 0
            yr_base_show = yr_base
            yr_end_show = 2200

    elif sim == 'historical':
        model = kwargs.get('model', 'MPI-ESM-LR')
        yr_span = kwargs.get('yr_span', '186001-200512')
        if 'ymonmean' not in timemean:
            yr_base = 1860
            yr_base_show = yr_base
            yr_end_show = 2200
        else:
            yr_base = 0
            yr_base_show = yr_base
            yr_end_show = 2200

    elif sim == 'echam':
        model = kwargs.get('model', 'rp000140')
        yr_span = kwargs.get('yr_span', '0001_0039')
        yr_span_ref = kwargs.get('yr_span_ref', '186001-201412')

        if model in ['rp000130', 'rp000131', 'rp000134', 'rp000135', 'rp000190','rp000190f']:
            yr_base = 20
            yr_base_show = yr_base
            yr_end_show = 2200
        elif model in ['rp000184', 'rp000185', 'rp000188', 'rp000189', 'rp000191', 'rp000191b','rp000191f']:
            yr_base = 1987
            yr_base_show = 1990
            yr_end_show = 2200

    elif sim == 'era5':
        model = None
        yr_span = kwargs.get('yr_span', '1979_2019')
        yr_base = 1979
        yr_base_show = yr_base
        yr_end_show = 2200

    return model, yr_span, yr_base, yr_span_ref, yr_base_show, yr_end_show
