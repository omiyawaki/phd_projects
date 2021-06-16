def translate_varname(varnames_in, **kwargs):
    # rename variable names to the CMIP5 standard

    # key are nonstandard names, value is the standard name
    trans_dict = {
        # RADIATION
        # ECHAM
        'srad0': 'srad0',
        'trad0': 'trad0',
        'srads': 'srads',
        'trads': 'trads',
        # ERA
        'tsr': 'srad0',
        'ttr': 'trad0',
        'ssr': 'srads',
        'str': 'trads',
        # CMIP
        'rsdt' : 'rsdt',
        'rsut' : 'rsut',
        'rlut' : 'rlut',
        'rsds' : 'rsds',
        'rsus' : 'rsus',
        'rlds' : 'rlds',
        'rlus' : 'rlus',

        # LATENT HEAT
        'hfls' : 'hfls', # CMIP
        'ahfl' : 'hfls', # ECHAM
        'slhf' : 'hfls', # ERA

        # SENSIBLE HEAT
        'hfss' : 'hfss', # CMIP
        'ahfs' : 'hfss', # ECHAM
        'sshf' : 'hfss', # ERA

        # TEMPERATURE
        'ta' : 'ta', # CMIP
        't' : 'ta', # ERA

        # GEOPOTENTIAL HEIGHT
        'zg' : 'zg', # CMIP
        'z' : 'zg', # ERA

        # SURFACE PRESSURE
        'ps' : 'ps', # CMIP
        'sp' : 'ps', # ERA

        # SURFACE HEIGHT
        'orog' : 'orog', # CMIP
        'zs' : 'orog', # ERA

        # SURFACE TEMPERATURE
        'ts' : 'ts', # CMIP
        'skt' : 'ts', # ERA

        # 2 M TEMPERATURE
        'tas' : 'tas', # CMIP
        't2m' : 'tas', # ERA

        # 2 M DEW POINT TEMP
        'd2m' : 'd2m', # ERA

        # SEA ICE CONCENTRATION
        'sic' : 'sic' # CMIP
    }

    # if input is a list of variables, output as list of variables
    if isinstance(varnames_in, list):
        varnames_out = []
        for i in range(len(varnames_in)):
            varnames_out[i] = trans_dict.get(varnames_in[i], "No translation available")

    # if input is a single string, output a string
    elif isinstance(varnames_in, str):
        varnames_out = trans_dict.get(varnames_in, "No translation available")

    return varnames_out

def translate_varcat(varnames_in, **kwargs):
    # category of variable for loading data

    # key are variable names, value is the category name
    trans_dict = {
        # RADIATION
        # ERA
        'ssr': 'rad',
        'str': 'rad',
        'tsr': 'rad',
        'ttr': 'rad',

        # SURFACE TURBULENT FLUXES
        # ERA
        'slhf' : 'stf',
        'sshf' : 'stf',

        # TEMPERATURE
        't' : 'temp', # ERA

        # GEOPOTENTIAL HEIGHT
        'z' : 'zg', # ERA

        # SURFACE
        # ERA
        'sp' : 'srfc',
        't2m' : 'srfc',
        'skt' : 'srfc',

        # OROGRAPHY
        'zs' : 'orog'
    }

    # if input is a list of variables, output as list of variables
    if isinstance(varnames_in, list):
        varcat = []
        for i in range(len(varnames_in)):
            varcat[i] = trans_dict.get(varnames_in[i], "No translation available")

    # if input is a single string, output a string
    elif isinstance(varnames_in, str):
        varcat = trans_dict.get(varnames_in, "No translation available")

    return varcat

def latetrans_grid(sim, gridname_in, **kwargs):
    # outputs non standard grid names for reanalyses/models
    if sim == 'era5':
        if gridname_in == 'lev':
            gridname_out = 'level'
        elif gridname_in == 'lat':
            gridname_out = 'lat'
        elif gridname_in == 'lon':
            gridname_out = 'lon'
    else:
        if gridname_in == 'lev':
            gridname_out = 'plev'
        elif gridname_in == 'lat':
            gridname_out = 'lat'
        elif gridname_in == 'lon':
            gridname_out = 'lon'

    return gridname_out