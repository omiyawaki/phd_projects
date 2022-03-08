def translate_varname(varnames_in, **kwargs):
    # rename variable names to the CMIP5 standard

    # key are nonstandard names, value is the standard name
    trans_dict = {
        # MSE storage
        'tend': 'tend',

        # LARGE SCALE CIRCULATION
        'ua': 'ua',
        'va': 'va',
        'wap': 'wap',
        'psi': 'psi',

        # DYNAMICS
        'aht': 'aht',
        'vmmmc': 'vmmmc',
        'vmse': 'vmse',
        'vmte': 'vmte',
        'qaht': 'qaht',
        'vqmmc': 'vqmmc',
        'vqse': 'vqse',
        'vqte': 'vqte',
        'saht': 'saht',
        'vsmmc': 'vsmmc',
        'vsse': 'vsse',
        'vste': 'vste',
        'daht': 'daht',
        'dvmmmc': 'dvmmmc',
        'dvmse': 'dvmse',
        'dvmte': 'dvmte',
        'dqaht': 'dqaht',
        'dvqmmc': 'dvqmmc',
        'dvqse': 'dvqse',
        'dvqte': 'dvqte',
        'dsaht': 'dsaht',
        'dvsmmc': 'dvsmmc',
        'dvsse': 'dvsse',
        'dvste': 'dvste',
        'aht_sm': 'aht_sm',
        'vmmmc_sm': 'vmmmc_sm',
        'vmse_sm': 'vmse_sm',
        'vmte_sm': 'vmte_sm',
        'qaht_sm': 'qaht_sm',
        'vqmmc_sm': 'vqmmc_sm',
        'vqse_sm': 'vqse_sm',
        'vqte_sm': 'vqte_sm',
        'saht_sm': 'saht_sm',
        'vsmmc_sm': 'vsmmc_sm',
        'vsse_sm': 'vsse_sm',
        'vste_sm': 'vste_sm',
        'daht_sm': 'daht_sm',
        'dvmmmc_sm': 'dvmmmc_sm',
        'dvmse_sm': 'dvmse_sm',
        'dvmte_sm': 'dvmte_sm',
        'dqaht_sm': 'dqaht_sm',
        'dvqmmc_sm': 'dvqmmc_sm',
        'dvqse_sm': 'dvqse_sm',
        'dvqte_sm': 'dvqte_sm',
        'dsaht_sm': 'dsaht_sm',
        'dvsmmc_sm': 'dvsmmc_sm',
        'dvsse_sm': 'dvsse_sm',
        'dvste_sm': 'dvste_sm',

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
        'rsutcs' : 'rsutcs',
        'rlutcs' : 'rlutcs',
        'rsdscs' : 'rsdscs',
        'rsuscs' : 'rsuscs',
        'rldscs' : 'rldscs',

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

        # SPECIFIC HUMIDITY
        'hus' : 'hus', # CMIP

        # RELATIVE HUMIDITY
        'hur' : 'hur', # CMIP

        # VERTICAL MEAN RELATIVE HUMIDITY
        'vhur' : 'vhur', # CMIP

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
        
        # PRECIPITATION
        'pr' : 'pr', # CMIP
        
        # EVAPORATION
        'evspsbl' : 'evspsbl', # CMIP
        
        # CONVECTIVE PRECIP
        'prc' : 'prc', # CMIP
        'cp' : 'prc', # ERA

        # PRECIPITABLE WATER
        'prw' : 'prw', # CMIP

        # total cloud fraction
        'clt' : 'clt', # CMIP

        # condensed water path
        'clwvi' : 'clwvi', # CMIP

        # condensed ice path
        'clivi' : 'clivi', # CMIP
        
        # LARGE-SCALE PRECIP
        'lsp' : 'prl', # ERA

        # SEA ICE CONCENTRATION
        'sic' : 'sic', # CMIP
        
        # CO2 MASS
        'co2mass' : 'co2mass', # CMIP

        # KINETIC ENERGY
        'tke': 'tke', # total KE
        'mke': 'mke', # mean meridional KE
        'teke': 'teke', # transient eddy KE
        'seke': 'seke', # stationary eddy KE
        'eke': 'eke', # eddy KE
        'reke': 'reke', # eddy KE (inferred)
        'res': 'res', # residual KE
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
        'zs' : 'orog',
        
        # HYDROLOGY
        # ERA
        'cp' : 'hydro',
        'lsp' : 'hydro',
        'e' : 'hydro',

        # SEA ICE
        # ERA
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

def translate_vardim(varname_std, **kwargs):
    # outputs whether a variable is 2D (horizontal) or 3D (+vertical) in space
    trans_dict = {
        'ps' : '2d',
        'ts' : '2d',
        'tas' : '2d',
        'orog' : '2d',

        'ta' : '3d',
        'zg' : '3d',
    }

    vardim = trans_dict.get(varname_std, "No translation available")

    return vardim

def translate_varsfc(varname_std, **kwargs):
    # outputs the surface variable name of a 3D (vertical) variable
    trans_dict = {
        'ta' : 'tas',
        'zg' : 'orog',
    }

    varsfc = trans_dict.get(varname_std, "No translation available")

    return varsfc

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
