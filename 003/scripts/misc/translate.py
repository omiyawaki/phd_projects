def translate_varname(varnames_in, **kwargs):
    # rename variable names to the CMIP5 standard

    # key are nonstandard names, value is the standard name
    trans_dict = {
        'srad0': 'srad0',
        'trad0': 'trad0',
        'srads': 'srads',
        'trads': 'trads',

        'rsdt' : 'rsdt',
        'rsut' : 'rsut',
        'rlut' : 'rlut',
        'rsds' : 'rsds',
        'rsus' : 'rsus',
        'rlds' : 'rlds',
        'rlus' : 'rlus',

        'hfls' : 'hfls',
        'ahfl' : 'hfls',

        'hfss' : 'hfss',
        'ahfs' : 'hfss',

        'sic' : 'sic'
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
        'ssr': 'rad',
        'str': 'rad',
        'tsr': 'rad',
        'ttr': 'rad',

        'slhf' : 'stf',
        'sshf' : 'stf',
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
