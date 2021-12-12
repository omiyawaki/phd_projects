def get_dimtype(plotvar, **kwargs):

    dimtype_dict = {

            # scalar
            'pr' : 'scalar',
            'prc' : 'scalar',
            'prfrac' : 'scalar',

            # lev
            'ga_dev' : 'lev',


            }

    dimtype = dimtype_dict.get(plotvar, "No dimtype entry is available for this variable. Please add one to dimtype.py")

    return dimtype
