import numpy as np

def calc_esat(T):
    esat = 611 * np.exp(6808 * (1/273 - 1/T) - 5.09 * np.log(T/273))
    return esat

def e2q(p,e):
    q = 0.622 * e / (p - e*(1-0.622))
    return q

def get_modelidx(model):
    idxtable = {
            'bcc-csm1-1': 0,
            'CCSM4': 1,
            'CNRM-CM5': 2,
            'CSIRO-Mk3-6-0': 3,
            'IPSL-CM5A-LR': 4,
            'HadGEM2-ES': 5,
            'MPI-ESM-LR': 6
            }

    idx = idxtable.get(model, '%s is not a valid model name.' % (model) )

    return idx

def get_modelstr(model):
    if model == 'mmm':
        modelstr = 'CMIP5 mean'
    else:
        modelstr = model

    return modelstr
