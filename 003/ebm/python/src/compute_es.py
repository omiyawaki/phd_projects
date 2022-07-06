import numpy as np

def compute_es(T):

    # This is a Python translation of Todd Mooring's MATLAB script

    # This function computes the saturation vapor pressure (es, units of Pa) as
    # a function of temperature (T, units of K). The formula used is equation
    # (1) of Frierson et al. 2007, which should be standard.

    # 1) load constants
    import par as p

    # 2) compute saturation vapor pressure
    es = p.esref * np.exp( -(p.Lv/p.Rv) * (1/T - 1/p.Tref) )

    return es
