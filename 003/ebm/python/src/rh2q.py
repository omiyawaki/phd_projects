import numpy as np

def rh2q(rh,T,p):

    # This is a Python translation of Todd Mooring's MATLAB script

    # This function converts relative humidities to specific humidities via a
    # formula derived in my 3/23/17 notes. 

    # INPUTS:
    # rh    relative humidity, unitless (values should range from 0 to 1)
    # T     temperature, units K
    # p     total air pressure, units Pa

    # OUTPUT:
    # q     specific humidity, unitless

    # 1) load constants
    import par as p
    ep = p.Rv/p.R

    # 2) get saturation vapor pressure
    from compute_es import compute_es
    es = compute_es(T);

    # 3) compute specific humidity
    q = 1/(1+(ep*(p/(rh*es)-1)))

    return q
