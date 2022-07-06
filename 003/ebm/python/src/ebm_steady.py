import par as p
import numpy as np

def ebm_steady(N, Do, **kwargs):

    # This is a Python translation of Todd Mooring's MATLAB script

    # This function is an energy balance model very similar to that used by
    # Hwang and Frierson 2010. Key differences are the following:

    # 1) The model explicitly computes a steady-state solution, rather than
    #    iterating until convergence is achieved.
    # 2) The net MSE generation rate with which the MSE flux divergence is in
    #    balance is directly specified in its entirety, rather than being a
    #    function of temperature. Thus there is no need for explicit reference
    #    to temperature or any other specific subcomponent of moist static energy
    #    anywhere in the code. However, to get the assumed steady state, the
    #    global mean MSE generation rate must be zero -- this is imposed by
    #    making a globally uniform adjustment to the input N before using it to
    #    drive the model. This also implies that the model's global mean MSE will
    #    remain fixed at whatever value it was initialized at. Since global mean
    #    MSE is an input parameter of this model, rather than an output, the MSE
    #    gradient is arguably more physically meaningful and the input global mean
    #    MSE must be chosen with care if the full MSE profile is to be physically
    #    interpreted.
    # 3) Diffusivity is allowed to vary with latitude.
    # 4) Numerics are more explicitly finite-volume, in order to obtain the desired
    #    conservation properties.

    # INPUTS (REQUIRED):
    # N     Zonal mean MSE generation rate in W/m^2. Defined on N evenly-spaced
    #       full latitudes stretching from pole to pole -- given this assumption,
    #       the grid boxes have meridional extend 180 deg N and the most southerly
    #       box has its center at -90 to +90 deg N.
    # Do    Zonal mean diffusivity in kg/s. ** This is the EBM definition of 
    #       diffusivity, not the D = -H/G diagnosed diffusivity. ** Must be
    #       specified on N-1 half latitudes evenly spaced between the N full
    #       latitudes. Can alternatively just pass a scalar, which will be
    #       automatically expanded.

    # INPUTS (OPTIONAL):
    # mi    Initial MSE profile, units J/kg. Can be defined on the N full
    #       latitudes, or by passing a scalar to specify the global mean MSE.
    #       Default value is 0.

    # OUPTUTS:
    # dm_dy MSE gradient, units J/kg/m. Defined on half latitudes.
    # lath  Half latitudes, deg N.
    # m     MSE, units J/kg. Defined on full latitudes.
    # latf  Full latitudes, deg N.
    # H     MSE transport, W. Defined on half latitudes.
    # Nim   Global mean MSE generation rate associated with the input N
    # Nout  MSE generation rate profile with the global mean removed to enable
    #       appropriate forcing of the model.
    # nc    Degree of non-closure in transport calculations. Two-element row
    #       vector, left (right) element is implied flux at north (south) pole
    #       from south-to-north (north-to-south) transport integral. Should be
    #       very small due to explicit adjustment of N to have zero global mean.
    #       Units are W/m^2.

    # 1) Validate inputs
    mi = kwargs.get('mi', 0) # J/kg

    try:
        if not ( (len(N) == len(Do)+1) or (len(Do)==1) ):
            raise DoLenError
        if not ( (len(N) == len(mi)) or (len(mi)==1) ):
            raise miLenError
    except DoLenError:
        print('Do (diffusivity) has incorrect length.') 
    except miLenError:
        print('mi (initial MSE profile) has incorrect length.') 

    # 2) Remove global mean MSE generation rate
    nlats = len(N) # Number of full latitudes
    lath = np.linspace(-90,90,nlats+1) # Create half latitudes
    latf = 1/2 * (lath[1:] + lath[:-1]) # Create full latitudes
    dsinl = np.sin(np.deg2rad(lath[1:])) - np.sin(np.deg2rad(lath[:-1]))
    Nim = 1/2 * np.sum(N * dsinl) # Global mean MSE generation rate, W/m^2
    Nout = N - Nim # Now has zero global mean for further use

    # 3) Compute MSE transports
    diff_H = 2*np.pi*p.a**2*(Nout*dsinl) # Difference between meridional MSE transports at successive half latitudes, W
    s2n = np.concatenate(([0], np.cumsum(diff_H))) # South to north
    n2s = np.concatenate((np.cumsum(diff_H[::-1]), [0])) # North to south
    nc = np.concatenate((s2n[-1], n2s[0]))/(4*np.pi*p.a**2) # Store degree of non-closure
    H = 1/2 * (s2n + n2s) # Average across both integration directions, W
    H = H[1:-1] # Do not return transports at the poles

    # 4) Compute MSE gradients
    lath = lath[1:-1] # Drop poles
    dm_dy = -H / (2*np.pi*p.a*(Do*np.cos(np.deg2rad(lath)))) # Employs definition of diffusivity to compute the gradient, J/kg/m

    # 5) Compute MSE profile
    dm = ((np.pi*p.a)/nlats)*dm_dy # MSE increment between successive boxes
    m = np.concatenate(([0], np.cumsum(dm))) # MSE increment relative to the value in the most southerly grid box, J/kg
    mm = 1/2 * np.sum(m*dsinl) # Global mean of the MSE increment profile
    m_im = 1/2 * np.sum(mi*dsinl) # Global mean of initial MSE profile
    m = m - mm + m_im # Adjust MSE increment profile so that it has the same mean as the initial MSE profile, i.e. could have been produced by rearranging the initial condition, J/kg.

    return dm_dy, lath, m, hatf, H, Nim, Nout, nc
