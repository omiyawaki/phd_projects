import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/ebm/python')
import time
import numpy as np
import src.par as p
from src.rh2q import rh2q
from scipy.interpolate import interp1d

def ebm_flux(Ns, lat, **kwargs):

    # This is a Python translation of Todd Mooring's MATLAB script

    # This function is an energy balance model very similar to that used by
    # Hwang and Frierson 2010. Key differences are the following:
    # 1) Overall numerical algorithm has been carefully formulated in flux form
    #    to improve the conservation properties of the model, in contrast
    #    to HF2010's finite difference formulation.
    # 2) The boundary condition (no flux through the poles) makes more sense
    #    than the HF2010 boundary condition which I don't fully understand.
    # 3) Diffusivity can vary with latitude and is input in units of kg/s.
    # 4) Relative Humidity can vary with latitude.
    # 5) A more accurate formula for computing specific humidity from relative
    #    humidity and temperature is used.
    # 6) A geopotential height contribution to the mean state moist static
    #    energy can be specified.
    # 7) Default first guess temperature profile is 286.65+13.5*cosd(2*lat),
    #    inspired by the Neale and Hoskins 2001 Qobs profile. Can
    #    alternatively be specified by the user as an input.
    # 8) Automated selection of the timestep used in the iterative process,
    #    based on making sure that temperature values never exceed the
    #    0-400 K range and restarting the iterative process with a smaller
    #    timestep if they do.

    # REQUIRED INPUTS:
    # Ns    Specific components of the MSE generation rate -- i.e., everything
    #       other than clear sky OLR. One way to break this down is as
    #       (net downward shortwave at the TOA + longwave CRF - net
    #       downward energy flux at the surface). Units are W/m^2. The specified
    #       values are assumed to be defined on the full latitude grid.
    # lat   Vector of full latitudes, units are deg N.

    # OPTIONAL INPUTS:
    # a     Coefficient in aT-b, the linearized formula for clear sky OLR. 
    #       A scalar. Default value is 2.07 W/m^2/K.
    # b     Coefficient in aT-b. A scalar. Default value is 332.4 W/m^2.
    # Do    Diffusivity. A scalar or vector. If the vector has the same number
    #       of elements as lat, Do is defined on the full latitudes and will be
    #       interpolated to the half latitudes via a spline. If the column
    #       vector has one fewer element than lat, it is assumed to be defined
    #       directly on the half latitudes. Default value is 1.06e10 kg/s. 
    #       ** This is the EBM definition of diffusivity, not the D = -H/G
    #       dianosed diffusivity. The two diffusivities are related by
    #       Do = D/cos(lat).**
    # ptot  Total air pressure to use in the conversion from relative to
    #       specific humidity. Default value is 92500 Pa.
    # rh    Relative humidity. A scalar or vector defined on full latitudes.
    #       Default value is 0.8.
    # zg    Specified geopotential height contribution to the mean state moist
    #       static energy. A scalar or vector defined on full latitudes.
    #       Default value is 0 m.
    # Too   First guess temperature profile used to initialize the iterative
    #       process. Default value is 286.65+13.5*cosd(2*lat) K.
    # dt_c  First guess of the ratio of the timestep to assumed system heat
    #       capacity, which controls the temperature change made in response
    #       to a given total heating rate. If iteration does not appear to be
    #       converging (i.e., departs from the 0-400 K temperature range),
    #       dt_c is halved and the iterative process is restarted. This
    #       process of halving is repeated until convergence is achieved or
    #       the maximum number of steps for a given value of dt_c is exceeded.
    #       Default value is 2^(-11) Km^2/W.
    # max_s First guess of the maximum allowed number of iteration steps, which
    #       is doubled for each necessary halving of dt_c (i.e., the time for
    #       which the integration is run doesn't change. The steps are just made
    #       smaller). If max_s is reached before converence to within the target
    #       accuracy is achieved, the iteration process will quit anyway.
    #       Default value is 2^22.
    # cnv_t Tolerance to within which net heating rate must converge everywhere
    #       in the domain to consider the model to have reached a steady state.
    #       Default value is 0.001 W/m^2.

    # OUTPUTS:
    # To    Final temperature profile, K.
    # qo    Final specific humidity profile, dimensionless.
    # mo    Final MSE profile, J/kg.
    # H     Energy transport on half latitudes, W.
    # lath  Assumed half latitudes, deg.
    # nl    Number of iterations required for model to converge. This will be
    #       negative if iteration process stopped before convergence was achieved.
    # dt_c  Ratio of timestep to assumed system heat capacity actually used for
    #       successful iteration, Km^2/W.
    # ll    Length of time in seconds required to run the final (successful)
    #       iteration attempt.
    # fi    Length of time required to run all of the iteration attempts.
    # na    Number of iteration attempts made (i.e., number of dt_c/max_s
    #       combinations used before convergence was achieved or max_s steps
    #       were taken without achieving convergence).

    # 1) Set default values for optional values
    a = kwargs.get('a', 2.07) # W/m^2/K
    b = kwargs.get('b', 332.4) # W/m^2
    Do = kwargs.get('Do', 1.06e10*np.ones(len(lat)-1)) # kg/s, half lat
    ptot = kwars.get('ptot', 92500) # Pa
    rh = kwargs.get('rh', 0.8*np.ones(len(lat))) # dimensionless, full lat
    zg = kwargs.get('zg', np.zeros(len(lat))) # full lat
    Too = kwargs.get('Too', 286.65+13.5*np.cos(np.deg2rad(2*lat))) # 27C at equator and 0 C at the poles. Same range as Neale and Hoskins 2001 Qobs, but much smoother.
    dt_c = kwargs.get('dt_c', 2^(-11)) # Proportionality constant controlling the relationship between the TOA energy imbalance and the adjustment to the low-level temperature, units Km^2/W. Larger value should speed integration, unless it causes the model to become unstable or to oscillate around the true solution. If oscillatory behavior is detected, reduce this value to enhance stability. This is basically a ratio of a timestep to a heat capacity.
    max_s = kwargs.get('max_s', 2^22) # Perform up to 2^22 iteration steps while seeking steady solution to the model.
    cnv_t = kwargs.get('cnv_t', 0.001) # Tolerance to within which net heating rate must converge everywhere in the domain to consider the model to have reached a steady state. Units W/m^2.
    zg = p.g * zg # Convert geopotential height units to J/kg

    # 2) Expand diffusivity to vector as needed
    lath = 0.5*(lat[:-1]+lat[1:])
    if lath[0] < 0:
        # south to north
        lath = np.concatenate(([-90],lath,[90]))
    else:
        lath = np.concatenate(([90],lath,[-90]))

    try:
        if not isinstance(Do, list)
            # diffusivity is globally uniform
            Do = Do*np.ones(len(lat)-1) # Diffusivity should be defined on the half latitudes that are not the poles
        elif len(Do) == len(lat) :
            # diffusivity is defined on N full latitudes
            fDo = interp1d(lat, Do, kind='cubic')
            Do = fDo(lath[1:-1]) # Do is now on N-1 half lats
        elif not (len(Do) == len(lat)-1) :
            # diffusivity is defined on N-1 half latitudes
            raise DoLenError
    except DoLenError:
        print('Unexpected number of Do latitudes!')

    # 3) Initialize EBM variables
    H = np.zeros(len(lath)) # MSE transport divided by 2*pi on N+1 half lats
    dlat = np.deg2rad(lat[1:]-lat[:-1]) # Full latitude spacings on N-1 half lats
    awt = np.sin(np.deg2rad(lath[1:]))-np.sin(np.deg2rad(lath[:-1])) # Differences of successive sines of half latitudes -- this quantity is defined on N full latitudes and is proportional to the surface area of each full latitude strip
    coslath = np.cos(np.deg2rad(lath[1:-1])) # Cosines of N-1 half latitudes
    H_pf = -((coslath*Do)/dlat) # Transform MSE differences to transports divided by 2*pi, N-1 half lats
    inv_dA = 1/(p.a**2*awt) # This quantity is proportional to the inverse of the surface area of each full latitude strip. N full latitudes, units m**-2.

    # 4) Run EBM
    ic = False # Iteration is not finished
    na = 0 # Number of attempts at finding iterative solution
    fi0 = time.time() # Time the entire iteration process

    while not ic:
        # Experiment with different values of dt_c until stable solution is found
        To = Too # Initialize model run with temperature profile Too
        na = na+1
        ll0 = time.time() # Time the final loop of the iteration process

        for ii in range(len(max_s)):
            # Take up to max_s iteration steps
            qo = rh2q(rh,To,ptot) # Compute specific humidity
            mo = cp*To + p.Lv*qo + zg # Compute low level MSE on N full lats
            H[1:-1] = H_pf*(mo[1:]-mo[:-1]) # Compute MSE transports divided by 2*pi from differences in adjacent MSE gridboxes, N-1 half lats
            Nh = inv_dA*(H[1:,...]-H[:-1,...]) # Heating rate of each full gridbox due to transport, W/m^2
            Ncs = b-a*To # Heating rate of each full gridbox due to clear sky OLR, W/m^2
            Nt = Ns+Ncs+Nh # Total heating rate of each full gridbox, W/m^2

            if np.max(np.abs(Nt)) < cnv_t
                # Max heating rate associated with this MSE profile is < cnv_t W/m^2
                ll = time.time() - ll0
                fi = time.time() - fi0
                ic = True
                break # When this happens, we consider the iteration converged. Break exits the for loop, while resetting ic causes exit from the while loop.

            To = To + dt_c*Nt # Update temperatures

            if np.any(To>400) or np.any(To<0)
                # Check to see if temperatures have become unreasonable
                dt_c = dt_c/2 # Shorten timestep
                max_s = 2*max_s # Boost maximum allowed number of iterations
                break # Leaves the for loop, reinitializes temperature profile

        if ii == max_s
            # Finished going through for loop but did not converge
            ll = time.time() - ll0
            fi = time.time() - fi0
            ic = True # Not really, but we want to terminate the run
            ii = -ii # Flip the sign of ii (and by extension nl) to indicate the lack of convergence

    # 5) Create variables for return
    H = 2*np.pi * H # Convert to actual transport in units of W on N+1 half lats
    H = H[1:-1] # Reduce to N-1 half lats
    lath = lath[1:-1]
    nl = ii # Number of passages through the loop required to converge

    return To, qo, mo, H, lath, nl, dt_c, ll, fi, na
