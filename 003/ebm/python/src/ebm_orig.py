import sys
sys.path.append('/project2/tas1/miyawaki/projects/003/ebm/python')
from tqdm import tqdm
import src.par as p
import numpy as np

def ebm_orig(SW, **kwargs):

    # This is a Python translation of Todd Mooring's MATLAB script

    # This function is intended as an easily runnable form of (basically) the
    # Hwang and Frierson EBM supplied by Tiffany in a 7/8/17 email. Some physical
    # constants differ slightly from their Hwang and Frierson values and others
    # are now set as inputs to the function, but the model structure and solution
    # algorithm are exactly the same.

    # REQUIRED INPUT: 
    # SW    Specified components of the MSE generation rate -- i.e., everything
    #       other than clear sky OLR. One way to break this down is as (net
    #       downward shortwave at the TOA + longwave CRF - net downward energy
    #       flux at the surface). Units are W/m^2. The specified values are assumed
    #       to be defined at the centers of gridboxes of uniform meridional size
    #       and the first element is assumed to be closest to the south pole.

    # OPTIONAL INPUTS:
    # a     Coefficient in aT-b, the linearized formula for clear sky OLR.
    #       A scalar. Default value is 2.07 W/m^2/K.
    # b     Coefficient in aT-b. A scalar. Default value is 332.4 W/m^2.
    # D     Diffusivity. A scalar. Default value is 1.06e6 m^2/s.
    # hs    Relative humidity. A scalar. Default is 0.8 (unitless).
    # sfc   Surface fluxes. Positive values cool the atmosphere.
    #       Defined in full latitudes. Default value is 0 (W/m^2).
    # stg   Atmospheric MSE storage. Positive values cool the atmosphere.
    #       Defined in full latitudes. Default is 0 (W/m^2).
    # vb    Logical, prints To every loop if true.

    # OUTPUTS:
    # To    Final temperature profile, K.
    # lat   Assumed latitudes, deg.
    # qo    Final specific humidity profile, dimensionless.
    # mo    Final MSE profile, J/kg.
    # trn   Energy transport on the half latitudes, W.
    # lath  Assumed half latitudes, deg.
    # nl    Number of iterations required for model to converge.

    # 1) Set default values
    a = kwargs.get('a', 2.07) # W/m^2/K
    b = kwargs.get('b', 332.4) # W/m^2
    D = kwargs.get('D', 1.06e6) # m^2/s
    hs = kwargs.get('hs', 0.8) # unitless
    sfc = kwargs.get('sfc', 0) # W/m^2 (surface flux, positive value cools atmos)
    stg = kwargs.get('stg', 0) # W/m^2 (atmospheric mse storage, positive value cools atmos)
    vb = kwargs.get('vb', False) # logical

    # 2) Set non-adjustable parameters
    nlat = len(SW)
    glat = np.pi/nlat # latitude grid spacing in radians
    lat = -1/2 * (glat+np.pi) + glat*(np.arange(1,nlat+1,1)) # np.arange(1,nlat+1,1) makes a vector of the correct length and spacing 1. Multiplying by glat changes spacing to glat. Subtracting pi/2 starts the array one gridbox north of the south pole. Subtracting glat/2 starts the array half a gridbox north of the south pole, as needed.
    To = 200*np.ones(nlat) # initial condition for temperature, K.
    dt = 10**-4 # time step, will assume it's years (i.e. 0.876 days)
    ps = 1e4*p.g # surface pressure, Pa
    Do = ps*D/p.g # Normalized diffusion coefficient, kg/s
    Ta = p.Tref # CC reference temp, K.
    eo = p.esref # CC reference vapor pressure, Pa.
    c = 1 # ocean heat capacity [Actually this parameter doesn't really have physical meaning except in the case where we are simulating a slab ocean model with no heat transport. However this doesn't matter if we are only interested in the model's steady-state behavior. Given our assumption about the units of dt, this has units of W yr/^2/K. 1 W yr/m^2/K = 3.1536e7 J/m^2/K.

    A = np.zeros([nlat,nlat]) # square array of nlat
    for k in np.arange(0,nlat-1,1):
        # Iterate over all but last latitude
        A[k,k] = -1
        A[k,k+1] = 1

    B = np.zeros([nlat,nlat]) # same dimensions as A
    for k in np.arange(1,nlat,1):
        # Iterate over all but the first latitude
        B[k,k] = 1
        B[k,k-1] = -1

    lat_ph = lat + glat/2 
    lat_nh = lat - glat/2 # For now let's assume that glat is the spacing of the latitude grid points in radians and that this is uniform. Thus lat_ph (lat_nh) can be considered a bector of northern (southern) edges of grid cells (assuming latitude is in radians N).

    # 3) Run EBM
    for n in tqdm(np.arange(1e6)):
        # Iterate over a very large number of steps
        es = eo*np.exp((-p.Lv/p.Rv)*(1/To - 1/Ta)) # Compute saturation vapor pressure given current temperature To
        mo = p.cpair*To + p.Lv*hs*p.Rd*es/(p.Rv*ps) # Moist static energy cp*T + L*q assuming that the specific humidity and mixing ratios are approximately equal and that the surface pressure is created entirely by dry air. z = 0 because we are evaluating this at ground level.
        OLR = a*To - b # OLR, W/m^2

        # general form of temperature tendency due to diffusion (horizontal transport)
        Diff = Do*p.a**(-2)*(np.cos(lat)**(-1)*(glat)**(-2)*(np.cos(lat_ph)*np.dot(A,mo)-np.cos(lat_nh)*np.dot(B,mo)))

        # For south boundary (first grid point)
        Diff[0] = Do*p.a**(-2)*(np.cos(lat[0]))**(-1)*glat**(-2)*(np.cos(lat[1])*(mo[2]-mo[1]) - np.cos(lat[0])*(mo[1]-mo[0]))

        # For north boundary (last grid point)
        Diff[-1] = Do*p.a**(-2)*(np.cos(lat[-1]))**(-1)*(glat)**(-2)*(np.cos(lat[-1])*(mo[-1]-mo[-2]) - np.cos(lat[-2])*(mo[-2]-mo[-3]))

        Tn = To + (dt*(SW-OLR+Diff-sfc-stg))/c # Technically this functional form is appropriate only to the case where column heat capacity is dominated by the ocean in the slab ocean/no ocean heat transport case. If this is not so (e.g., when seeking to emulate an AGCM), the moist static energy capacity dm/dT that relates the moist static energy change rate dm/dt to the temperature change rate dT/dt has a well-defined (albeit in principle temperature dependent) value which will not generally equal the specified c. However, heat capacity is not relevant to the equilibrium behavior of the model and thus this equation can be considered a heuristic adjustment that moves the system closer to equilibrium -- warmer if MSE is being generated/converged into the latitude band, cooler if MSE is being destroyed/diverged from the latitude band.
        To = Tn # updating surface temperature

        if vb:
            # print current To if in verbose mode
            print(To)

        # break the loop if equilibrium is reached
        W = (np.max(np.abs(SW-OLR+Diff-sfc-stg))<0.001 and n>0) # Set to 1 if the magnitude of the MSE change rate at any given latitude is < 0.001 W/m^2 and this is not the first passage through the loop, 0 otherwise.

        if W==1:
            break # stop iterating if convergence has been achieved
        
    # 4) Create output variables
    lat = np.rad2deg(lat) # change unit to deg
    es = eo*np.exp((-p.Lv/p.Rv)*(1/To - 1/Ta)) # compute saturation vapor pressure given current temperature To
    qo = hs*p.Rd*es/(p.Rv*ps) # compute specific humidity
    mo = p.cpair*To + p.Lv*qo # compute MSE
    lath = 1/2 * (lat[1:]+lat[:-1]) # half lats at which transport is defined
    trn = -2*np.pi*Do*np.cos(np.deg2rad(lath)) * (mo[1:]-mo[:-1])/glat # compute transports on the half latitudes, units of W
    nl = n

    return To, lat, qo, mo, trn, lath, nl

