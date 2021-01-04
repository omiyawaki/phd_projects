import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import splrep, splev

def calc_ma(ps, tas, hurs, plev, frz=None, ma_type=None, dpa=None, pa_top=None, eps=None, cpd=None, cpv=None, cpl=None, cpi=None,
            Rd=None, Rv=None, g=None):
    """
    This function calculates the moist adiabatic temperature profile ta_s given the following inputs:
        ps: surface pressure (Pa), type: float
        tas: 2 m temperature (K), type: float
        hurs: 2 m humidity (%), type: float

    The fourth argument, plev (type: 1D array), specifies the vertical pressure levels where ta_s should be reported.

    The fifth argument, (type: struct), specifies meters and constants that are used when computing the moist adiabat. These include:
        frz: specifies whether or not to include the latent heat of fusion from freezing.
            0 = no freezing (default behavior)
            1 = freezing happens immediately at 273 K
        ma_type: type of moist adiabat
            'std' = standard adiabat; simplified form of the pseudoadiabat where the water vapor mixing ratio is assumed to be << 1. This is the standard definition in the AMS glossary.
            'pseudo' = pseudoadiabat; all condensates precipitate out immediately.
            'reversible' = reversible adiabat; all condensates remain in the el. Conserves entropy
        dpa: the step size when integrating the dry adiabat up to LCL. The smaller the step size, the more accurate the LCL is calculated. By default, this is set to -10 Pa.
        pa_top: the topmost pressure level to calculate the moist adiabat up to. By default, 100 hPa.
        eps: the ratio of Rd to Rv. Set to 0.622 by default.
        cpd: isobaric specific heat capacity of dry air. Set to 1005.7 J/kg/K by default.
        cpv: isobaric specific heat capacity of water vapor. Set to 1870 J/kg/K by default.
        cpl: isobaric specific heat capacity of liquid water. Set to 4186 J/kg/K by default.
        cpi: isobaric specific heat capacity of ice. Set to 2108 J/kg/K by default.
        Rd: gas constant of dry air. Set to 287 J/kg/K by default.
        Rv: gas constant of water vapor. Set to 461 J/kg/K by default.
        g: gravitational acceleration constant. Set to 9.81 m/s**2 by default.
    """

    # if any meters are unspecified, set them to the following defaults.
    if frz is None: frz = 0
    if ma_type is None: ma_type = 'std'
    if dpa is None: dpa = -10.
    if pa_top is None: pa_top = 100e2
    if eps is None: eps = 0.622
    if cpd is None: cpd = 1005.7
    if cpv is None: cpv = 1870.
    if cpl is None: cpl = 4186.
    if cpi is None: cpi = 2108.
    if Rd is None: Rd = 287.
    if Rv is None: Rv = 461.
    if g is None: g = 9.81

    if np.isnan(ps) | np.isnan(tas) | np.isnan(hurs):
        ta_s = nan(plev.shape); # if any of the initial values are nan, output nan profile
    else:
        pa_cd = []
        ta_cd = []
        esat_cd = []
        qsat_cd = []
        rh = []
        r = []
        deriv_dry = []

        pa_cd.append(ps); # set initial pressure
        ta_cd.append(tas); # set initial temperature
        esat_cd.append(calc_esat(ta_cd[0], frz)); # calculate initial saturation vapor pressure
        qsat_cd.append(calc_q(pa_cd[0], esat_cd[0], eps)); # calculate initial specific humidity
        rh.append(hurs/100); # set initial relative humidity

        r.append(eps * rh[0] * esat_cd[0] / (pa_cd[0] - rh[0] * esat_cd[0]));
        deriv_dry.append(dry(pa_cd[0], ta_cd[0], Rd, cpd));
        sol=solve_ivp(dry, (ps, pa_top), ta_cd[0].reshape(1,), args=(Rd, cpd, ));
        pa_d = sol.t; ta_pa_d = sol.y; sol = None;
        i = 1;
        while rh[i-1] < 1:
            idx_lcl = i;
            pa_cd.append(pa_cd[i-1] + dpa);
            ta_cd.append(ta_cd[i-1] + dpa * deriv_dry[i-1]);
            esat_cd.append(calc_esat(ta_cd[i], frz));
            qsat_cd.append(calc_q(pa_cd[i], esat_cd[i], eps));
            r.append(r[i-1]);
            rh.append(r[i] * pa_cd[i] / (esat_cd[i] * (r[i] + eps)));
            deriv_dry.append(dry(pa_cd[i], ta_cd[i], Rd, cpd));
            i = i + 1;

        pa_lcl = pa_cd[-1];
        peval = plev[np.logical_and(plev<pa_lcl, plev>pa_top)]
        peval[::-1].sort() # sort pressure levels to descending order
        if ma_type == 'reversible':
            r_t = r[-2]; # moisture at LCL is the total moisture of rising el
            sol = solve_ivp(sat_rev, [pa_lcl, pa_top], ta_cd[-1].reshape(1,), args=(Rd, cpd, cpv, cpl, cpi, eps, frz, ));
            pa_cs = sol.t; ta_cs = sol.y; sol = None;
        elif ma_type == 'pseudo':
            sol = solve_ivp(sat_pseudo, [pa_lcl, pa_top], ta_cd[-1].reshape(1,), args=(Rd, cpd, cpv, eps, frz, ));
            pa_cs = sol.t; ta_cs = sol.y; sol = None;
        elif ma_type == 'std':
            sol = solve_ivp(sat_std, [pa_lcl, pa_top], ta_cd[-1].reshape(1,), t_eval=peval, args=(Rd, Rv, cpd, cpv, eps, frz, ));
            pa_cs = sol.t; ta_cs = sol.y; sol = None;
        else:
            error(sprintf('The chosen moist adiabat type, #s, is not available. Choose between reversible, pseudo, or std.', ma_type));

        deriv_cs = np.empty_like(pa_cs)
        deriv_csd = np.empty_like(pa_cs)
        qsat_cs = np.empty_like(pa_cs)
        for i in np.arange(0,len(pa_cs)):
            if ma_type == 'reversible':
                deriv_cs[i], qsat_cs[i] = sat_rev(pa_cs[i], ta_cs[0,i], Rd, cpd, cpv, cpl, cpi, eps, frz, diagq=1);
            elif ma_type == 'pseudo':
                deriv_cs[i], qsat_cs[i] = sat_pseudo(pa_cs[i], ta_cs[0,i], Rd, cpd, cpv, eps, frz, diagq=1);
            elif ma_type == 'std':
                deriv_cs[i], qsat_cs[i] = sat_std(pa_cs[i], ta_cs[0,i], Rd, Rv, cpd, cpv, eps, frz, diagq=1);
            else:
                print('The chosen moist adiabat type, %s, is not available. Choose between reversible, pseudo, or std.' % (ma_type))
            deriv_csd[i] = dry(pa_cs[i], ta_cs[0,i], Rd, cpd);

        ta_s = np.append(ta_cd, ta_cs[0,1:-1]);
        pa_s = np.append(pa_cd, pa_cs[1:-1]);
        qsat_s = np.append(qsat_cd, qsat_cs[1:-1]);
        dtadpa_s = np.append(deriv_dry, deriv_cs[1:-1]);

        ta_s_int = splrep(np.flip(pa_s), np.flip(ta_s));
        qsat_s_int = splrep(np.flip(pa_s), np.flip(qsat_s));
        dtadpa_s_int = splrep(np.flip(pa_s), np.flip(dtadpa_s));

        ta_s = splev(plev, ta_s_int, ext=1);
        ta_s[ta_s==0] = np.nan
        qsat_s = splev(plev, qsat_s_int, ext=1);
        qsat_s[qsat_s==0] = np.nan
        dtadpa_s = splev(plev, dtadpa_s_int, ext=1);
        dtadpa_s[dtadpa_s==0] = np.nan

    return ta_s


# auxiliary functions
def calc_esat(T, frz):
    # Emanuel (1994) p.116
    # esat = 100 * np.exp(53.67957 - 6743.769 / T - 4.8451 * np.log(T));

    # Bolton (1980) as in Emanuel (1994) p.117
    #esat = 611.2 * np.exp(17.67 * (T - 273.15) / (T - 273.15 + 243.5));

    # Bohren and Albrecht (1998)
    if frz == 0:
        esat = 611 * np.exp(6808 * (1/273 - 1/T) - 5.09 * np.log(T/273));
    elif frz == 1:
        if T > 273:
            esat = 611 * np.exp(6808 * (1/273 - 1/T) - 5.09 * np.log(T/273));
        else:
            esat = 611 * np.exp(6293 * (1/273 - 1/T) - 0.555 * np.log(T/273));
    elif frz == 1.2:
        if T > 273 - 40:
            esat = 611 * np.exp(6808 * (1/273 - 1/T) - 5.09 * np.log(T/273));
        else:
            esat = 611 * np.exp(6293 * (1/273 - 1/T) - 0.555 * np.log(T/273));

    return esat

def calc_q(p, e, eps):
    q = eps * e / p;

    # from Emanuel
    # q = eps * e / (p - e * (1 - eps));
    return q

def calc_L(T, frz=None):

    """
    Options for "frz"
        0   = no freezing; use Lv for all temperatures
        0.1 = 0, but with fixed Lv
        1   = freezing happens immediately; step function
        1.1 = 1, but with fixed Lv
        1.2 = freezing happens immediately at -40 C
        2   = gradual freezing; ramp function
        2.1 = 2, but with fixed Lv
        3   = gradual freezing; ramp function (wide)
        3.1 = 3, but with fixed Lv
        4   = ramp from 0 to -20 C
        4.1 = 4, but with fixed Lv
        5   = ramp from 0 to -35 C
        5.1 = 5, but with fixed Lv
    """

    # set freezing effect off by default
    if frz is None: frz = 0

    # freezing temperature
    Ti = 273.15;

    # half width of transition (default is delT)
    delT = 1;
    wide_delT = 3;

    # define Ls (sublimation) and Lv (vaporation)
    Ls = 2.835e6;
    Lv = 2.501 * 10**6 - (4190. - 1870.) * (T - 273.15);   # from emanuel (1994), p.115
    Lvf = 2.501 * 10**6 - (4190. - 1870.) * (15.);  # use 15 deg C as fixed Lv

    if frz == 0:
        L = Lv;

    elif frz == 0.1:
        L = Lvf;

    # step
    elif frz == 1:
        if T < Ti:
            L = Ls;
        else:
            L = Lv;

    elif frz == 1.1:
        if T < Ti:
            L = Ls;
        else:
            L = Lvf;

    elif frz == 1.2:
        if T < Ti - 40:
            L = Ls;
        else:
            L = Lv;

    # ramp
    elif frz == 2:
        if T < Ti - delT:
            L = Ls;
        elif T > Ti + delT:
            L = Lv;
        else:
            slope = (Lv - Ls)/(2 * delT);
            L = Ls + slope * (T - (Ti - delT));

    elif frz == 2.1:
        if T < Ti - delT:
            L = Ls;
        elif T > Ti + delT:
            L = Lvf;
        else:
            slope = (Lvf - Ls)/(2 * delT);
            L = Ls + slope * (T - (Ti - delT));

    # wide ramp
    elif frz == 3:
        delT = wide_delT;

        if T < Ti - delT:
            L = Ls;
        elif T > Ti + delT:
            L = Lv;
        else:
            slope = (Lv - Ls)/(2 * delT);
            L = Ls + slope * (T - (Ti - delT));

    elif frz == 3.1:
        delT = wide_delT;

        if T < Ti - delT:
            L = Ls;
        elif T > Ti + delT:
            L = Lvf;
        else:
            slope = (Lvf - Ls)/(2 * delT);
            L = Ls + slope * (T - (Ti - delT));

    # shifted ramp
    elif frz == 4:
        delT = 10;

        if T < Ti - 10 - delT:
            L = Ls;
        elif T > Ti -10 + delT:
            L = Lv;
        else:
            slope = (Lv - Ls)/(2 * delT);
            L = Ls + slope * (T - (Ti - 10 - delT));

    elif frz == 4.1:
        delT = 10;

        if T < Ti - 10 - delT:
            L = Ls;
        elif T > Ti - 10 + delT:
            L = Lvf;
        else:
            slope = (Lvf - Ls)/(2 * delT);
            L = Ls + slope * (T - (Ti - 10 - delT));

    # shifted ramp
    elif frz == 5:
        delT = 35/2;

        if T < Ti - 35/2 - delT:
            L = Ls;
        elif T > Ti - 35/2 + delT:
            L = Lv;
        else:
            slope = (Lv - Ls)/(2 * delT);
            L = Ls + slope * (T - (Ti - 35/2 - delT));

    elif frz == 5.1:
        delT = 35/2;

        if T < Ti - 35/2 - delT:
            L = Ls;
        elif T > Ti - 35/2 + delT:
            L = Lvf;
        else:
            slope = (Lvf - Ls)/(2 * delT);
            L = Ls + slope * (T - (Ti - 35/2 - delT));

    return L

def dry(p, T, Rd, cpd):
    gamma = Rd * T / (p * cpd);
    return gamma

def sat_rev(p, T, Rd, cpd, cpv, cpl, cpi, eps, frz=None, diagq=None):
    if frz is None: frz=0
    if diagq is None: diagq=0

    L = calc_L(T, frz);
    esat = calc_esat(T, frz);
    qsat = calc_q(p, esat, eps);
    rsat = eps * esat / (p - esat);
    gamma_dry = dry(p, T, Rd, cpd);
    r_l = r_t - rsat;
    if frz:
        gamma = gamma_dry * ((1 + r_t) *(1 + (L * rsat) / (Rd * T))) / (1 + rsat * cpv / cpd + r_l * cpi / cpd + (L**2 * rsat * (eps + rsat)) / (Rd * T**2 * cpd));
    else:
        gamma = gamma_dry * ((1 + r_t) *(1 + (L * rsat) / (Rd * T))) / (1 + rsat * cpv / cpd + r_l * cpl / cpd + (L**2 * rsat * (eps + rsat)) / (Rd * T**2 * cpd));
    # gamma = gamma_dry * ((1 + r_t) / (1 + rsat * cpv / cpd)) * ((1 + (L * rsat)/(Rd * T)) / (1 + r_l * cpl / (cpd + rsat * cpv) + (L**2 * rsat * (1 + rsat / eps)) / (Rv * T**2 * (cpd + rsat * cpv))));

    if diagq == 0:
        return gamma
    elif diagq == 1:
        return gamma, qsat

def sat_pseudo(p, T, Rd, cpd, cpv, eps, frz=None, diagq=None):
    if frz is None: frz=0
    if diagq is None: diagq=0

    L = calc_L(T, frz);
    esat = calc_esat(T, frz);
    qsat = calc_q(p, esat, eps);
    rsat = eps * esat / (p - esat);
    gamma_dry = dry(p, T, Rd, cpd);

    gamma = gamma_dry * ((1 + rsat) * (1 + (L * rsat) / (Rd * T))) / (1 + rsat*cpd/cpv + (L**2 * rsat) / (Rd * T**2 * cpd) * (eps + rsat));

    if diagq == 0:
        return gamma
    elif diagq == 1:
        return gamma, qsat

def sat_std(p, T, Rd, Rv, cpd, cpv, eps, frz=None, diagq=None):
    if frz is None: frz=0
    if diagq is None: diagq=0

    L = calc_L(T, frz);
    esat = calc_esat(T, frz);
    qsat = calc_q(p, esat, eps);
    rsat = eps * esat / (p - esat);
    gamma_dry = dry(p, T, Rd, cpd);

    gamma = gamma_dry * (1 + (L * qsat) / (Rd * T)) / (1 + (L**2 * qsat) / (Rv * T**2 * cpd));

    if diagq == 0:
        return gamma
    elif diagq == 1:
        return gamma, qsat
