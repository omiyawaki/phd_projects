import sys
from tqdm import tqdm
import numpy as np
from scipy import interpolate

# Converts data to various vertical grids

def pa_to_sigma(varname, var_pa, varname_sfc, var_sfc, var_ps, si):
    # converts pressure grid data to sigma coordinates
    # varname is the name of the variable to be converted (e.g. ta for temperature)
    # var_pa is the original data to be converted in pressure coordinates
    # varname_sfc is the name of the surface companion of the variable to be converted
    # varname_sfc is the surface companion of the data to be converted
    # si is the sigma grid that the data should be converted to

    # rename input data for easier access
    vpa = var_pa[varname]	
    vsfc = var_sfc[varname_sfc]
    ps = var_ps['ps']

    lat2d = var_sfc['grid']['lat']
    lat3d = var_pa['grid']['lat']
    pa = var_pa['grid']['lev']

    # initialize output array in sigma coord
    vsi = np.empty([vpa.shape[0], si.size, vpa.shape[2], vpa.shape[3]])
    var_si = {}
    var_si['grid'] = {}
    var_si['grid']['lat'] = var_pa['grid']['lat']
    var_si['grid']['lon'] = var_pa['grid']['lon']
    var_si['grid']['lev'] = si

    var_pa = None; var_sfc = None; var_ps = None;

    # create surface mask
    vpa = vpa.filled(fill_value=np.nan)

    # check if 2d and 3d lat grids are the same and interpolate 2d data to 3d grid if different
    if not np.array_equal(lat2d, lat3d):
        ps_filled = ps.filled(fill_value=np.nan)
        lat2d_filled = lat2d.filled(fill_value=np.nan)
        lat3d_filled = lat3d.filled(fill_value=np.nan)
        f = interpolate.interp1d(lat2d_filled, ps_filled, axis=1)
        ps = f(lat3d_filled)
        ps_filled = None; f = None;

    # replace subsurface data with nans
    ps3d = np.tile(ps, [pa.size, 1, 1, 1])
    ps3d = np.transpose(ps3d, [1, 0, 2, 3])
    pa3d = np.tile(pa, [ps.shape[0], ps.shape[1], ps.shape[2], 1])
    pa3d = np.transpose(pa3d, [0, 3, 1, 2])
    idx_subsrfc = pa3d > 0.9961*ps3d
    ps3d = None; pa3d = None;

    vpa[idx_subsrfc] = np.nan

    pa_dec = 1 if pa[2]-pa[1] < 0 else 0 # check if pa is increasing or decreasing

    for ilon in tqdm(range(vpa.shape[3])):
        for ilat in range(vpa.shape[2]):
            for itime in range(vpa.shape[0]):
                idx_atm_local = ~idx_subsrfc[itime,:,ilat,ilon]
                si_local = (pa/ps[itime,ilat,ilon])[idx_atm_local]
                vpa_local = (vpa[itime,:,ilat,ilon])[idx_atm_local]

                if len(vsfc.shape)==2: # if surface data is not a function of time, only read using lat and lon indices (e.g. for surface orography)
                    # vsfc_local = vsfc[ilat,ilon]
                    
                    # UPDATE: instead of using orography, extrapolate surface geopotential
                    vsfc_int = interpolate.interp1d(pa[idx_atm_local], vpa_local, kind='linear', fill_value='extrapolate')
                    vsfc_local = vsfc_int(ps[itime,ilat,ilon])
                else:
                    vsfc_local = vsfc[itime,ilat,ilon]

                if pa_dec: # flip data if pressure is decreasing (required for interpolation)
                    si_local = np.flip(si_local)
                    vpa_local = np.flip(vpa_local)

                # add surface value
                si_local = np.append(si_local, 1)
                vpa_local = np.append(vpa_local, vsfc_local)

                # f = interpolate.interp1d(si_local, vpa_local, kind='cubic')
                f = interpolate.interp1d(si_local, vpa_local, kind='linear')
                vsi[itime,:,ilat,ilon] = f(si)

                # tck = interpolate.splrep(si_local, vpa_local)
                # vsi[itime,:,ilat,ilon] = interpolate.splev(si, tck)

    var_si[varname] = vsi

    return var_si
