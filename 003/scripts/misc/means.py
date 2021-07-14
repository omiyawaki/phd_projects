import sys
import numpy as np
from scipy.interpolate import interp1d

def lat_mean(var, grid, lat, **kwargs):

    dim = kwargs.get('dim')

    new_shape = np.empty(len(var.shape))
    for i in range(len(var.shape)):
        if i == dim:
            new_shape[i] = lat.size
        else:
            new_shape[i] = var.shape[i]
    new_shape = tuple(new_shape.astype(int))

    var_latmean = np.empty(new_shape)

    clat = np.cos(np.radians(lat))
    
    var_fint = interp1d(grid['lat'], var, axis=dim, bounds_error=False)
    var_int = var_fint(lat)
    nanfilt = np.copy(var_int)
    nanfilt[~np.isnan(nanfilt)] = 1

    var_latmean = np.nansum(var_int*clat,dim)/np.nansum(nanfilt*clat,dim)

    return var_latmean 

def global_int(var, grid, a, **kwargs):

    dim_lat = kwargs.get('dim_lat')
    dim_lon = kwargs.get('dim_lon')

    var_lonint = np.trapz(var, np.radians(grid['lon']), axis=dim_lon)

    clat = np.cos(np.radians(grid['lat']))

    var_gint = a**2 * np.trapz(var_lonint*clat, np.radians(grid['lat']), axis=dim_lat)

    return var_gint

def mmm_mon_lat(var, grid_in, grid_out, **kwargs):
    
    var_list = np.empty([len(var), var[0].shape[0], len(grid_out['lat'])])
    for i in range(len(var)):
        var_fint_lat = interp1d(grid_in[i]['lat'], var[i], axis=1, bounds_error=False)
        var_list[i,:,:] = var_fint_lat(grid_out['lat'])
    
    var_out = {}    
    var_out['mmm'] = np.squeeze(np.mean(var_list, axis=0))
    var_out['prc25'] = np.squeeze(np.percentile(var_list, 25, axis=0))
    var_out['prc75'] = np.squeeze(np.percentile(var_list, 75, axis=0))
    var_out['min'] = np.squeeze(np.min(var_list, axis=0))
    var_out['max'] = np.squeeze(np.max(var_list, axis=0))
    
    return var_out