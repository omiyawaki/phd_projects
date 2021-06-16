import sys
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

	vpa = var_pa[varname]	
	vsfc = var_sfc[varname_sfc]
	ps = var_ps['ps']

	lat2d = var_sfc['grid']['lat']
	lat3d = var_pa['grid']['lat']
	pa = var_pa['grid']['lev']

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
	print(pa3d.shape)

	sys.exit()



	var_si = None

	return var_si