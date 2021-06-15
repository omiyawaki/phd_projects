import glob
from netCDF4 import Dataset
from misc.translate import translate_varcat

def filenames_raw(sim, varname, **kwargs):

    model = kwargs.get('model') # name of model
    timemean = kwargs.get('timemean', '.yearmean') # do annual mean? (bool)
    yr_span = kwargs.get('yr_span') # considered span of years

    if sim == 'rcp85':
        fname = glob.glob('/project2/tas1/miyawaki/projects/003/data/raw/%s/%s/%s_*_%s_%s_r1i1p1_%s%s.nc' % (sim, model, varname, model, sim, yr_span, timemean))
        filename = Dataset(fname[0], 'r')
    elif sim == 'longrun':
        filename = Dataset('/project2/tas1/miyawaki/projects/003/data/raw/%s/%s/%s_mon_%s_%s%s.nc' % (sim, model, varname, model, yr_span, timemean), 'r')
    elif sim == 'echam':
        filename = Dataset('/project2/tas1/miyawaki/projects/003/data/raw/%s/%s/%s_%s_%s%s.nc' % (sim, model, varname, model, yr_span, timemean), 'r')
    elif sim == 'era5':
        filename = Dataset('/project2/tas1/miyawaki/projects/002/data/raw/era5c/%s/era5c_%s_%s%s.nc' % (translate_varcat(varname), translate_varcat(varname), yr_span, timemean), 'r')

    return filename