import sys
import glob
from netCDF4 import Dataset
from misc.translate import translate_varcat

def filenames_raw(sim, varname, **kwargs):

    model = kwargs.get('model') # name of model
    timemean = kwargs.get('timemean', 'yearmean') # do annual mean? (bool)
    yr_span = kwargs.get('yr_span') # considered span of years
    
    if sim in ['hist+ssp585', 'ssp585', 'hist+rcp85', 'rcp85', 'historical']:
        if varname in ['orog', 'sftlf']:
            fname = glob.glob(remove_repdots('/project2/tas1/miyawaki/projects/003/data/raw/%s/%s/%s_*_%s_%s_*.nc' % (sim, model, varname, model, sim)))
        else:
            try:
                if varname in ['sic']:
                    fname = glob.glob(remove_repdots('/project2/tas1/miyawaki/projects/003/data/raw/%s/%s/%s_OImon_%s_%s_*_%s.%s.nc' % (sim, model, varname, model, sim, yr_span, timemean)))
                elif varname in ['siconc']:
                    fname = glob.glob(remove_repdots('/project2/tas1/miyawaki/projects/003/data/raw/%s/%s/%s_SImon_%s_%s_*_%s.%s.nc' % (sim, model, varname, model, sim, yr_span, timemean)))
                elif varname in ['diffv92500', 'tdiffv92500', 'vmte', 'aht', 'gmse92500', 'mse92500']:
                    fname = glob.glob(remove_repdots('/project2/tas1/miyawaki/projects/003/data/raw/%s/%s/%s_Amon_%s_%s_*_%s.zonmean.shsmooth.%s.nc' % (sim, model, varname, model, sim, yr_span, timemean)))
                else:
                    fname = glob.glob(remove_repdots('/project2/tas1/miyawaki/projects/003/data/raw/%s/%s/%s_Amon_%s_%s_*_%s.%s.nc' % (sim, model, varname, model, sim, yr_span, timemean)))
            except:
                print('/project2/tas1/miyawaki/projects/003/data/raw/%s/%s/%s_Amon_%s_*_%s.%s.nc' % (sim, model, varname, model, sim, yr_span, timemean))
        try:
            filename = Dataset(fname[0], 'r')
        except:
            print('\nERROR: The following file cannot be read: /project2/tas1/miyawaki/projects/003/data/raw/%s/%s/%s_Amon_%s_%s_*_%s.%s.nc\n' % (sim, model, varname, model, sim, yr_span, timemean))

    elif sim == 'longrun':
        filename = Dataset(remove_repdots('/project2/tas1/miyawaki/projects/003/data/raw/%s/%s/%s_mon_%s_%s.%s.nc' % (sim, model, varname, model, yr_span, timemean)), 'r')

    elif sim == 'echam':
        print('/project2/tas1/miyawaki/projects/003/data/raw/%s/%s/%s_%s_%s.%s.nc' % (sim, model, varname, model, yr_span, timemean))
        fname = glob.glob(remove_repdots('/project2/tas1/miyawaki/projects/003/data/raw/%s/%s/%s_%s_%s.%s.nc' % (sim, model, varname, model, yr_span, timemean)))
        print(fname)
        filename = Dataset(fname[0], 'r')

    elif sim == 'era5':
        filename = Dataset(remove_repdots('/project2/tas1/miyawaki/projects/002/data/raw/era5c/%s/era5c_%s_%s.%s.nc' % (translate_varcat(varname), translate_varcat(varname), yr_span, timemean)), 'r')

    return filename

def remove_repdots(filename_in):
    # removes repeating periods from filename

    import re

    remove_str = '.'
    pattern = "(?P<char>[" + re.escape(remove_str) + "])(?P=char)+"
    filename_out = re.sub(pattern, r"\1", filename_in)

    return filename_out
