import os

def get_datadir(sim, **kwargs):

    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')

    datadir = '/project2/tas1/miyawaki/projects/003/data/proc/%s' % (suffix_str(sim, model, yr_span))

    if not os.path.isdir(datadir):
        os.makedirs(datadir)

    return datadir

def get_plotdir(sim, **kwargs):

    model = kwargs.get('model')
    yr_span = kwargs.get('yr_span')
    categ = kwargs.get('categ')

    plotdir = '/project2/tas1/miyawaki/projects/003/plot/%s/%s' % (suffix_str(sim, model, yr_span), categ)

    if not os.path.isdir(plotdir):
        os.makedirs(plotdir)

    return plotdir

def suffix_str(sim, model, yr_span, **kwargs):

    if sim in ['hist+ssp585', 'ssp585', 'hist+rcp85', 'rcp85', 'historical', 'echam']:
       suffix_str = '%s/%s/%s' % (sim, model, yr_span)
    elif sim == 'longrun':
       suffix_str = '%s/%s' % (sim, model)
    elif sim == 'era5':
       suffix_str = '%s/%s' % (sim, yr_span)

    return suffix_str
