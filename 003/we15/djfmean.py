import sys
import pickle
from tqdm import tqdm
import numpy as np

f_in = float(sys.argv[1])
plotpath = sys.argv[2]

# load data
raw=pickle.load(open('%s/raw_data.pickle' % (plotpath), 'rb'))

#-------------------------------------------------------------------------
#output only converged, final year
td = raw['td']
ntd = len(td) # total number of time steps
nyr = int(np.ceil(td[-1])+1) # total number of years
tmon = int(ntd/nyr/12) # number of time steps per month

djfmean = {}
djfmean['x'] = raw['x']
djfmean['Lf'] = raw['Lf']
djfmean['varnames'] = raw['varnames']
for varname in raw['varnames']:
    djfmean[varname] = np.zeros([len(raw['x']), nyr])

# first season indices
sel0 = np.concatenate((np.arange(tmon)-tmon, np.arange(tmon), np.arange(tmon)+tmon))

for seas in tqdm(np.arange(nyr)):
    isel = sel0 + seas*tmon*12
    null = np.concatenate((np.where(isel<0)[0], np.where(isel >= ntd)[0])) # discard indices outside of range (i.e. remove dec from first season and jan and feb from last season)
    isel = np.delete(isel, null)

    for varname in raw['varnames']:
        djfmean[varname][:,seas] = np.mean(np.take(raw[varname], isel, axis=1), axis=1)

# save data
pickle.dump(djfmean, open('%s/djfmean.pickle' % (plotpath), 'wb'))

