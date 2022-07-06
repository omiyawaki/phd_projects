import sys
import pickle
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (4, 3)

f_in = float(sys.argv[1])
plotpath = sys.argv[2]

# load data
raw=pickle.load(open('%s/raw_data.pickle' % (plotpath), 'rb'))

td = raw['td']
nyr = int(np.ceil(td[-1]))

yearmean = {}
yearmean['x'] = raw['x']
yearmean['varnames'] = raw['varnames']
for varname in raw['varnames']:
    yearmean[varname] = np.zeros([len(raw['x']), nyr])

for yr in tqdm(np.arange(nyr)):
    iyr = np.where(td >=yr)  and np.where(td < yr+1)

    for varname in raw['varnames']:
        yearmean[varname][:,yr] = np.mean(np.take(raw[varname], iyr[0], axis=1), axis=1)

# save data
pickle.dump(yearmean, open('%s/yearmean.pickle' % (plotpath), 'wb'))
