#!/usr/bin/env python3
import sys,os
import numpy as np
from scipy.stats import bootstrap

#get user arguement
if len(sys.argv) != 3:
    print('[1] = input directory containing analysis files from binding_v1.1.py')
    print('[2] = output file')
    quit()

#get files
data = []
for file in os.listdir(sys.argv[1]):
        if file.endswith("_binding.txt"):

            file = os.path.join(sys.argv[1], file)

            file_data = np.loadtxt(file, usecols=(-1))
            data.append(file_data)

#reformat data to 1D array
data = np.stack(data).flatten()[:,None].T

#bootstrap
mean = np.mean(data)
res = bootstrap(data, np.mean, confidence_level=0.95, method='percentile', n_resamples=10000)
outdata = np.asarray([mean, res.confidence_interval[0], res.confidence_interval[1]])[:,None].T

#save outdata
np.savetxt(sys.argv[2], outdata, header='mean lower_95%ci upper_95%ci', fmt='%f1')
print(f'SAVED: {sys.argv[2]}')
