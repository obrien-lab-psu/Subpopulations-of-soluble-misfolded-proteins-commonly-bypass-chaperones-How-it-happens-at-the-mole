#!/usr/bin/env python3
import numpy as np
import sys,os
import glob
from scipy import stats

if len(sys.argv) != 3:
    print('[1] path to GQ files')
    print('[2] comma separated list of tags to find')
    quit()

xstal_files = glob.glob(f'{sys.argv[1]}*_xstal_*GQ.txt')
print(xstal_files)

labels = sys.argv[2].split(',')
print(labels)

native_data = []
for f in xstal_files:
    native_data.append(np.loadtxt(f)[:,-1:])

native_data = np.vstack(native_data).flatten()
print(native_data, native_data.shape)
native_mode = stats.mode(native_data).mode
native_std = np.std(native_data)
print(native_mode, native_std)

outdata = {}
for label in labels:
    outdata[label] = 0
    for t in range(20):
        f = f'entanglement_analysis_2.0/output/p3hwo_{label}_t{t+1}_GQ.txt'
        print(t, label, f)

        qdata = np.loadtxt(f)[-100:,-1:]
        qdata_mode = stats.mode(qdata).mode
        gdata = np.loadtxt(f)[-100:,-2]
        gdata_mode = stats.mode(gdata).mode

        if qdata_mode > native_mode-3*native_std:
            folded = True
            outdata[label] += 1

        else:
            folded = False
        print(qdata_mode, gdata_mode, folded)

for label,frac in outdata.items():
    print(f'{label}: {frac/20.0}')
